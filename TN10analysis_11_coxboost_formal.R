library(CoxBoost) # ver 1.4 
library(mboost) # ver 2.9-5
library(tidyverse)
library(magrittr)
library(caret)

source('modified_model_based_functions.R')
source('modified_likelihood_based_functions.R')
nCV<-10 # used 10 times repeated 10 fold CV for caret so also using it here 

tn10<-read.csv("TN10full_smallgate.csv", header = T, as.is = T)
tn10$Date_of_Draw<-as.Date(tn10$Date_of_Draw)
tn10$Date.IDDM<-as.Date(tn10$Date.IDDM)

str(tn10)
theme_set(
  theme_bw(base_size = 15)
)

tn10%<>%fill(T1D.event,Time.to.T1D,Age,Gender,Treatment.Arm)%>% # for 3 samples w/o ACAb info but has Cpep and other metabolic measures
  mutate(Sex = factor(Gender))%>%select(-Gender)%>%filter(grepl("SENT",DR4))%>%
  mutate(Date_of_Draw = as.Date(Date_of_Draw),
         Date.IDDM = as.Date(Date.IDDM),
         Sex = ifelse(Sex == "Male", -1, 1),
         DR3 = ifelse(DR3 == "ABSENT", -1, 1),
         DR4 = ifelse(DR4 == "ABSENT", -1, 1),
         Treatment.Arm = ifelse(Treatment.Arm =="Placebo", -1,1),
         GADA = ifelse(GADA == 0, -1, 1),
         IA.2A = ifelse(IA.2A == 0, -1, 1),
         mIAA = ifelse(mIAA == 0, -1, 1),
         ZnT8A = ifelse(ZnT8A == 0, -1, 1),
         glycemia = factor(glycemia, levels = c("normal","dysglycemia","hyperglycemia")))

# didn't factor Visit here because need to merge with 3 month data for time-dependent cox
tn10[is.na(tn10$Barcode),]
tn10[tn10$MaskID == "245390" & tn10$Visit == 12,"Date.IDDM"]<-tn10[tn10$MaskID == "245390" & tn10$Visit == 0,"Date.IDDM"]
# check that all individuals whose T1D.event == 0 have na Date.IDDM
tn10%>%count(Date.IDDM,T1D.event)%>%spread(T1D.event,n)

tn10_ppc<-preProcess(tn10%>%select(-MaskID,-Visit,-ID,-T1D.event,-Time.to.T1D,-Sex,-matches("A$",ignore.case = F)),method = c("center", "scale","YeoJohnson"))
tn10.scaled<- predict(tn10_ppc, tn10)

# functions to implement a repeated cross-validation ----------------------
# for deciding the number of boosting iterations m_stop
# 15 repeated ten-fold CV
rep_cv.mb<-function(seed,mod)
{
  set.seed(seed)
  apply(cvrisk(mod,folds=cv(model.weights(mod),type='kfold',B=10),papply = lapply),2,mean) # from mboost, offset
}
rep_cv.mbN<-function(seed,y,xx,mandatory,maxstep)
{
  set.seed(seed)
  cvpl.l2boostN(maxstep,y=y,xx=xx,mandatory=mandatory)$cvllik # from modified_model_based_functions.R
}
rep_cv.lb<-function(seed,maxstep,time,status,xx,penalty,unpen.index=NULL)
{
  set.seed(seed)
  cv.CoxBoost(time,status,xx,penalty=penalty,maxstepno=maxstep,unpen.index=unpen.index)$mean.logplik # from CoxBoost, fovor strategy
}
rep_cv.lbN<-function(seed,maxstep,y,xx,offset,lambda){
  cvpl.offBoost(maxstep,y,xx,offset=offset,lambda=lambda,foldCV=10,seed=seed)$cvllik # modified_likelihood_based_functions.R
}

base<-tn10%>%filter(Visit == 0)%>%
  select(MaskID,Time.to.T1D,T1D.event,Age,Sex,BMI,DR3,DR4,Treatment.Arm,mu_cpep,GADA,IA.2A,mIAA,ZnT8A,contains(c("IgG1","IgG2","IgA")))%>%
  filter(!is.na(Bd.IgG2.RI) )%>%distinct(MaskID, .keep_all = T)
rownames(base)<-paste0("X",base$MaskID)
set.seed(42)
inTrain <- createDataPartition( ## caret package for splitting train-test set
  y = base$T1D.event, ## the outcome data are needed
  p = .8, ## The percentage of data in the training set
  list = FALSE
)
base.bin<-base%>%select(Sex,DR3,DR4,Treatment.Arm,GADA,IA.2A,mIAA,ZnT8A)
base.bin["X787396","IA.2A"]<--1
base.cont<-base%>%select(Age,BMI,mu_cpep,contains(c("IgG1","IgG2","IgA")))


outcome<-Surv(base$Time.to.T1D,base$T1D.event)

models.baseline<-ctime.base<-NULL

lambda<-(length(inTrain)-1)*9

# ------------------- unscaled,baseline ----------------------
### offset strategy ###
# split base table into clinical (madatory) features and optional features
base.clin<-cbind(select(base.bin,Sex,DR3,DR4,Treatment.Arm),select(base,Age,BMI,mu_cpep))
base.acab<-cbind(base.bin[,grep("A$",colnames(base.bin))],select(base,contains(c("IgG1","IgG2","IgA"))))
# compute the offset
mod.clin<-coxph(outcome[inTrain,]~.,data=base.clin[inTrain,])
tmp<-mod.clin$linear.predictors
clin.predictor<-rep(0,length(inTrain))
clin.predictor[inTrain]<-tmp

# model-based boosting 
mod1<-glmboost(y=outcome[inTrain,],x=as.matrix(base.acab)[inTrain,],family=CoxPH(),offset=tmp)
par1<-which.min(apply(sapply(1:nCV,rep_cv.mb,mod=mod1),1,mean))
data1<-data.frame(base$Time.to.T1D,base$T1D.event,clin.predictor,base.acab[,names(coefficients(mod1[par1]))])[inTrain,]
if(ncol(data1)==4) colnames(data1)[4]<-names(coefficients(mod1[par1]))
models.baseline$mb.offset<-coxph(Surv(base.Time.to.T1D,base.T1D.event)~.,data=data1,init = c(1,coefficients(mod1[par1])),iter.max=0,x=T)

### favoring strategy ###
p0<-ncol(base.clin)
mod2<-CoxBoost(base$Time.to.T1D[inTrain],base$T1D.event[inTrain],as.matrix(cbind(base.clin,base.acab)[inTrain,]),penalty=lambda,stepno=par2,unpen.index=1:p0)
par2<-which.min(apply(sapply(1:nCV,rep_cv.mb,mod=mod2),1,mean))
data2<-data.frame(base$Time.to.T1D,base$T1D.event,cbind(base.clin,base.acab)[,mod2$coefficients[par1+1,]!=0])[inTrain,]
if(ncol(data2)==3) colnames(data2)[3]<-colnames(cbind(clin,gene))[mod2$coefficients[par2+1,]!=0]

# ------------------- scaled, baseline ----------------------
base.scaled<- tn10.scaled%>%filter(Visit == 0)%>%distinct(MaskID,.keep_all = T)%>%column_to_rownames(var = "MaskID")%>%select(colnames(base.cont))
### offset strategy ###
# split base table into clinical (madatory) features and optional features
base.clin<-cbind(select(base.bin,Sex,DR3,DR4,Treatment.Arm),select(base.scaled,Age,BMI,mu_cpep))
base.acab<-cbind(base.bin[,grep("A$",colnames(base.bin))],select(base.scaled,contains(c("IgG1","IgG2","IgA"))))
# compute the offset
mod.clin<-coxph(outcome[inTrain,]~.,data=base.clin[inTrain,])
tmp<-mod.clin$linear.predictors
clin.predictor<-rep(0,length(inTrain))
clin.predictor[inTrain]<-tmp

# model-based boosting 
mod3<-glmboost(y=outcome[inTrain,],x=as.matrix(base.acab)[inTrain,],family=CoxPH(),offset=tmp)
par3<-which.min(apply(sapply(1:nCV,rep_cv.mb,mod=mod3),1,mean))
data3<-data.frame(base$Time.to.T1D,base$T1D.event,clin.predictor,base.acab[,names(coefficients(mod3[par3]))])[inTrain,]
if(ncol(data3)==4) colnames(data3)[4]<-names(coefficients(mod3[par3]))
models.baseline$mb.offset<-coxph(Surv(base.Time.to.T1D,base.T1D.event)~.,data=data3,init = c(1,coefficients(mod3[par3])),iter.max=0,x=T)

### favoring strategy ###
par4<-which.max(apply(sapply(1:nCV,rep_cv.lbN,maxstep=100,y=outcome[inTrain,],xx=as.matrix(base.acab[inTrain,]),offset=tmp,lambda=lambda),1,mean))-1
mod4<-CoxBoost(base$Time.to.T1D[inTrain],base$T1D.event[inTrain],as.matrix(cbind(base.clin,base.acab)[inTrain,]),penalty=lambda,stepno=par4,unpen.index=1:p0)
data4<-data.frame(base$Time.to.T1D,base$T1D.event,cbind(base.clin,base.acab)[,mod4$coefficients[par4+1,]!=0])[inTrain,]
if(ncol(data4)==3) colnames(data4)[3]<-colnames(cbind(clin,gene))[mod4$coefficients[par4+1,]!=0]

# ------------------- scaled, 6 month ----------------------
mon6.scaled<- tn10.scaled%>%filter(Visit == 6)%>%select(MaskID,colnames(base.cont),-BMI,-mu_cpep)
oldnames<-colnames(mon6.scaled)[grepl("RI",colnames(mon6.scaled))]
newnames<-paste(oldnames,"T6",sep=".")
mon6.scaled%<>%rename_at(vars(oldnames), ~ newnames)

base.scaled$MaskID<-as.integer(rownames(base.scaled))
all.scaled<-full_join(base.scaled,mon6.scaled)%>%column_to_rownames("MaskID")%>%select(contains("RI"))
### offset strategy ###
mod5<-glmboost(y=outcome[inTrain,],x=as.matrix(all.scaled)[inTrain,],family=CoxPH(),offset=tmp)
par5<-which.min(apply(sapply(1:nCV,rep_cv.mb,mod=mod5),1,mean))
data5<-data.frame(base$Time.to.T1D,base$T1D.event,clin.predictor,all.scaled[,names(coefficients(mod5[par5]))])[inTrain,]

### favoring strategy ###
inTrain2 <- inTrain[!inTrain %in% c(46,50)]
par6<-which.max(apply(sapply(1:nCV,rep_cv.lbN,maxstep=100,y=outcome[inTrain2,],xx=as.matrix(all.scaled[inTrain2,]),offset=tmp[-c(36,40)],lambda=lambda),1,mean))-1
mod6<-CoxBoost(base$Time.to.T1D[inTrain2],base$T1D.event[inTrain2],as.matrix(cbind(base.clin,all.scaled)[inTrain2,]),penalty=lambda,stepno=par6,unpen.index=1:p0)
data6<-data.frame(base$Time.to.T1D,base$T1D.event,cbind(base.clin,all.scaled)[,mod6$coefficients[par4+1,]!=0])[inTrain,]



