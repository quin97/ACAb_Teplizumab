library(survival)
library(tidyverse)
library(magrittr)
library(caret)
library(survminer)
library(ghibli)
library(broom)
library(ggrepel)
library(readxl)
theme_set(
  theme_classic(base_size = 30)
)

diblRaw<-read_xlsx("../TN10RIdata/Data file S1.xlsx",sheet = "RI_validation")
tn01meta.full<-read.csv("20230329_TN01_Sample_Metadata_vialCount.csv",as.is = T)


tn01meta<-tn01meta.full%>%select(Sex,DR3,DR4,Sample_Age,time_to_T1D,Status,flow_id)%>%right_join(diblRaw)
tn01meta%<>%mutate(
  Time.to.T1D = time_to_T1D*365,
  Status = 1,
  DR3 = factor(DR3, levels = c("ABSENT","PRESENT")),
  DR4 = factor(DR4, levels = c("ABSENT","PRESENT")),
  Sex = ifelse(Sex == "Male", -1, 1))

bl.igg2.bu<-tn01meta%>%
  coxph(Surv(Time.to.T1D,Status) ~Bl.IgG2.RI*DR4+Sample_Age+Sex+DR3, data = .)
bl.igg2.bu
cox.zph(bl.igg2.bu) 
ggplot(tn01meta, aes(x=Bl.IgG2.RI))+geom_histogram(bins = 13, fill ="#92BBD9FF",color="#4D6D93FF")+labs(x=expression(italic("B.longum")*"-IgG2"), y = "Sample count")+
  geom_vline(xintercept = 65,linetype = "dashed")# +scale_x_continuous(limits = c(0,200))
summary(tn01meta$Bl.IgG2.RI)


tn01.bligg2.b<-tn01meta%>%mutate(bligg2=ifelse(Bl.IgG2.RI >65, "high","low"))
table(tn01.bligg2.b$bligg2,tn01.bligg2.b$DR4)



(bligg2.bu.survfitsplots<-tn01.bligg2.b%>%split(.$bligg2)%>%imap(~survfit(Surv(Time.to.T1D,Status) ~ DR4, data = .))%>%
    imap(~ggsurvplot(fit=., risk.table = F,conf.int = T,censor.shape = NA,pval = F, break.time.by = 200,
                     data=tn01.bligg2.b%>%filter(bligg2 == .y),ggtheme = theme_classic(base_size = 30),
                     legend = "right",palette = c("#E7B800", "#2E9FDF"),legend.title = bquote(italic("B.longum")*"-IgG2"~ .(.y)),
                     xlab = "On-Study (Days)", ylab = "Proportion T1D-Free")))

di.iga.bu<-tn01meta%>%
  coxph(Surv(Time.to.T1D,Status) ~Di.IgA.RI*DR4+Sample_Age+strata(Sex)+DR3, data = .)
cox.zph(di.iga.bu)
ggplot(tn01meta, aes(x=Di.IgA.RI,y=Di.IgG2.RI))+geom_point()+labs(x=expression(italic("D.invisus")*"-IgA"), y = expression(italic("D.invisus")*"-IgG2"))+geom_smooth(method = "lm")+scale_x_continuous(limits = c(0,100))
lm(Di.IgG2.RI~Di.IgA.RI,data = tn01meta%>%filter(Di.IgA.RI< 200)
)%>%summary()
ggplot(tn01meta, aes(x=Di.IgA.RI))+geom_histogram(bins = 22, fill ="#92BBD9FF",color="#4D6D93FF")+geom_vline(xintercept  = 30,linetype = "dashed")+labs(x=expression(italic("D. invisus")*"-IgA"), y = "Sample count")+scale_x_continuous(limits = c(0,100))


tn01.diiga.b<-tn01meta%>%mutate(diiga=ifelse(Di.IgA.RI > 30, "high","low"))
table(tn01.diiga.b$diiga,tn01meta$DR4)
(diiga.bs.survfitsplots<-tn01.diiga.b%>%split(.$diiga)%>%imap(~survfit(Surv(Time.to.T1D,Status) ~ DR4, data = .))%>%
    imap(~ggsurvplot(fit=., risk.table = F,conf.int = T,censor.shape = NA,pval = F, break.time.by = 200,
                     data=tn01.diiga.b%>%filter(diiga == .y),ggtheme = theme_classic(base_size = 30),
                     legend = "right", palette = c("#E7B800", "#2E9FDF"),legend.title = bquote(italic("D.invisus")*"-IgA"~ .(.y)),
                     xlab = "On-Study (Days)", ylab = "Proportion T1D-Free")))




# pool data ---------------------------------------------------------------

tn10<-read.csv("../TN10analysis_20210827/TN10full_smallgate.csv", header = T, as.is = T)
tn10%<>%fill(T1D.event,Time.to.T1D,Age,Gender,Treatment.Arm)%>% # for 3 samples w/o ACAb info but has Cpep and other metabolic measures
  mutate(Sex = factor(Gender))%>%select(-Gender)%>%
  mutate(Date_of_Draw = as.Date(Date_of_Draw),
         Date.IDDM = as.Date(Date.IDDM),
         Time.landmark = Time.to.T1D-45,
         DR3 = factor(DR3, levels = c("ABSENT","PRESENT")),
         DR4 = factor(DR4, levels = c("ABSENT","PRESENT")),
         Sex = ifelse(Sex == "Male", -1, 1),
         Treatment.Arm = factor(Treatment.Arm), 
         GADA = ifelse(GADA == 0, -1, 1),
         IA.2A = ifelse(IA.2A == 0, -1, 1),
         mIAA = ifelse(mIAA == 0, -1, 1),
         ZnT8A = ifelse(ZnT8A == 0, -1, 1),
         glycemia = factor(glycemia, levels = c("normal","dysglycemia","hyperglycemia")))

tn10[is.na(tn10$Barcode),]
tn10[tn10$MaskID == "245390" & tn10$Visit == 12,"Date.IDDM"]<-tn10[tn10$MaskID == "245390" & tn10$Visit == 0,"Date.IDDM"]
# check that all individuals whose T1D.event == 0 have na Date.IDDM
tn10%>%count(Date.IDDM,T1D.event)%>%spread(T1D.event,n)



tn10.subset<-tn10%>%filter(Visit == 0)%>%distinct(MaskID, .keep_all = T)%>%filter(!is.na(DR3))%>%
  select(MaskID,Sex,DR3,DR4,Age,Time.to.T1D,T1D.event, Treatment.Arm,contains(c("Bl.","Di")))
tn10.subset$cohort = "TN10"
tn10.subset%>%group_by(DR3,DR4)%>%count()

tn01meta$Treatment.Arm<-factor("Placebo",levels = c("Placebo","Teplizumab"))
tn01.subset<-tn01meta%>%select(MaskID = flow_id,Sex,DR3,DR4,Age = Sample_Age,Time.to.T1D, #= time_to_T1D,
                               T1D.event= Status,
                               Treatment.Arm,contains(c("Bl.","Di")))
tn01.subset$cohort = "TN01"
tn.scaled<-rbind(tn10.subset,tn01.subset)
ggplot(tn.scaled, aes(x=Time.to.T1D))+geom_histogram(position="identity",bins = 17,alpha = .75,aes(fill = grepl("_",MaskID)))+labs(x="Time to T1D", y = "Sample count")+
  scale_fill_manual(name = "Cohort",labels = c("TN10","TN01"),values = c("blue","violetred4"))
IQR(tn01.subset$Time.to.T1D)
IQR(tn10.subset$Time.to.T1D)
wilcox.test(Time.to.T1D~cohort,data = tn.scaled)
ggplot(tn.scaled, aes(x=Age))+geom_histogram(position="identity",bins = 17,alpha = .75,aes(fill = grepl("_",MaskID)))+labs(x="Age", y = "Sample count")+scale_fill_manual(name = "Cohort",labels = c("TN10","TN01"),values = c("blue","violetred4"))
mean(tn01.subset$Age)
mean(tn10.subset$Age)
t.test(Age~cohort,data = tn.scaled)


