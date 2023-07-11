library(readxl)
library(tidyverse)
library(magrittr)
library(caret)
library(ggrepel)
library(ghibli)

theme_set(
  theme_classic(base_size = 18)
)

tn10<-read.csv("TN10full_smallgate.csv", header = T, as.is = T)
str(tn10)

tn10%<>%fill(T1D.event,Time.to.T1D,Age,Gender,Treatment.Arm)%>% # for 3 samples w/o ACAb info but has Cpep and other metabolic measures
  mutate(Sex = factor(Gender))%>%select(-Gender,-GADA,-IA.2A,-mIAA,-ZnT8A)%>%
  mutate(Date_of_Draw = as.Date(Date_of_Draw),
         Date.IDDM = as.Date(Date.IDDM),
         Time.landmark = Time.to.T1D-45,
         DR3 = factor(DR3, levels = c("ABSENT","PRESENT")),
         DR4 = factor(DR4, levels = c("ABSENT","PRESENT")),
         Sex = ifelse(Sex == "Male", -1, 1),
         Treatment.Arm = factor(Treatment.Arm, levels=c("Teplizumab", "Placebo")), 
         glycemia = factor(glycemia, levels = c("normal","dysglycemia","hyperglycemia")),
         Visit = factor(Visit))%>%
  filter(!is.na(Barcode))

# check that all individuals whose T1D.event == 0 have na Date.IDDM
tn10%>%count(Date.IDDM,T1D.event)%>%spread(T1D.event,n)

# IAB overtime ------------------------------------------------------------

iabtable<-read_excel("TN10_autoabdata_K_Herold_08_2021.xlsx", skip=1)

iabtable<-iabtable%>%select(MaskID:TGA)%>%filter(!MaskID %in% c(839473,983719))%>%arrange(MaskID, `Draw Date`) # 839473,983719 dont have baseline value
str(iabtable)
iabtable$`Draw Date`<-as.Date(iabtable$`Draw Date`)
iabtable[,5:9]<-apply(iabtable[,5:9], 2, as.numeric)
apply(iabtable[,5:9], 2, summary)
iabtable%>%group_by(Visit)%>%count()
iabtable%<>%
  mutate(Visit = factor(Visit,levels = c("Baseline","12 months","24 months", "36 months")),
         GADA.bin = ifelse(GADA > 20, 1, 0), 
         IA2A.bin = ifelse(`IA-2A` > 5, 1, 0), 
         mIAA.bin = ifelse(mIAA > .01, 1, 0), 
         ZnT8A.bin=ifelse(ZnT8A > .02, 1, 0),
         TGA.bin=ifelse(TGA > .05, 1, 0))

iabtable<-tn10%>%filter(Visit == 0)%>%distinct(MaskID,.keep_all = T)%>%
  select(MaskID, Treatment.Arm)%>%right_join(iabtable)

iabtable$iab.no<-rowSums(select(iabtable,GADA.bin:ZnT8A.bin),na.rm = T)
iabtable%<>%filter(!(is.na(Visit)|iab.no == 0))%>%
  mutate(iab.no = factor(iab.no))
iabtable%>%group_by(Visit)%>%
  count(iab.no)%>%
  mutate(total = sum(n), perc_iabnum = n/total)%>%
  ggplot(aes(x = Visit, y = perc_iabnum)) + geom_bar(aes(fill = iab.no),stat="identity") +
  # geom_text(aes(label = scales::percent(perc_iabnum)), position = position_stack(vjust = .5)) +
  geom_text(aes(y = 1, label = paste("n = ",total)), size = 5,vjust = -.2)+
  scale_fill_ghibli_d("LaputaMedium", direction = -1)+
  guides(fill = guide_legend(title = "No. of IABs"))+
  labs( x = "Visit", y = "% of Participants Positive for \n Given No. of IABs")
library(MASS)
mod_iab<-polr(iab.no ~ Visit+Treatment.Arm, data = iabtable, Hess=TRUE)
summary(mod_iab)
table_iab<- coef(summary(mod_iab))
p_iab <- pnorm(abs(table_iab[, "t value"]), lower.tail = FALSE) * 2
(table_iab <- cbind(table_iab, "p value" = round(p_iab,2)))
exp(coef(mod_iab))
detach("package:MASS")

acab.base<-tn10%>%filter(Visit == 0)%>%
  select(MaskID, Treatment.Arm,T1D.event,Time.to.T1D,BLbase = Bl.IgG2.RI,EFbase = Ef.IgG2.RI,DIbase = Di.IgG2.RI)

iabtable.toplot<-iabtable%>%left_join(select(acab.base, MaskID,Treatment.Arm))%>%
  mutate(Visit = factor(Visit, levels = c("Baseline","6 months","12 months","18 months","24 months","36 months"), labels = c(0,6,12,18,24,36)))%>%
  group_by(Treatment.Arm,Visit)%>%
  summarise(GADA.mean=mean(GADA,na.rm = T),
            GADA.sd=sd(GADA,na.rm = T),
            GADA.count = sum(GADA!=0, na.rm = T),
            IA2A.mean=mean(`IA-2A`,na.rm = T),
            IA2A.sd=sd(`IA-2A`,na.rm = T),
            IA2A.count = sum(`IA-2A`!=0, na.rm = T),
            mIAA.mean=mean(mIAA,na.rm = T),
            mIAA.sd=sd(mIAA,na.rm = T),
            mIAA.count = sum(mIAA!=0, na.rm = T),
            ZnT8A.mean=mean(ZnT8A,na.rm = T),
            ZnT8A.sd=sd(ZnT8A,na.rm = T),
            ZnT8A.count = sum(ZnT8A!=0, na.rm = T))%>%
  filter_at(vars(contains("count")), all_vars(.>5))%>%
  pivot_longer(-c(Treatment.Arm,Visit),names_to = c("IAB", ".value"),names_sep = "\\.")%>% # pivot_longer multiple sets of columns
  mutate(mean = log10(mean),
         sd = log10(sd),
         # moe=qlnorm(0.975,mean,sd),
         upper = mean + sd,
         lower = mean - sd)
ggplot(iabtable.toplot, aes(x =Visit,y =mean,linetype = Treatment.Arm,color = IAB,group = interaction(IAB,Treatment.Arm)))+
  geom_errorbar(aes(ymin = lower , ymax = upper),alpha = .7,position = position_dodge(width = .2))+
  geom_line(size =1)+
  ylab(expression(IAB~values~"["*scale*":"*log[10]*"]"))+xlab("On-Study (Months)")+
  # geom_text_repel(aes(label = count,group = Treatment.Arm),hjust =2,vjust=3)+
  # theme(axis.text.x = element_text(angle = 30,hjust = 1))+
  scale_colour_manual(values = ghibli_palette("SpiritedMedium")[3:7])


# ACAb overtime -----------------------------------------------------------

acab.toplot<-tn10%>%distinct(MaskID,Visit,.keep_all=T)%>%filter(!is.na(Bl.IgG2.RI))%>%select(Treatment.Arm,Visit,contains("RI"))%>%
  pivot_longer(cols = !c(Treatment.Arm,Visit))%>%split(.$name)%>%
  imap(~group_by(.,Treatment.Arm,Visit)%>%summarise(meanRI = log(mean(value, na.rm = T)),
                                                    sdRI = log(sd(value, na.rm = T)),
                                                    upper = meanRI + sdRI,
                                                    lower = meanRI - sdRI,
                                                    name = .y))%>%
  reduce(rbind)%>%mutate(Isotype = unlist(lapply((str_split(name,pattern = "\\.")), "[",2)))

acab.toplot%>%split(.$Isotype)%>%imap(~ggplot(.,aes(x =Visit,y = meanRI,linetype = Treatment.Arm,color = name,group = interaction(name,Treatment.Arm)))+
                                         geom_errorbar(aes(ymin = lower , ymax = upper),alpha = .7,position = position_dodge(width = .2))+
                                         geom_line(size =1)+
                                         labs(y = paste(.y,"responses [scale:log]"),x="On-Study (Months)",color = "Bacteria")+
                                         # geom_text_repel(aes(label = count,group = Treatment.Arm),hjust =2,vjust=3)+
                                         # scale_color_manual(labels = c(ghibli_palette("MarnieMedium1"),ghibli_palette("MarnieMedium2")[3:7],ghibli_palette("LaputaMedium")[3:7]))+
                                        scale_color_discrete(labels = c("Anaerotruncus colihominis", "Bifidobacterium animalis","Bifidobacterium bifidum",
                                                                        "Bacteroides dorei", "Bacteroides fragilis","Bifidobacterium longum",
                                                                        "Bacteroides vulgatus","Coprococcus eutactus","Dialister invisus",
                                                                        "Escherichia coli","Enterococcus faecalis","Lactobacillus acidophilus",
                                                                        "Roseburia faecis","Ruminococcus gnavus","Streptococcus gallolyticus"))+
                                        theme(legend.text = element_text(face="italic"))
)


acab.toplot2<-tn10%>%distinct(MaskID,Visit,.keep_all=T)%>%filter(!is.na(Bl.IgG2.RI))%>%select(Treatment.Arm,Visit,contains("RI"))%>%
  pivot_longer(cols = !c(Treatment.Arm,Visit))%>%separate(name,c("Bacteria","Isotype"),sep = "\\.")%>%arrange(Bacteria)
acab.toplot2$Bacteria<-rep(c("Anaerotruncus colihominis", "Bifidobacterium animalis","Bifidobacterium bifidum",
  "Bacteroides dorei", "Bacteroides fragilis","Bifidobacterium longum",
  "Bacteroides vulgatus","Coprococcus eutactus","Dialister invisus",
  "Escherichia coli","Enterococcus faecalis","Lactobacillus acidophilus",
  "Roseburia faecis","Ruminococcus gnavus","Streptococcus gallolyticus"),each = 678)
isotype_color<-c("IgA"="#00AFBB","IgG1"="#FC4E07", "IgG2"="#E7B800")
acab.toplot2%>%split(.$Isotype)%>%imap(~ggplot(.,aes(x =Visit,y =log(value),shape = Treatment.Arm, color = Isotype))+
                                        geom_boxplot()+
                                        geom_point(position = position_jitterdodge(),alpha = .5)+
                                        labs(y = paste(.y,"responses [scale:log]"),x="On-Study (Months)",color = "Isotype")+
                                        scale_color_manual(values = isotype_color[.y])+
                                        theme(legend.text = element_text(face="italic"))+facet_wrap(~Bacteria,ncol = 3)
    
)
acab.toplot2%>%
  ggplot(aes(x =Isotype,y =log(value),shape = Treatment.Arm, color = Isotype))+
    geom_boxplot()+
    geom_point(position = position_jitterdodge(),alpha = .5)+
    labs(y = "Response Index [scale:log]",x="Isotype",color = "Isotype")+
    scale_color_manual(values = isotype_color)+
    theme(legend.text = element_text(face="italic"))+facet_wrap(~Bacteria,ncol = 5)
  


tn10_ppc<-preProcess(tn10%>%select(-MaskID,-Visit,-ID,-T1D.event,-Time.to.T1D,-Time.landmark,-Sex,-matches("A$",ignore.case = F)),method = c("center", "scale","YeoJohnson"))
tn10.scaled<- predict(tn10_ppc, tn10)%>%mutate(Visit = as.numeric(levels(Visit))[Visit])

acab.scaled.long<-tn10.scaled%>%dplyr::select(MaskID,Visit,Treatment.Arm,DR3,DR4,Sex,Age,BMI,T1D.event,Time.to.T1D,Time.landmark,Bl.IgG2.RI,Ef.IgG2.RI,Di.IgG2.RI)%>%filter(!is.na(DR3))%>%pivot_longer(cols = contains("RI"))
# joint modeling with lmm -------------------------------------------------
library(nlme)
library(JM) 
ctrl <- lmeControl(opt='optim')
fit.bl <- lme(fixed = value ~ Visit, random = ~ 1 | MaskID, data = acab.scaled.long%>%filter(name == "Bl.IgG2.RI"),na.action = "na.omit")
summary(fit.bl)
# can only use cont. Visit for JM, and fit2 has lower AIC
fitSURV <- coxph(Surv(Time.landmark, T1D.event) ~ Treatment.Arm+Age+Sex+BMI+DR3+DR4, data = acab.scaled.long%>%filter(Visit == 0)%>%distinct(MaskID, .keep_all = T), x = T)
summary(fitSURV)
fitJM <- jointModel(fit.bl, fitSURV, timeVar = "Visit", method = "piecewise-PH-GH")
summary(fitJM) # HR: exp(0.6886), p = .002

fit.ef <- lme(fixed = value ~ Visit, random = ~ 1 | MaskID, data = acab.scaled.long%>%filter(name == "Ef.IgG2.RI"),na.action = "na.omit")
summary(fit.ef)
fitJM.ef <- jointModel(fit.ef, fitSURV, timeVar = "Visit", method = "piecewise-PH-GH")
summary(fitJM.ef) # HR: exp(0.5627), p = .008

fit.di <- lme(fixed = value ~ Visit, random = ~ 1| MaskID, data = acab.scaled.long%>%filter(name == "Di.IgG2.RI"),na.action = "na.omit")
summary(fit.di)
fitJM.di <- jointModel(fit.di, fitSURV, timeVar = "Visit", method = "piecewise-PH-GH",
                       interFact = list(value = ~ DR4, data = acab.scaled.long%>%filter(Visit == 0)%>%distinct(MaskID, .keep_all = T)))
summary(fitJM.di) # HR: exp(1.0269), p = 0.015

# select samples: acab.scaled.long%>%distinct(DR4,Treatment.Arm, .keep_all = TRUE)
# 211730,224910,378864,557085
set.seed(42) # it uses Monte Carlo samples
pred.bl <- survfitJM(fitJM, newdata = acab.scaled.long%>%filter(MaskID%in%c(211730,224910,378864,557085)), idVar = "MaskID")  
set.seed(42) # it uses Monte Carlo samples
pred.ef <- survfitJM(fitJM.ef, newdata = acab.scaled.long%>%filter(MaskID%in%c(211730,224910,378864,557085)), idVar = "MaskID")  
set.seed(42)
pred.di <- survfitJM(fitJM.di, newdata = acab.scaled.long%>%
                       filter(MaskID%in%c(232794,245390,378864,557085,715944,798801,873459,926463)), idVar = "MaskID")  


detach("package:JM")
detach("package:MASS")

pred.meta<-tn10%>%filter(Visit == 0)%>%distinct(MaskID, .keep_all = T)%>%dplyr::select(MaskID,Treatment.Arm,DR4, Bl.IgG2.RI,Ef.IgG2.RI,Di.IgG2.RI)%>%
  mutate(Bl.IgG2.bin = ifelse(Bl.IgG2.RI > 140,"High","Low"),
         Ef.IgG2.bin = ifelse(Ef.IgG2.RI > 159.6,"High","Low"),
         Di.IgG2.bin = ifelse(Di.IgG2.RI > 25.1,"High","Low"))
# 211730,224910,232794,378864,420795,488675,557085,672212
pred.meta%>%distinct(DR4,Treatment.Arm, Di.IgG2.bin,.keep_all = TRUE)

pred.bl$summaries%>%imap(~data.frame(.)%>%mutate(MaskID=.y))%>%bind_rows(.)%>%
  merge(pred.meta,"MaskID")%>%
  ggplot(aes(x= times, y = Mean,linetype = Bl.IgG2.bin,color = Treatment.Arm,fill=Treatment.Arm))+geom_line(size = 1.5)+
  geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=0.25)+theme_classic(base_size = 30)+
  labs(x = "On-Study (Days)", y ="Probabilities of T1D-free")+
  scale_linetype_manual(values =c("twodash","solid"),name = expression(italic("B.longum")*"-IgG2"))+
  scale_colour_manual(values = c("#2E9FDF","#E7B800"))+
  scale_fill_manual(values = c("#2E9FDF","#E7B800"))

pred.ef$summaries%>%imap(~data.frame(.)%>%mutate(MaskID=.y))%>%bind_rows(.)%>%
  merge(pred.meta,"MaskID")%>%
  ggplot(aes(x= times, y = Mean,color = Treatment.Arm,fill=Treatment.Arm,linetype = Ef.IgG2.bin))+geom_line(size = 1.5)+
  geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=0.25)+theme_classic(base_size = 30)+
  labs(x = "On-Study (Days)", y ="Probabilities of T1D-free")+
  scale_linetype_manual(values =c("twodash","solid"),name = expression(italic("E.faecalis")*"-IgG2"))+
  scale_colour_manual(values = c("#2E9FDF","#E7B800"))+
  scale_fill_manual(values = c("#2E9FDF","#E7B800"))

# pred.meta%>%distinct(DR4,Treatment.Arm, Di.IgG2.bin,.keep_all = TRUE)
p<-pred.di$summaries%>%imap(~data.frame(.)%>%mutate(MaskID=.y))%>%bind_rows(.)%>%
  merge(pred.meta)%>%split(.$DR4)%>%
  imap(~ggplot(data = .,aes(x= times, y = Mean,linetype = Di.IgG2.bin,color = Treatment.Arm,fill = Treatment.Arm))+geom_line(size = 1.5)+
  geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=0.25)+theme_classic(base_size = 30)+
  labs(x = "On-Study (Days)", y ="Probabilities of T1D-free",subtitle = bquote(italic("DR4")~.(.y)))+
  scale_linetype_manual(values =c("twodash","solid"),name = expression(italic("D.invisus")*"-IgG2"))+
  scale_colour_manual(values = c("#2E9FDF","#E7B800"))+
  scale_fill_manual(values = c("#2E9FDF","#E7B800")))

p$ABSENT+theme(plot.subtitle = element_text(size = 25))
p$PRESENT+theme(plot.subtitle = element_text(size = 25))

# total Ab over time ------------------------------------------------------
# setwd("~/Downloads/MBP/TN10/ELISA")
# change directory for elisa data
tn10elisa<-read.csv("TN10_withELISA_final.csv",as.is = T)
tnmerge<-tn10elisa%>%select("Barcode","ID","IgG1_mgml","IgG2_mgml","IgA_mgml")%>%full_join(tn10)
# tnmerge_ppc<-preProcess(tnmerge%>%select(-MaskID,-Visit,-ID,-T1D.event,-Time.to.T1D,-Sex,-matches("A$",ignore.case = F)),method = c("center", "scale","YeoJohnson"))
# tnmerge.scaled<- predict(tnmerge_ppc, tnmerge)
summary(tnmerge$IgG1_mgml)
sd(tnmerge$IgG1_mgml)
summary(tnmerge$IgG2_mgml)
sd(tnmerge$IgG2_mgml)
summary(tnmerge$IgA_mgml,na.rm = T)
sd(tnmerge$IgA_mgml,na.rm = T)

elisa.long<-tnmerge%>%select(MaskID,Visit,Treatment.Arm,T1D.event,Time.to.T1D,IgG1_mgml,IgG2_mgml,IgA_mgml)%>%
  mutate(Visit=factor(Visit))%>%pivot_longer(cols = contains("mgml"),names_to = "Isotype", values_to = "Titre")%>%mutate(Isotype = factor(Isotype, levels = c("IgG1_mgml","IgG2_mgml","IgA_mgml"), labels = c("IgG1","IgG2","IgA")))

ggplot(elisa.long,aes(x=Isotype,y = Titre, group = interaction(Isotype,Treatment.Arm),color = Isotype))+
  geom_violin(scale = "width",aes(linetype = Treatment.Arm) )+geom_point(aes(fill = Treatment.Arm),position = position_jitterdodge(1),alpha = .75,size = 2)+
  scale_color_manual(values = c("#FC4E07", "#E7B800","#00AFBB"))+labs(y = "Ab titres (mg/mL)")

elisa.toplot<-tnmerge%>%
  group_by(Treatment.Arm,Visit)%>%
  summarise(IgG1.mean=mean(IgG1_mgml,na.rm = T),
            IgG1.sd=sd(IgG1_mgml,na.rm = T),
            IgG2.mean=mean(IgG2_mgml,na.rm = T),
            IgG2.sd=sd(IgG2_mgml,na.rm = T),
            IgA.mean=mean(IgA_mgml,na.rm = T),
            IgA.sd=sd(IgA_mgml,na.rm = T))%>%
  pivot_longer(-c(Treatment.Arm,Visit),names_to = c("Isotype", ".value"),names_sep = "\\.")%>% # pivot_longer multiple sets of columns
  mutate(upper = mean + sd,
         lower = mean - sd,
         Isotype = factor(Isotype, levels = c("IgG1","IgG2","IgA")))
ggplot(elisa.toplot, aes(x =as.factor(Visit),y =mean,linetype = Treatment.Arm,color = Isotype,group = interaction(Isotype,Treatment.Arm)))+
  geom_errorbar(aes(ymin = lower , ymax = upper),alpha = .7,position = position_dodge(width = .2))+
  geom_line(size =1)+
  labs(y="Ab titres (mg/mL)",x="On-Study (Months)")+
  scale_color_manual(values = c("#FC4E07", "#E7B800","#00AFBB"))

require(nlme)
fit<-lme(IgG2_mgml ~ Visit+Treatment.Arm,random = ~ 1|MaskID,data=tnmerge)
require(lme4)
fit<-lmer(IgG2_mgml ~ Visit+Treatment.Arm+(1|MaskID),tnmerge)
summary(fit)

