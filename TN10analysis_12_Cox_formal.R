library(survival)
library(tidyverse)
library(magrittr)
library(caret)
library(survminer)
library(ghibli)
library(broom)
library(ggrepel)
library(survAUC)
tn10<-read.csv("TN10full_smallgate.csv", header = T, as.is = T)
str(tn10)
theme_set(
  theme_classic(base_size = 30)
)
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

# didn't factor Visit here because need to merge with 3 month data for time-dependent cox
tn10[is.na(tn10$Barcode),]
tn10[tn10$MaskID == "245390" & tn10$Visit == 12,"Date.IDDM"]<-tn10[tn10$MaskID == "245390" & tn10$Visit == 0,"Date.IDDM"]
# check that all individuals whose T1D.event == 0 have na Date.IDDM
tn10%>%count(Date.IDDM,T1D.event)%>%spread(T1D.event,n)

tn10_ppc<-preProcess(tn10%>%select(-MaskID,-Visit,-ID,-T1D.event,-Time.to.T1D,-Time.landmark,-Sex,-matches("A$",ignore.case = F)),method = c("center", "scale","YeoJohnson"))
tn10.scaled<- predict(tn10_ppc, tn10)


# Sanity check ------------------------------------------------------------
# survival curve analysis by group
dat1<-tn10%>%filter(Visit == 0)%>%distinct(MaskID, .keep_all = T)%>%mutate(Time.to.T1D.mon = Time.landmark/30)
dat1%>%filter(!is.na(DR3))%>%coxph(Surv(Time.to.T1D.mon,T1D.event)~Treatment.Arm+Age+BMI+Sex+DR3+DR4, data=.)

fit <- dat1%>%survfit(Surv(Time.to.T1D.mon, T1D.event) ~ Treatment.Arm, data =.)
print(fit) # for median time to event
ggsurvplot(fit,dat1,
           pval = FALSE, 
           censor.shape = NA,
           risk.table = F, # Add risk table
           # risk.table.col = "strata", # Change risk table color by groups
           # fontsize=7, # change font size of text in the table
           font.tickslab = c(20),
           break.x.by = 12,
           surv.median.line = "h", # Specify median survival
           xlab = "On-Study (Months)",
           ylab = "Proportion T1D-Free",
           ggtheme = theme_classic(base_size = 20), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"),
           tables.theme = clean_theme(),
           legend = "none"
           # legend.labs = c("Placebo", "Teplizumab"),
           # legend.title = ""
           )
tn10%>%filter(Visit == 0)%>%distinct(MaskID, .keep_all = T)%>%coxph(Surv(Time.to.T1D,T1D.event)~Treatment.Arm+Age+BMI+Sex, data = .)
exp(-0.68878+0.34733) # treatment effect only

dr3<-tn10%>%filter(Visit == 0)%>%distinct(MaskID, .keep_all = T)%>%filter(!is.na(DR3))%>%split(.$DR3)%>%
  imap(~coxph(Surv(Time.to.T1D,T1D.event)~Treatment.Arm+Age+BMI+Sex, data = .))
print(dr3$ABSENT) # check HR of treatment
print(dr3$PRESENT) # check HR of treatment
# no subsetting: DR3 p=0.0394, Treatment.Arm p=0.0189, no interaction
dr3.survfits<-tn10%>%filter(Visit == 0)%>%distinct(MaskID, .keep_all = T)%>%filter(!is.na(DR3))%>%
  split(.$DR3)%>%imap(~survfit(Surv(Time.to.T1D,T1D.event) ~ Treatment.Arm, data = .))
# survfit looks at combinatorial risk of variables, can't include secondary variables
(dr3.survfitsplots<-dr3.survfits%>%imap(~ggsurvplot(fit=., risk.table = F,pval = T,conf.int = T,
                                                    data=tn10%>%filter(Visit == 0)%>%distinct(MaskID, .keep_all = T)%>%filter(DR3 == .y),# when using imap(), .x is the value, .y is the position/index
                                                    palette = c("#E7B800", "#2E9FDF"),legend.labs = c("Placebo", "Teplizumab"),
                                                    xlab = "On-Study (days)", ylab = "Proportion T1D-Free", title = paste("DR3", .y),ggtheme = theme_classic2(base_size=30)))
  
)
dr4<-tn10%>%filter(Visit == 0)%>%distinct(MaskID, .keep_all = T)%>%filter(!is.na(DR4))%>%split(.$DR4)%>%
  imap(~coxph(Surv(Time.to.T1D,T1D.event) ~Treatment.Arm+Age+BMI+Sex, data = .))
print(dr4$ABSENT) # check HR of treatment
exp(0.77087-0.69473)
print(dr4$PRESENT) # check HR of treatment
exp(-1.62999-0.47913)
# no subsetting: DR4 p=0.0356, DR4*Treatment.Arm p=0.0029 (if no interaction DR4 p=0.9169, Treatment.Arm p=0.0423)
dr4.survfits<-tn10%>%filter(Visit == 0)%>%distinct(MaskID, .keep_all = T)%>%filter(!is.na(DR4))%>%
  split(.$DR4)%>%imap(~survfit(Surv(Time.to.T1D,T1D.event) ~ Treatment.Arm, data = .))
(dr4.survfitsplots<-dr4.survfits%>%imap(~ggsurvplot(fit=., risk.table = F,pval = T,conf.int = T,
                                                    data=tn10%>%filter(Visit == 0)%>%distinct(MaskID, .keep_all = T)%>%filter(DR4 == .y),# when using imap(), .x is the value, .y is the position/index
                                                    palette = c("#E7B800", "#2E9FDF"),legend.labs = c("Placebo", "Teplizumab"),
                                                    xlab = "On-Study (days)", ylab = "Proportion T1D-Free", title = paste("DR4", .y),ggtheme = theme_classic2(base_size=30)))
)

znt8<-tn10%>%filter(Visit == 0)%>%distinct(MaskID, .keep_all = T)%>%filter(!is.na(ZnT8A))%>%coxph(Surv(Time.to.T1D,T1D.event) ~ZnT8A*Treatment.Arm+Age+BMI+Sex, data = .)
# Treatment.Arm p=0.0103, ZnT8A*Treatment.Arm = 0.0413
(znt8.survfitsplots<-tn10%>%filter(Visit == 0)%>%distinct(MaskID, .keep_all = T)%>%
    split(.$ZnT8A)%>%imap(~survfit(Surv(Time.to.T1D,T1D.event) ~ Treatment.Arm, data = .))%>%
    imap(~ggsurvplot(fit=., risk.table = F,pval = T,conf.int = T,
                     data=tn10%>%filter(Visit == 0)%>%distinct(MaskID, .keep_all = T)%>%filter(ZnT8A == .y),
                     palette = c("#E7B800", "#2E9FDF"),legend.labs = c("Placebo", "Teplizumab"),
                     xlab = "On-Study (days)", ylab = "Proportion T1D-Free", title = paste("ZnT8A", .y),,ggtheme = theme_classic2(base_size=30)))
)
cpep.base<-tn10%>%filter(Visit == 0)%>%distinct(MaskID, .keep_all = T)%>%filter(!is.na(mu_cpep)&!is.na(DR3))%>%
  coxph(Surv(Time.to.T1D,T1D.event) ~mu_cpep+Treatment.Arm+Age+BMI+Sex+DR3+DR4, data = .)
# mu_cpep p=0.0118, Treatment.Arm p=0.0177 (only Age and BMI sign. if interaction)
tn10.bincpep<-tn10%>%filter(Visit == 0)%>%distinct(MaskID, .keep_all = T)%>%mutate(cpep=ifelse(mu_cpep > median(mu_cpep), "> median","< median"))
(cpep.survfitsplots<-tn10.bincpep%>%split(.$cpep)%>%imap(~survfit(Surv(Time.to.T1D,T1D.event) ~ Treatment.Arm, data = .))%>%
    imap(~ggsurvplot(fit=., risk.table = F,pval = T,conf.int = T,
                     data=tn10.bincpep%>%filter(cpep == .y),
                     palette = c("#E7B800", "#2E9FDF"),legend.labs = c("Placebo", "Teplizumab"),
                     xlab = "On-Study (days)", ylab = "Proportion T1D-Free", title = paste("C peptide", .y),ggtheme = theme_classic2(base_size=30)))
)


# Baseline Cox ------------------------------------------------------------
# Bl.IgG2 response
bl.igg2.bu<-tn10%>%filter(Visit == 0)%>%distinct(MaskID, .keep_all = T)%>%
  filter(!is.na(DR3))%>%
  coxph(Surv(Time.landmark,T1D.event) ~Bl.IgG2.RI+Treatment.Arm+Age+BMI+Sex+DR3+DR4, data = .)
# Treatment.Arm p=0.03, Bl.IgG2 p=0.01
cox.zph(bl.igg2.bu) 
plot(cox.zph(bl.igg2.bu)[1])
summary(bl.igg2.bu)

bl.igg2.bu.add<- mgcv::gam(Time.to.T1D~ s(Bl.IgG2.RI) + DR3+DR4+Treatment.Arm+Age+BMI+Sex, 
                              data = tn10%>%filter(Visit == 0)%>%distinct(MaskID, .keep_all = T)%>%
                             filter(!is.na(DR3)), family = "cox.ph", weights = T1D.event)
summary(bl.igg2.bu.add)

ggplot(tn10, aes(x=Bl.IgG2.RI))+geom_histogram(bins = 25, fill ="#92BBD9FF",color="#4D6D93FF")+geom_vline(xintercept  = 140,linetype = "dashed")+labs(x=expression("All-time"~italic("B.longum")*"-IgG2"), y = "Sample count")
tn10%>%filter(Visit == 0)%>%distinct(MaskID, .keep_all = T)%>%
  ggplot(., aes(x=Bl.IgG2.RI))+geom_histogram(bins = 18, fill ="#92BBD9FF",color="#4D6D93FF")+geom_vline(xintercept  = 140,linetype = "dashed")+labs(x= expression(Baseline~italic("B.longum")*"-IgG2"), y = "Sample count")
tn10%>%filter(Visit == 6)%>%distinct(MaskID, .keep_all = T)%>%
  ggplot(., aes(x=Bl.IgG2.RI))+geom_histogram(bins = 13, fill ="#92BBD9FF",color="#4D6D93FF")+geom_vline(xintercept  = 140,linetype = "dashed")+labs(x= expression("6-month"~italic("B.longum")*"-IgG2"), y = "Sample count")

tn10.bligg2.b<-tn10%>%mutate(bligg2=ifelse(Bl.IgG2.RI > 140, "high","low"))%>%filter(Visit == 0)%>%distinct(MaskID, .keep_all = T)
table(tn10.bligg2.b$bligg2,tn10.bligg2.b$Treatment.Arm)
# for forest plot
tn10.bligg2.b%>%split(.$bligg2)%>%imap(~coxph(Surv(Time.landmark,T1D.event) ~ Treatment.Arm+Age+BMI+Sex+DR3+DR4, data = .))%>%walk(print)
exp(-1.71636-0.60443) # high Bl.IgG2
exp(-0.19115-0.55873) # low Bl.IgG2
# KM
(bligg2.bu.survfitsplots<-tn10.bligg2.b%>%split(.$bligg2)%>%imap(~survfit(Surv(Time.landmark,T1D.event) ~ Treatment.Arm, data = .))%>%
    imap(~ggsurvplot(fit=., risk.table = F,conf.int = T,censor.shape = NA,break.time.by = 500,
                     data=tn10.bligg2.b%>%filter(bligg2 == .y),ggtheme = theme_classic(base_size = 30),
                     palette = c("#E7B800", "#2E9FDF"),legend = "right",legend.labs = c("Placebo", "Teplizumab"),legend.title = bquote("Baseline "*italic("B.longum")*"-IgG2"~ .(.y)),
                     xlab = "On-Study (Days)", ylab = "Proportion T1D-Free"))
)
# all four strata on same plot
p<-tn10.bligg2.b%>%survfit(Surv(Time.landmark,T1D.event) ~ Treatment.Arm+bligg2, data = .)%>%
  ggsurvplot(fit=., risk.table = F,pval = T,conf.int = F,censor.shape = NA,
             data=tn10.bligg2.b,break.time.by = 500,# size = 2,
             linetype = c("bligg2"),
             ggtheme = theme_classic(base_size = 30),
             legend = "right",legend.labs = rep(c("Placebo","Teplizumab"),each = 2),
             xlab = "On-Study (days)", ylab = "Proportion T1D-Free")+guides(colour = guide_legend(nrow = 4))
p$plot+scale_linetype_manual(name= expression(paste("Baseline ",italic("B.longum"),"-IgG2")),values = c(6,1))+
  scale_color_manual(values = rep(c("#E7B800", "#2E9FDF"),2),name = "Treatment Arm")# +theme(plot.margin = unit(x = c(0, 0, 0, 5), units = "mm"))
pairwise_survdiff(Surv(Time.landmark,T1D.event) ~ Treatment.Arm+bligg2,p.adjust.method = "BH",
                         data = tn10.bligg2.b)

# Ef.IgG2 response
ef.igg2.bu<-tn10%>%filter(Visit == 0)%>%distinct(MaskID, .keep_all = T)%>%filter(!is.na(DR3))%>%
  coxph(Surv(Time.landmark,T1D.event) ~Ef.IgG2.RI+Treatment.Arm+Age+BMI+Sex+DR3+DR4, data = .)
cox.zph(ef.igg2.bu) 
# Treatment.Arm p=0.023, Ef.IgG2 p=0.003, DR3 p = 0.015

ggplot(tn10, aes(x=Ef.IgG2.RI))+geom_histogram(bins = 25, fill ="#92BBD9FF",color="#4D6D93FF")+geom_vline(aes(xintercept  = median(Ef.IgG2.RI, na.rm = T)),linetype = "dashed")+labs(x=expression("All-time"~italic("E.faecalis")*"-IgG2"), y = "Sample count")
tn10%>%filter(Visit == 0)%>%distinct(MaskID, .keep_all = T)%>%
  ggplot(., aes(x=Ef.IgG2.RI))+geom_histogram(bins = 18, fill ="#92BBD9FF",color="#4D6D93FF")+geom_vline(xintercept = 159.6,linetype = "dashed")+labs(x= expression(Baseline~italic("E.faecalis")*"-IgG2"), y = "Sample count")
tn10%>%filter(Visit == 6)%>%distinct(MaskID, .keep_all = T)%>%
  ggplot(., aes(x=Ef.IgG2.RI))+geom_histogram(bins = 13, fill ="#92BBD9FF",color="#4D6D93FF")+geom_vline(xintercept  = 140,linetype = "dashed")+labs(x= expression("6-month"~italic("E.faecalis")*"-IgG2"), y = "Sample count")

tn10.efigg2.b<-tn10%>%mutate(efigg2=ifelse(Ef.IgG2.RI > median(Ef.IgG2.RI,na.rm = T), "high","low"))%>%filter(Visit == 0)%>%distinct(MaskID, .keep_all = T)
table(tn10.efigg2.b$efigg2)
# for forest plot
tn10.efigg2.b%>%split(.$efigg2)%>%imap(~coxph(Surv(Time.landmark,T1D.event) ~ Treatment.Arm+Age+BMI+Sex+DR3+DR4, data = .))%>%walk(print)
exp(-0.85320-0.53339) # high Ef.IgG2
exp(-0.24962-0.66601) # low Ef.IgG2
# KM
(efigg2.bu.survfitsplots<-tn10.efigg2.b%>%split(.$efigg2)%>%imap(~survfit(Surv(Time.landmark,T1D.event) ~ Treatment.Arm, data = .))%>%
    imap(~ggsurvplot(fit=., risk.table = F,pval = F,conf.int = T,censor.shape = NA,break.time.by = 500,
                     data=tn10.efigg2.b%>%filter(efigg2 == .y),ggtheme = theme_classic(base_size = 30),
                     palette = c("#E7B800", "#2E9FDF"),legend = "right",legend.labs = c("Placebo", "Teplizumab"),legend.title = bquote("Baseline "*italic("E.faecalis")*"-IgG2" ~ .(.y)),
                     xlab = "On-Study (Days)", ylab = "Proportion T1D-Free"))
)


# 6 month + baseline Cox  -------------------------------------------------
### unscaled ###
base.unscaled<-tn10%>%filter(Visit == 0)%>%distinct(MaskID, .keep_all = T)%>%filter(!is.na(DR3))%>%select(-ID,-Barcode,-Date.IDDM,-Visit)
time.indep.varname<-c("MaskID", "Time.landmark","Time.to.T1D", "T1D.event", "Age", "Sex","Treatment.Arm", "DR3", "DR4")
oldnames<-colnames(base.unscaled)[!(colnames(base.unscaled)%in%time.indep.varname)]
newnames<-paste(oldnames,"T0",sep=".")
base.unscaled%<>%rename_at(all_of(oldnames), ~ newnames)


month6.unscaled<-tn10%>%filter(Visit == 6)%>%filter(!is.na(DR3))%>%select(-ID,-Barcode,-Date.IDDM,-Visit)
newnames.6<-paste(oldnames,"T6",sep=".")
month6.unscaled%<>%rename_at(vars(oldnames), ~ newnames.6)

base.month6.unscaled<-left_join(month6.unscaled,base.unscaled)

# Bl.IgG2 
bl.igg2.6u<-coxph(Surv(Time.landmark,T1D.event) ~Bl.IgG2.RI.T6+Treatment.Arm+DR3+DR4+Age+BMI.T6+Sex, data = base.month6.unscaled)
# Treatment.Arm p=0.03, Bl.IgG2.RI6=0.0069
cox.zph(bl.igg2.6u)
bl.igg2.6u.add<- mgcv::gam(Time.landmark~s(Bl.IgG2.RI.T6)+DR3+DR4+Treatment.Arm+Age+BMI.T6+Sex, 
                              data = base.month6.unscaled, family = "cox.ph", weights = T1D.event)
summary(bl.igg2.6u.add)

# Ef.IgG2 
ef.igg2.6u<-coxph(Surv(Time.landmark,T1D.event) ~Ef.IgG2.RI.T6+Treatment.Arm+DR3+DR4+Age+BMI.T0+Sex, data = base.month6.unscaled)
# Treatment.Arm p=0.03, Ef.IgG2.RI6=0.0009)
cox.zph(ef.igg2.6u)



### scaled ###
base.scaled<-tn10.scaled%>%filter(Visit == 0)%>%distinct(MaskID, .keep_all = T)%>%filter(!is.na(DR3))%>%select(-ID,-Barcode,-Date.IDDM,-Visit)
# newnames<-paste(oldnames,"T0",sep=".")
base.scaled%<>%rename_at(vars(oldnames), ~ newnames)

month6.scaled<-tn10.scaled%>%filter(Visit == 6)%>%filter(!is.na(DR3))%>%select(-ID,-Barcode,-Date.IDDM,-Visit)
# newnames.6<-paste(oldnames,"T6",sep=".")
month6.scaled%<>%rename_at(vars(oldnames), ~ newnames.6)

base.month6.scaled<-left_join(month6.scaled,base.scaled)


# plotting 6 month ACAbs survival analysis
tn10.bligg2<-full_join(month6.unscaled,base.unscaled)%>%mutate(bligg2.0=ifelse(Bl.IgG2.RI.T0 > 140, "High","Low"),
                                                               bligg2.6=ifelse(Bl.IgG2.RI.T6 > 140, "High","Low"))
table(tn10.bligg2$bligg2.0,tn10.bligg2$bligg2.6)
(bligg2.6u.survfitsplots<-tn10.bligg2%>%split(.$bligg2.6)%>%imap(~survfit(Surv(Time.landmark,T1D.event) ~ Treatment.Arm, data = .))%>%
    imap(~ggsurvplot(fit=., risk.table = F,conf.int = T,censor.shape = NA,break.time.by = 500,
                     data=tn10.bligg2%>%filter(bligg2.6 == .y),ggtheme = theme_classic(base_size = 30),
                     palette = c("#E7B800","#2E9FDF"),legend = "right",legend.labs = c("Placebo", "Teplizumab"),legend.title = bquote("6-month "*italic("B.longum")*"-IgG2"~ .(.y)),
                     xlab = "On-Study (Days)", ylab = "Proportion T1D-Free"))
)
# KM analysis of surv prob in high vs low group (not adjusted for other var)
# ggsurvplot(survfit(Surv(Time.to.T1D,T1D.event) ~ bligg2.0, data = tn10.bligg2),pval=T,conf.int = T,palette = c("#E7B800","#2E9FDF"),
#            title = "Baseline Bl-IgG2 association with T1D risk", xlab = "On-Study (Days)", ylab = "Proportion T1D-Free",
#            legend.labs = c("Bl-IgG2 > 140", "Bl-IgG2 < 140"),ggtheme = theme_classic2(base_size=20))
# ggsurvplot(survfit(Surv(Time.to.T1D,T1D.event) ~ bligg2.6, data = tn10.bligg2),pval=T,conf.int = T,palette = c("#E7B800","#2E9FDF"),
#            title = "6 month Bl-IgG2 association with T1D risk", xlab = "On-Study (Days)", ylab = "Proportion T1D-Free",
#            legend.labs = c("Bl-IgG2 > 140", "Bl-IgG2 < 140"),ggtheme = theme_classic2(base_size=20))


ef2.median<-median(tn10$Ef.IgG2.RI,na.rm = T)
tn10.efigg2<-full_join(month6.unscaled,base.unscaled)%>%mutate(efigg2.0=ifelse(Ef.IgG2.RI.T0 > ef2.median, "High","Low"),
                                                               efigg2.6=ifelse(Bl.IgG2.RI.T6 > ef2.median, "High","Low"))

(efigg2.6u.survfitsplots<-tn10.efigg2%>%split(.$efigg2.6)%>%imap(~survfit(Surv(Time.landmark,T1D.event) ~ Treatment.Arm, data = .))%>%
    imap(~ggsurvplot(fit=., risk.table = F,conf.int = T,censor.shape = NA,break.time.by = 500,
                     data=tn10.efigg2%>%filter(efigg2.6 == .y),ggtheme = theme_classic(base_size = 30),
                     palette = c("#E7B800", "#2E9FDF"),legend = "right",legend.labs = c("Placebo", "Teplizumab"),
                     legend.title = bquote("6-month "*italic("E.faecalis")*"-IgG2"~ .(.y)),
                     xlab = "On-Study (Days)", ylab = "Proportion T1D-Free"))
)

# ggsurvplot(survfit(Surv(Time.to.T1D,T1D.event) ~ efigg2.0, data = tn10.efigg2),pval=T, conf.int = T, 
#            palette = c("#2E9FDF","#E7B800"),title = "Baseline Ef IgG2 association with T1D risk",xlab = "On-Study (days)", ylab = "Proportion T1D-Free",
#            legend.labs = c("Ef IgG2 < median", "Ef IgG2 > median"),ggtheme = theme_classic2(base_size=20))
# ggsurvplot(survfit(Surv(Time.to.T1D,T1D.event) ~ efigg2.6, data = tn10.efigg2),pval=T, conf.int = T, 
#            palette = c("#2E9FDF","#E7B800"),title = "6 month Ef IgG2 association with T1D risk",xlab = "On-Study (days)", ylab = "Proportion T1D-Free",
#            legend.labs = c("Ef IgG2 < median", "Ef IgG2 > median"),ggtheme = theme_classic2(base_size=20))


models.6s<-full_join(month6.scaled,base.scaled)%>%select(-Date_of_Draw.T6,-Date_of_Draw.T0,-glycemia.T0,-glycemia.T6,-contains("A.T"))%>%
  pivot_longer(cols = !c(time.indep.varname,"BMI.T0","BMI.T6"), names_to = c("feature","time"), names_pattern = "(.+)\\.(T\\d)")%>%
  pivot_wider(names_from = time, values_from = value )%>%
  nest_by(feature)%>%
  mutate(cox6s=list(coxph(Surv(Time.landmark,T1D.event) ~T0+T6*DR4+Treatment.Arm+DR3+Age+BMI.T6+Sex,data=data)))%>%
  summarise(tidy(cox6s))
models.6s%>%group_by(term)%>%filter(p.value < .05)%>%count()
models.6s%>%filter(term == "T6" & p.value < .05/45)
models.6s%>%filter(term == "T6:DR4PRESENT" & p.value < .05/41.6573)
models.6s%>%filter(feature == "Di.IgG2.RI")

di2.median <-median(tn10.scaled$Di.IgG2.RI,na.rm = T)
ggplot(tn10.scaled, aes(x=Di.IgG2.RI))+geom_histogram(bins = 25, fill ="#92BBD9FF",color="#4D6D93FF")+geom_vline(xintercept  = di2.median,linetype = "dashed")+labs(x=expression("All-time"~italic("D.invisus")*"-IgG2"), y = "Sample count")
tn10.scaled%>%filter(Visit == 6)%>%
  ggplot(., aes(x=Di.IgG2.RI))+geom_histogram(bins = 18, fill ="#92BBD9FF",color="#4D6D93FF")+geom_vline(xintercept  = di2.median,linetype = "dashed")+labs(x= expression("6-month"~italic("D.invisus")*"-IgG2"), y = "Sample count")

tn10.diigg2<-full_join(month6.scaled,base.scaled)%>%mutate(diigg2.6=ifelse(Di.IgG2.RI.T6 > di2.median, "High","Low"),
                                                           diigg2.0=ifelse(Di.IgG2.RI.T0 > di2.median, "High","Low"))

table(tn10.diigg2$diigg2.0,tn10.diigg2$DR4)

# forest plot
# tn10.diigg2%>%split(.$diigg2.0)%>%map(~split(.,.$DR4))%>%imap(~print(.y))
tn10.diigg2%>%split(.$DR4)%>%imap(~coxph(Surv(Time.landmark,T1D.event) ~ Treatment.Arm+Age+BMI.T0+Sex+DR3, data = filter(.,diigg2.0=="Low")))%>%walk(print)
exp(-0.84001-1.27708) # DR4-absent, low Di.IgG2
exp(-1.5832+1.2281) # DR4-present, low Di.IgG2
tn10.diigg2%>%split(.$DR4)%>%imap(~coxph(Surv(Time.landmark,T1D.event) ~ Treatment.Arm+Age+BMI.T0+Sex+DR3, data = filter(.,diigg2.0=="High")))%>%walk(print)
exp(1.2564-1.8197) # DR4-absent, high Di.IgG2
exp(-2.90287+1.07136) # DR4-present, high Di.IgG2

cols <-  rep(c("#E7B800", "#2E9FDF"), 2)
lines <-  rep(c("twodash","solid"), each = 2)
(diigg2.6s.survfitsplots<-tn10.diigg2%>%split(.$DR4)%>%imap(~survfit(Surv(Time.landmark,T1D.event) ~ diigg2.6+Treatment.Arm, data = .))%>%
    imap(~ggsurvplot(fit=., risk.table = F,pval = F,conf.int = F,censor.shape = NA,break.time.by = 500,
                     data=tn10.diigg2%>%filter(DR4 == .y),ggtheme = theme_classic(base_size = 30),
                     linetype = c("strata"),legend = "right",
                     legend.labs = c("High,Placebo", "High,Teplizumab","Low,Placebo","Low,Teplizumab"),
                     legend.title = expression("6-month "*italic("D.invisus")*"-IgG2"), 
                     xlab = "On-Study (Days)", ylab ="Proportion T1D-Free")+guides(color = guide_legend(nrow = 4)))
)
diigg2.6s.survfitsplots$ABSENT$plot+
  scale_linetype_manual(values = lines)+
  scale_colour_manual(values = cols)+
  labs(subtitle = substitute(italic("HLA-DR4")~x, list(x="Absent")))
  # guides(color=guide_legend(title=bquote(italic("HLA-DR")~"Absent, 6-month "*italic("D.invisus")*"-IgG2")),
  #        linetype=guide_legend(title=bquote(italic("HLA-DR")~"Absent, 6-month "*italic("D.invisus")*"-IgG2")))

diigg2.6s.survfitsplots$PRESENT$plot+
  scale_linetype_manual(values = lines)+
  scale_colour_manual(values = cols)+
  labs(subtitle = substitute(italic("HLA-DR4")~x, list(x="Present")))
models.6s%>%filter(feature == "Di.IgG2.RI")

# coxph(Surv(Time.to.T1D,T1D.event) ~ Di.IgG2.RI.T0*DR4, data = tn10.diigg2)%>%broom::tidy()
# (diigg2.0s.survfitsplots<-tn10.diigg2%>%survfit(Surv(Time.to.T1D,T1D.event) ~ diigg2.0+DR4+Treatment.Arm, data = .)%>%
#     ggsurvplot(fit=., risk.table = F,pval = T,conf.int = F,
#                data=tn10.diigg2,facet.by = "DR4",
#                palette = ghibli_palette("YesterdayMedium")[c(5,6,8,7)],legend.labs = c("< median,Placebo","< median,Teplizumab","> median,Placebo", "> median,Teplizumab"),legend.title = "Di IgG2",
#                xlab = "On-Study (days)", ylab = "Proportion T1D-Free", ggtheme = theme_classic2(base_size=30))+guides(colour = guide_legend(nrow = 2)))
# 
# coxph(Surv(Time.to.T1D,T1D.event) ~ Di.IgG2.RI.T6*DR4+Treatment.Arm, data = tn10.diigg2)%>%summary()
# (diigg2.6s.survfitsplots2<-tn10.diigg2%>%survfit(Surv(Time.to.T1D,T1D.event) ~ diigg2.6+DR4, data = .)%>%
#     ggsurvplot(fit=., risk.table = F,pval = T,conf.int = T,
#                data=tn10.diigg2,facet.by = "DR4",
#                palette = c("#E7B800", "#2E9FDF"),legend.title = "Di IgG2",
#                xlab = "On-Study (days)", ylab = "Proportion T1D-Free", ggtheme = theme_classic2(base_size=30))+guides(colour = guide_legend(nrow = 2)))


# Sanity check, mu_cpep
cpep.median<-median(tn10$mu_cpep,na.rm = T)
tn10.cpep<-full_join(month6.unscaled,base.unscaled)%>%mutate(cpep.0=ifelse(mu_cpep.T0 > cpep.median, "High","Low"),
                                                             cpep.6=ifelse(mu_cpep.T6 > cpep.median, "High","Low"))

ggplot(tn10, aes(x=mu_cpep))+geom_histogram(bins = 25, fill ="#92BBD9FF",color="#4D6D93FF")+geom_vline(xintercept  = cpep.median,linetype = "dashed")+labs(x= "All-time C-peptide AUC", y = "Sample count")
tn10%>%filter(Visit == 0)%>%distinct(MaskID, Visit, .keep_all =T)%>%
  ggplot(., aes(x=mu_cpep))+geom_histogram(bins = 18, fill ="#92BBD9FF",color="#4D6D93FF")+geom_vline(xintercept  = cpep.median,linetype = "dashed")+labs(x= "Baseline C-peptide AUC", y = "Sample count")
tn10%>%filter(Visit == 6)%>%distinct(MaskID, Visit, .keep_all =T)%>%
  ggplot(., aes(x=mu_cpep))+geom_histogram(bins = 18, fill ="#92BBD9FF",color="#4D6D93FF")+geom_vline(xintercept  = cpep.median,linetype = "dashed")+labs(x= "6-month C-peptide AUC", y = "Sample count")

cpep.6u<-coxph(Surv(Time.landmark,T1D.event) ~mu_cpep.T6+DR3+DR4+Treatment.Arm+Age+BMI.T6+Sex, data = base.month6.unscaled)
(cpep.6u.survfitsplots<-tn10.cpep%>%split(.$cpep.6)%>%imap(~survfit(Surv(Time.landmark,T1D.event) ~ Treatment.Arm, data = .))%>%
    imap(~ggsurvplot(fit=., risk.table = F,conf.int = T,censor.shape = NA,break.time.by = 500,
                     data=tn10.cpep%>%filter(cpep.6 == .y),ggtheme = theme_classic(base_size = 30),
                     palette = c("#E7B800", "#2E9FDF"),legend = "right",legend.labs = c("Placebo", "Teplizumab"),
                     xlab = "On-Study (Days)", ylab = "Proportion T1D-Free", legend.title = paste("6-month C-peptide AUC", .y)))
)
(cpep.0u.survfitsplots<-tn10.cpep%>%split(.$cpep.0)%>%imap(~survfit(Surv(Time.landmark,T1D.event) ~ Treatment.Arm, data = .))%>%
    imap(~ggsurvplot(fit=., risk.table = F,conf.int = T,censor.shape = NA,break.time.by = 500,ggtheme = theme_classic(base_size = 30),
                     data=tn10.cpep%>%filter(cpep.0 == .y),
                     palette = c("#E7B800", "#2E9FDF"),legend = "right",legend.labs = c("Placebo", "Teplizumab"),
                     xlab = "On-Study (Days)", ylab = "Proportion T1D-Free", legend.title = paste("Baseline C-peptide AUC", .y)))
)



# Time dependent Cox (tdc) --------------------------------------------
# This analysis is for time dependent covariates, not time dependent coeffcients -- The proportional hazards model estimates an average hazard over time
# an alternative to landmark analysis, using a covariate that is measured after follow-up time begins: https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html#Landmark_approach
# manipulate original table to calculate days from dates
tn10fortdc.unscaled<-tn10%>%filter(!is.na(DR3))%>%
  distinct(MaskID,Visit,.keep_all = T)%>%
  mutate(Date.start = as.Date(ifelse(Visit == 0, Date_of_Draw, NA),origin="1970-01-01"))%>%
  arrange(MaskID,Visit)%>%
  fill(Date.start,ID,Age,Sex,DR3,DR4)%>%
  mutate(futime = Date.IDDM - Date.start,
         futime = ifelse(is.na(futime),Time.to.T1D+2, futime),
         Date_of_Draw = as.Date(ifelse(is.na(Date_of_Draw), Date.start+30*Visit,Date_of_Draw),origin="1970-01-01"),
         days = Date_of_Draw - Date.start)%>% # 3 individuals don't have ACAb but have cpep data so no date of draw
  select(MaskID, Date_of_Draw,Date.start,futime,days,Date.IDDM,Time.to.T1D,T1D.event,everything())

tn10.tdcu<-subset(tn10fortdc.unscaled, select=c(time.indep.varname,"futime"))%>%
  distinct(MaskID, futime,.keep_all =T)%>% # only keep one obervation per subject
  tmerge(.,.,id=MaskID,death=event(futime,T1D.event)) # set range
tn10.tdcu<-tmerge(tn10.tdcu,tn10fortdc.unscaled,id=MaskID,BMI = tdc(days,BMI))%>%
  tmerge(tn10fortdc.unscaled, id=MaskID, Bl.IgG2.RI = tdc(days, Bl.IgG2.RI))%>%
  tmerge(tn10fortdc.unscaled, id=MaskID, Ef.IgG2.RI = tdc(days, Ef.IgG2.RI))%>%
  tmerge(tn10fortdc.unscaled, id=MaskID, mu_cpep = tdc(days, mu_cpep))

plot_res_and_HR<-function(df,var_name,other_var=c("DR3","DR4","Treatment.Arm","Age","BMI","Sex")){
  df%<>%select(var_name,other_var,tstart,tstop,death)
  tdc<-coxph(Surv(tstart, tstop, death) ~ .,data=df)
  var_mean<-tdc$means[var_name]
  residual<-cox.zph(tdc)
  print(residual)
  output<-vector(mode = "list", length = 3)
  output[[1]]<-tdc
  output[[2]]<-ggplot(data.frame(cbind(residual$y[,var_name],time = residual$time)), aes(x= time, y = V1))+
    geom_point(size = 3, col = "#92BBD9FF")+geom_smooth(method="lm",formula =y ~ x,col = "#4D6D93FF",fill = "#4D6D93FF")+ 
    labs(x="On-Study(days)",y=paste("Coefficient for", var_name,"Overtime"))+
    geom_hline(aes(yintercept = 0), linetype = 3,size = 1,color = "#505050") +
    scale_x_continuous(trans=scales::modulus_trans(.5),breaks = seq(0,2500,500))# +theme_linedraw(base_size =20) 
  tdc_termplot<-termplot(tdc,se = T, plot = F)
  var_ref<-tdc_termplot[[var_name]]$y[which.min(abs(tdc_termplot[[var_name]]$x - var_mean))]
 output[[3]]<-ggplot(tdc_termplot[[var_name]],aes(x=x)) +
    geom_line(aes(y = exp(y - var_ref)),size = 1.5,color = "#4D6D93FF") +
    geom_line(aes(y = exp(y-1.96*se - var_ref)), linetype = 2,color = "#4D6D93FF",size = .1) +
    geom_line(aes(y = exp(y+1.96*se - var_ref)), linetype = 2,color = "#4D6D93FF",size = .1) +
    geom_ribbon(aes(ymin=exp(y-1.96*se - var_ref),ymax=exp(y+1.96*se - var_ref)), fill = "#4D6D93FF",alpha=0.5)+
    geom_hline(aes(yintercept = 1), linetype = 3,size = 1,color = "#505050") +
    geom_vline(aes(xintercept = var_mean), linetype = 2,size = .5,color = "#505050")+ 
    geom_rug(data = tn10fortdc.unscaled, aes(x = tn10fortdc.unscaled[,var_name]), sides = "b",color="#4D6D93FF") +
    labs(x = var_name,
         y = paste("Hazard Ratio compared to mean",var_name)# ,title = paste("T1D Hazard Ratio as a function of",var_name)
         )# + theme_linedraw(base_size =20)
  return(output)
}
bl.tdcu.model<-plot_res_and_HR(tn10.tdcu, "Bl.IgG2.RI")
bl.tdcu.model[[3]]+labs(y = "Hazard Ratio", x = expression("All-time"~italic("B.longum")*"-IgG2"))+theme(plot.margin = unit(c(.5,.5,.5,.5), "cm"))
ef.tdcu.model<-plot_res_and_HR(tn10.tdcu, "Ef.IgG2.RI")
ef.tdcu.model[[3]]+labs(y = "Hazard Ratio", x = expression("All-time"~italic("E.faecalis")*"-IgG2"))
cpep.tdcu.model<-plot_res_and_HR(tn10.tdcu, "mu_cpep")
cpep.tdcu.model[[3]]+labs(y = "Hazard Ratio", x = "All-time C-peptide AUC")

# residual.cpep<-cox.zph(cpep.tdcu.model[[1]])
# cox.zph(cpep.tdcu.model[[1]], transform= function(time) log(time/365.25))
# 
# cpep.tdcu.model2<-coxph(Surv(tstart, tstop, death) ~ mu_cpep+tt(mu_cpep)+Treatment.Arm+Age+BMI+Sex+DR3+DR4,
#                         data=tn10.tdcu,
#                         tt =function(x,t,...) x*log(t/365.25))
# plot(residual.cpep[1]);
# abline(0,0,col=2)
# abline(coef(cpep.tdcu.model2)[1:2],col=3)
# 
# ggplot(data.frame(cbind(residual.cpep$y[, "mu_cpep"],time = residual.cpep$time)), aes(x= time, y = V1))+
#   geom_point(size = 3, col = "#92BBD9FF")+geom_smooth(method="lm",formula =y ~ x,col = "#4D6D93FF",fill = "#4D6D93FF")+ 
#   labs(x="On-Study(days)",y=paste("Coefficient for mu_cpep Overtime"))+
#   geom_hline(aes(yintercept = 0), linetype = 3,size = 1,color = "#505050") +
#   scale_x_continuous(trans=scales::trans_new("by_year",function(x) log(x/365.25),function(x) exp(x)*365.25),breaks = seq(0,2500,500))+
#   theme_linedraw(base_size =20) 

# cpep.pterms<-predict(cpep.tdcu.model2,se.fit = T, type="terms") 
# cpep.mean<-cpep.tdcu.model2$means["mu_cpep"]
# cpep.termplot <- data.frame(cpep = cpep.pterms[["fit"]][,1],cpep.se = cpep.pterms[["se.fit"]][,1])%>%rownames_to_column()%>%
#   mutate(rowname = ifelse(grepl("\\.",rowname),rowname,paste0(rowname,".0")))%>%
#   separate(rowname,c("ID","time"))%>%mutate(ID = as.integer(ID), time = as.integer(time))
# cpep.termplot <-cbind(cpep.termplot,x=tn10.tdcu$mu_cpep[cpep.termplot$ID])
# cpep.ref <- cpep.termplot$cpep[which.min(abs(cpep.termplot$x - cpep.mean))]
# 
# cpep.termplot%>%filter(time <= 4)%>% mutate(time = paste("timepoint",time))%>%
# ggplot(aes(x=x)) +
#   geom_line(aes(y = exp(cpep - cpep.ref)),size = 1.5,color="#4D6D93FF") +
#   geom_line(aes(y = exp(cpep-1.96*cpep.se -cpep.ref)), linetype = 2,size = .1) +
#   geom_line(aes(y = exp(cpep+1.96*cpep.se -cpep.ref)), linetype = 2,size = .1) +
#  geom_ribbon(aes(ymin=exp(cpep-1.96*cpep.se - cpep.ref),ymax=exp(cpep+1.96*cpep.se - cpep.ref)),alpha=0.25,fill = "#4D6D93FF")+
#   geom_hline(aes(yintercept = 1), linetype = 3,size = 1,color = "#505050") +
#   geom_vline(aes(xintercept = cpep.mean), linetype = 2,size = .5,color = "#505050")+ 
#   geom_rug(data = tn10.tdcu, aes(x = mu_cpep), sides = "b",color="#4D6D93FF") +
#   labs(x = "C-peptide AUC",
#        y = "Hazard Ratio compared to mean C-peptide AUC")

tn10fortdc.scaled<-tn10.scaled%>%filter(!is.na(DR3))%>%
#  left_join(month3.scaled)%>%
  distinct(MaskID,Visit,.keep_all = T)%>%
  mutate(Date.start = as.Date(ifelse(Visit == 0, Date_of_Draw, NA),origin="1970-01-01"))%>%
  arrange(MaskID,Visit)%>%
  fill(Date.start,ID,Age,Sex,DR3,DR4)%>%
  mutate(futime = Date.IDDM - Date.start,
         futime = ifelse(is.na(futime),Time.to.T1D+2, futime),
         Date_of_Draw = as.Date(ifelse(is.na(Date_of_Draw), Date.start+30*Visit,Date_of_Draw),origin="1970-01-01"),
         days = Date_of_Draw - Date.start)%>% # 3 individuals don't have ACAb but have cpep data so no date of draw
  select(MaskID, Date_of_Draw,Date.start,futime,days,Date.IDDM,Time.to.T1D,T1D.event,everything())
# check that things are filled in properly, one individual should have the same Treatment.Arm, Date.IDDM, Age, Sex ...
nSamplePerSub<-tn10fortdc.scaled%>%group_by(MaskID,Treatment.Arm,Date.IDDM,Time.to.T1D,Date.start,ID,Age,Sex,DR3,DR4)%>%count()

tn10.tdcs<-subset(tn10fortdc.scaled, select=c(time.indep.varname,"futime"))%>%
  distinct(MaskID, futime,.keep_all =T)%>% # only keep one obervation per subject
  tmerge(.,.,id=MaskID,death=event(futime,T1D.event)) # set range
tn10.tdcs<-tmerge(tn10.tdcs,tn10fortdc.scaled,id=MaskID,BMI = tdc(days,BMI))%>%
  tmerge(tn10fortdc.scaled, id=MaskID, Di.IgG2.RI = tdc(days, Di.IgG2.RI))

nSamplePerSub2<-tn10.tdcs%>%group_by(MaskID)%>%count()
rmMaskID<-nSamplePerSub2[which(nSamplePerSub2$n - nSamplePerSub$n !=0),]$MaskID
tn10fortdc.scaled%<>%filter(!(MaskID%in%rmMaskID & (is.na(Di.IgG2.RI)| days>Time.to.T1D)))

di2.tdcs.dr4<-tn10.tdcs%>%coxph(Surv(tstart, tstop, death) ~ Di.IgG2.RI*DR4+Treatment.Arm+Age+BMI+Sex,data=.)
di2.tdcs.dr4.zp<-cox.zph(di2.tdcs.dr4)
di2.mean<-di2.tdcs.dr4$means["Di.IgG2.RI"]
plot(di2.tdcs.dr4.zp[7]) # looks at Di.IgG2:DR4
plot(di2.tdcs.dr4.zp[1]) # looks at Di.IgG2

# log transform y axis to convert coefficient to hazards ratio
# grab DR4 info from original table
di2.tdcs.dr4.survtime<-data.frame(tstart = di2.tdcs.dr4$y[,1], tstop = di2.tdcs.dr4$y[,2] )
DR4.1<-tn10.tdcs%>%filter(death==1)%>%arrange(futime)%>%select(DR4)
DR4.2<-di2.tdcs.dr4.zp$y[,"DR4"]>0
(DR4.1 == "PRESENT") == DR4.2

data.frame(cbind(di2.tdcs.dr4.zp$y,time = di2.tdcs.dr4.zp$time, DR4.1))%>%
  ggplot(aes(x= time, y = Di.IgG2.RI.DR4,group = DR4.1,color = DR4.1))+
  geom_point(size = 3)+geom_smooth(method="lm",formula =y ~ x,alpha = .3,aes(fill=DR4.1))+
  labs(x="On-Study(days)",y="Coefficient for Di.IgG2-DR4 Interaction")+
  geom_hline(aes(yintercept = 0), linetype = 3,size = 1,color = "#505050") +
  scale_x_continuous(trans=scales::modulus_trans(.5),breaks = seq(0,2500,500))+
  # scale_y_continuous(limits = c(-.1,2))+ # this is for when y = exp(tdcvar), which looks really nice but don't want to exclude outliers
  scale_colour_ghibli_d("MarnieMedium1", direction = -1,label = c("Absent","Present"))+
  scale_fill_ghibli_d("MarnieMedium1", direction = -1,label = c("Absent","Present"))+
  guides(color=guide_legend(title="DR4"),fill=guide_legend(title="DR4"))+
  theme_linedraw(base_size =20)

di2.pterms<-predict(di2.tdcs.dr4,se.fit = T, type="terms") # calculate ßixi for each individual for each term
di2.termplot <- data.frame(di2.pterms[["fit"]],Di.IgG2.se =di2.pterms[["se.fit"]][,1],Di.IgG2.DR4.se = di2.pterms[["se.fit"]][,7] )%>%
  mutate(Di.IgG2.DR4.se = ifelse(DR4==0, Di.IgG2.se,Di.IgG2.DR4.se),
         DR4 = ifelse(DR4==0, "ABSENT","PRESENT"))
di2.termplot <-cbind(di2.termplot,x=tn10fortdc.scaled$Di.IgG2.RI)


di2.ref <- di2.termplot$Di.IgG2.RI.DR4[which.min(abs(di2.termplot$x -di2.mean))]


ggplot(di2.termplot,aes(x=x, group = DR4, color = DR4,fill = DR4)) +
  geom_line(aes(y = exp(Di.IgG2.RI.DR4 - di2.ref)),size = 1.5) +
  geom_line(aes(y = exp(Di.IgG2.RI.DR4-1.96*Di.IgG2.DR4.se - di2.ref )), linetype = 2,size = .1) +
  geom_line(aes(y = exp(Di.IgG2.RI.DR4+1.96*Di.IgG2.DR4.se - di2.ref)), linetype = 2,size = .1) +
  geom_ribbon(aes(ymin=exp(Di.IgG2.RI.DR4-1.96*Di.IgG2.DR4.se - di2.ref ),ymax=exp(Di.IgG2.RI.DR4+1.96*Di.IgG2.DR4.se - di2.ref )),alpha=0.5)+
  geom_hline(aes(yintercept = 1), linetype = 3,size = 1,color = "#505050") +
  geom_vline(aes(xintercept = di2.mean), linetype = 2,size = .5,color = "#505050")+ 
  scale_colour_ghibli_d("MarnieMedium1", direction = -1,label = c("Absent","Present"))+
  scale_fill_ghibli_d("MarnieMedium1", direction = -1,label = c("Absent","Present"))+
  geom_rug()+guides(colour = guide_legend(expression(italic("HLA-DR4"))),fill=guide_legend(expression(italic("HLA-DR4"))))+
  labs(x = expression("All-time"~italic("D.invisus")*"-IgG2"),
       y = "Hazard Ratio"# ,title = "T1D Hazard Ratio as a function of\n Di.IgG2-DR4 interaction"
       )


# Model performance -------------------------------------------------------

tn10selected<-full_join(month6.unscaled,base.unscaled)%>%dplyr::select(MaskID,Time.landmark,T1D.event,DR4,DR3,
                                                                       contains(c("Bl.IgG2.RI","Ef.IgG2.RI","mu_cpep","BMI")),
                                                                       Age,Sex,Treatment.Arm)%>%
  full_join(dplyr::select(month6.scaled,MaskID,Di.IgG2.RI.T6))%>%
  full_join(dplyr::select(base.scaled,MaskID,Di.IgG2.RI.T0))


library(nonnestcox) # https://github.com/thomashielscher/nonnestcox
(coxph.null<-coxph(Surv(Time.landmark,T1D.event)~ 1, data = tn10selected,x=T))
(coxph.trt<-coxph(Surv(Time.landmark,T1D.event)~ Treatment.Arm, data = tn10selected,x=T))
lmtest::lrtest(coxph.null,coxph.trt)
AIC(coxph.null); AIC(coxph.trt)

(coxph.age<-coxph(Surv(Time.landmark,T1D.event)~ Treatment.Arm+Age, data = tn10selected,x=T))
AIC(coxph.age)
plrtest(coxph.age,coxph.trt,nested = T)

(coxph.DR4<-coxph(Surv(Time.landmark,T1D.event) ~Treatment.Arm*DR4+Age,data = tn10selected,x=T))
AIC(coxph.DR4)
lmtest::lrtest(coxph.DR4,coxph.age)
plrtest(coxph.DR4,coxph.age,nested = T)
(coxph.DR4.di0<-coxph(Surv(Time.landmark,T1D.event) ~ Treatment.Arm*DR4+Age+Di.IgG2.RI.T0*DR4,data = tn10selected,x=T))
(coxph.DR4.di6<-coxph(Surv(Time.landmark,T1D.event) ~ Treatment.Arm*DR4+Age+Di.IgG2.RI.T6*DR4,data = tn10selected,x=T))
lmtest::lrtest(coxph.DR4,coxph.DR4.di0) # p < .1
plrtest(coxph.DR4,coxph.DR4.di0,nested = T)
AIC(coxph.DR4.di0)
lmtest::lrtest(coxph.DR4,coxph.DR4.di6) # p < .05

(coxph.cpep0<-coxph(Surv(Time.landmark,T1D.event) ~ mu_cpep.T0+Treatment.Arm+Age+BMI.T0, data = tn10selected,x=T)) # only sign. if account for Age and BMI
AIC(coxph.cpep0)
(coxph.cpep.ef0<-coxph(Surv(Time.landmark,T1D.event) ~ mu_cpep.T0+Treatment.Arm+Age+BMI.T0+Ef.IgG2.RI.T0, data = tn10selected,x=T))
AIC(coxph.cpep.ef0)
lmtest::lrtest(coxph.cpep0,coxph.age)
plrtest(coxph.cpep0,coxph.age,nested = T)
lmtest::lrtest(coxph.cpep0,coxph.cpep.ef0)
plrtest(coxph.cpep0,coxph.cpep.ef0,nested = T)

(coxph.cpep.bl0<-coxph(Surv(Time.landmark,T1D.event) ~ mu_cpep.T0+Treatment.Arm+Age+BMI.T0+Bl.IgG2.RI.T0, data = tn10selected,x=T))
AIC(coxph.cpep.bl0)
lmtest::lrtest(coxph.cpep0,coxph.cpep.bl0)
plrtest(coxph.cpep0,coxph.cpep.bl0,nested = T)

# (coxph.cpep.bl.ef0<-coxph(Surv(Time.to.T1D,T1D.event) ~ mu_cpep.T0+Treatment.Arm+Age+BMI.T0+Bl.IgG2.RI.T0+Ef.IgG2.RI.T0, data = tn10selected))
# (coxph.cpep.bl.ef6<-coxph(Surv(Time.to.T1D,T1D.event) ~ mu_cpep.T0+Treatment.Arm+Age+BMI.T0+Bl.IgG2.RI.T6+Ef.IgG2.RI.T6, data = tn10selected))
# lmtest::lrtest(coxph.cpep.bl.ef0,coxph.cpep.bl0)
# lmtest::lrtest(coxph.cpep.bl.ef0,coxph.cpep0)

(coxph.cpep.tdc<-coxph(Surv(tstart, tstop, death) ~ mu_cpep+Treatment.Arm+Age+BMI+Sex,data=tn10.tdcu))
(coxph.cpep.ef.tdc<-coxph(Surv(tstart, tstop, death) ~ mu_cpep+Treatment.Arm+Age+BMI+Sex+Ef.IgG2.RI,data=tn10.tdcu))
lmtest::lrtest(coxph.cpep.tdc,coxph.cpep.ef.tdc) # p< .1
(coxph.cpep.bl.tdc<-coxph(Surv(tstart, tstop, death) ~ mu_cpep+Treatment.Arm+Age+BMI+Sex+Bl.IgG2.RI,data=tn10.tdcu))
lmtest::lrtest(coxph.cpep.bl0,coxph.cpep.bl.tdc) #p <.05

(coxph.di.tdc<-coxph(Surv(tstart, tstop, death) ~ Di.IgG2.RI*DR4+Treatment.Arm*DR4+Age,data=tn10.tdcs))
lmtest::lrtest(coxph.di.tdc,coxph.DR4.di0)# not as good as single time point

# ROC curves for survival model
# https://rpubs.com/kaz_yos/survival-auc
## Put linear predictors ("lp") into pbc dataset
tn10selected$model1<-predict(coxph.trt, type = "lp")
tn10selected$model2<-predict(coxph.age, type = "lp")
tn10selected$model3<-predict(coxph.DR4, type = "lp")

tn10selected$model4<-predict(coxph.DR4.di0, type = "lp")

tn10selected$model5<-predict(coxph.cpep0, type = "lp")

tn10selected$model6<-predict(coxph.cpep.ef0, type = "lp")

tn10selected$model7<-predict(coxph.cpep.bl0, type = "lp")


library(risksetROC) # this one is smoothed and gives better result
# uses MASS so don't import it from the begining

multiPlotRisksetROC <-function(time,df,colSelect,colTime,colStatus,baseSize,labSize,ref = T,orderLevels=NA,orderLabels=waiver()){
  res.base <- apply(df%>%dplyr::select(all_of(colSelect)),2, function(x)
    risksetROC(Stime        = df[,colTime],
               status       = df[,colStatus],
               marker       = x,
               predict.time = time,
               plot         = FALSE))
  output<-vector(mode = "list", length = 2)
  toplot<-res.base%>%imap(~data.frame(TruePositive=.$TP,FalsePositive=.$FP,AUC=.$AUC,modelname = .y))%>%reduce(rbind)
  if(ref){
    refdata <- data.frame(TruePositive=seq(0,1,1/100),FalsePositive=seq(0,1,1/100),AUC = rep(.5,101), modelname = rep("reference",101))
    toplot<-rbind(toplot,refdata)
  }
  toplot$modelname<-str_remove(toplot$modelname,"lp\\.")
  if(is.na(orderLevels)){
    toplot$modelname<-factor(toplot$modelname)
  }else{
    toplot$modelname<-factor(toplot$modelname,levels = orderLevels)
  }
  output[[1]]<-toplot
  
  p<-ggplot(toplot,aes(x=FalsePositive,y=TruePositive,group = modelname,color=modelname,linetype = modelname),show.legend = FALSE)+
    geom_smooth(se = F)+
    labs(x = "False Positive Rate", y = "True Positive Rate")
  output[[2]]<-p
  
  return(output)
}

# orderLevels<-c("reference","age","trt","DR4","DR4.di0","cpep0","cpep.ef0","cpep.bl0")
orderLevels<-c("reference",paste0("model",1:7))

dr4models<-multiPlotRisksetROC(550,tn10selected,paste0("model",1:4),"Time.landmark","T1D.event", 15,4,T,orderLevels)
dr4models[[1]]
dr4models[[2]]+scale_color_manual(name = "",values = c("grey",ghibli_palette("SpiritedMedium")[6:3]))+scale_linetype_manual(name = "",values = c("dashed",rep("solid",4)))

cpepmodels<-multiPlotRisksetROC(550,tn10selected,paste0("model",c(1,2,5:7)),"Time.landmark","T1D.event", 15,4,T,orderLevels)
cpepmodels[[1]]
cpepmodels[[2]]+scale_color_manual(name = "",values = c("grey",ghibli_palette("SpiritedMedium")[6:5],ghibli_palette("MarnieMedium2")[5:1]))+scale_linetype_manual(name = "",values = c("dashed",rep("solid",5)))

multiPlotRisksetAUC <-function(time,df,colSelect,colTime,colStatus,cutTime,baseSize=30){
  res.base <- apply(df%>%dplyr::select(all_of(colSelect)),2, function(x)
    risksetAUC(Stime        = df[,colTime],
               status       = df[,colStatus],
               marker       = x,
               tmax = time,
               plot         = FALSE))
  output<-vector(mode = "list", length = 2)
  toplot<-res.base%>%imap(~data.frame(Time=.$utimes,AUC=.$AUC,Cind =.$Cindex,modelname = .y))%>%reduce(rbind)
  
  toplot$modelname<-str_remove(toplot$modelname,"lp\\.")
  output[[1]]<-toplot
  
  p<-ggplot(toplot,aes(x=Time,y=AUC,group = modelname,color=modelname),show.legend = FALSE)+
    geom_vline(xintercept = cutTime,size = 1,linetype = "dashed",color = "grey4")+
    geom_line(size = 1.2)+labs(x = "Time (Days)")+
    # geom_label_repel(data = toplot%>%group_by(modelname)%>%arrange(abs(FalsePositive - median(FalsePositive))) %>% slice(1),
    #                  size = labSize,aes(label = paste(modelname,round(AUC,2),sep = ":")))
    scale_x_continuous(limits = c(0,time))+scale_y_continuous(limits = c(0.5,.8))+theme_classic(base_size = baseSize)
  output[[2]]<-p
  
  return(output)
}
dr4models<-multiPlotRisksetAUC(1500,tn10selected,paste0("model",1:4),"Time.landmark","T1D.event", 550,30)
dr4models[[1]]
dr4models[[2]]+scale_color_manual(name = "",values = c(ghibli_palette("SpiritedMedium")[6:3]))
rcorrcens(formula = Surv ~ I(-1 * model7), data = tn10selected)
cpepmodels<-multiPlotRisksetAUC(1500,tn10selected,paste0("model",c(1,2,5:7)),"Time.landmark","T1D.event",550,30)
cpepmodels[[1]]
cpepmodels[[2]]+scale_color_manual(name = "",values = c(ghibli_palette("SpiritedMedium")[6:5],ghibli_palette("MarnieMedium2")[5:1]))
