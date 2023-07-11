library(tidyverse)
library(UpSetR)

tn10<-read.csv("TN10full_smallgate.csv", header = T, as.is = T)
str(tn10)
theme_set(
  theme_bw(base_size = 20)
)
tn10%<>%fill(T1D.event,Time.to.T1D,Age,Gender,Treatment.Arm)%>% # for 3 samples w/o ACAb info but has Cpep and other metabolic measures
  mutate(Sex = factor(Gender))%>%select(-Gender)%>%
  mutate(Date_of_Draw = as.Date(Date_of_Draw),
         Date.IDDM = as.Date(Date.IDDM),
         DR3 = factor(DR3, levels = c("ABSENT","PRESENT")),
         DR4 = factor(DR4, levels = c("ABSENT","PRESENT")),
         Sex = ifelse(Sex == "Male", -1, 1),
         Treatment.Arm = factor(Treatment.Arm), 
         GADA = ifelse(GADA == 0, -1, 1),
         IA.2A = ifelse(IA.2A == 0, -1, 1),
         mIAA = ifelse(mIAA == 0, -1, 1),
         ZnT8A = ifelse(ZnT8A == 0, -1, 1),
         glycemia = factor(glycemia, levels = c("normal","dysglycemia","hyperglycemia")))%>%
  filter(!is.na(Barcode))

test<-tn10%>%select(MaskID,Visit)%>%distinct(MaskID,Visit, .keep_all = T)%>%mutate(bin = 0)%>%
  pivot_wider(id_cols = MaskID,names_from = Visit,values_from = bin,values_fill=list(bin=1))
invertmiss<-apply(test,2,function(x)ifelse(x==0,NA,x))%>%data.frame(.)
colnames(invertmiss)[2:5]<-c("Baseline", "6 months", "12 months", "18 months")
gg_miss_upset(invertmiss)

test2<-tn10%>%select(MaskID,Visit)%>%distinct(MaskID,Visit, .keep_all = T)%>%mutate(bin = 1)%>%
  pivot_wider(names_from = Visit,values_from = bin)
test2ed<-apply(test2[2:5],2,function(x) x*test2$MaskID)%>%data.frame()
colnames(test2ed)<-c("Baseline", "6 months", "12 months", "18 months")
test2list <- apply(test2ed, 2, as.list)
upset(fromList(test2list), order.by = "freq",mainbar.y.label = "No. of individuals",sets.x.label = "Timepoint", text.scale = 2)
