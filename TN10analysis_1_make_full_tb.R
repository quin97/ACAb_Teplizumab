library(readxl)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(DescTools)
library(lubridate)
library(magrittr)
library(readxl)
# tn10ri<-read.csv("TN10RI.csv")
# str(tn10ri)
# tn10ri.small<-select(tn10ri,colnames(tn10ri)[!grepl("large",colnames(tn10ri))])
tn10ri.small<-read_xlsx("../TN10RIdata/Data file S1.xlsx",sheet = "RI_discovery")

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X)) # guess_max = 60
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
tn10meta<-read_excel_allsheets("TN10 Clinical Dataset_Danksa_20210827.xlsx")
#============================ Serum List============================
serum.list<-tn10meta$`Serum List`%>%select(-`Random ID`)
# 2 repeated samples @ baseline: 420795,886154
serum.list%>% filter(Visit == "Baseline")%>%
  count(MaskID,Visit) %>% 
  filter(n > 1)
serum.list$Date_of_Draw<-as.Date(serum.list$Date_of_Draw)
serum.list$Visit<-gsub("Baseline",0,serum.list$Visit)
serum.list$Visit<-parse_number(serum.list$Visit)

tn10all.small<-merge(tn10ri.small, serum.list, by = "Barcode")%>%arrange(MaskID,Visit)
#============================ Treatment Outcome ============================
tx.outcome<-tn10meta$Dem_Tx_Outcome
str(tx.outcome)
# Date.Followup and Date.IDDM are in Excel date format (days since 1900-01-01) 
# (Unix Epoch is from 1970-01-01)
tx.outcome$Date.Followup<-parse_number(tx.outcome$Date.Followup)
tx.outcome$Date.Followup<-as.Date(tx.outcome$Date.Followup,origin="1900-01-01")
tx.outcome$Date.IDDM<-parse_number(tx.outcome$Date.IDDM)
tx.outcome$Date.IDDM<-as.Date(tx.outcome$Date.IDDM,origin="1900-01-01")

tx.outcome$`Treatment Start Date`<-as.Date(tx.outcome$`Treatment Start Date`)
cutdate<-mdy("03/19/2020")

sum(tx.outcome$T1D.event == 1)
tx.outcome%>% group_by(T1D.event,`Treatment Arm`)%>%summarise(count = n())
summary(tx.outcome[tx.outcome$`Treatment Arm`=="Teplizumab",]$Age)
summary(tx.outcome[tx.outcome$`Treatment Arm`=="Placebo",]$Age)
dim(tx.outcome[tx.outcome$`Treatment Arm`=="Teplizumab" & tx.outcome$Age <= 18,])[1]
dim(tx.outcome[tx.outcome$`Treatment Arm`=="Placebo" & tx.outcome$Age <= 18,])[1]
dim(tx.outcome[tx.outcome$`Treatment Arm`=="Teplizumab" & tx.outcome$Gender=="Female",])[1]
dim(tx.outcome[tx.outcome$`Treatment Arm`=="Placebo" & tx.outcome$Gender=="Female",])[1]

# check if time to T1D is the same as date difference between diagnosis date and treatment start date
tx.outcome<-tx.outcome%>%mutate(cutoff= pmin(Date.IDDM,cutdate) , date.diff = (cutoff -`Treatment Start Date`))
tx.outcome$date.diff<-as.numeric(tx.outcome$date.diff)
tx.outcome%<>%select(-`Treatment Start Date`,-Date.Followup,-cutoff,-date.diff)
tn10all.small<-tx.outcome%>%left_join(tn10all.small, by = "MaskID")

#============================ BMI ============================
height_weight<-tn10meta$`Ht-Wt`[tn10meta$`Ht-Wt`$Visit %in% c("Baseline","6 months", "12 months","18 months"),]
# need to convert here to get NA
height_weight$HeightCM<-parse_number(height_weight$HeightCM)
height_weight$Impute_Height<-parse_number(height_weight$Impute_Height)
height_weight$HeightCM2<-height_weight$HeightCM
# for those without height, need to use imputed height
height_weight$HeightCM2[is.na(height_weight$HeightCM2)]<-height_weight$Impute_Height[is.na(height_weight$HeightCM2)]
height_weight$Visit<-gsub("Baseline",0,height_weight$Visit)
height_weight$Visit<-parse_number(height_weight$Visit)
height_weight$WeightKG<-parse_number(height_weight$WeightKG)
height_weight<-height_weight%>%mutate(BMI=WeightKG/HeightCM2^2*10000)%>%select(MaskID,Visit, BMI)
sum(is.na(height_weight$BMI))
tn10all.small<-left_join(tn10all.small,height_weight, by = c("MaskID","Visit"))%>%
  select(MaskID,Visit,BMI,Barcode,ID,Date_of_Draw,everything())%>%
  arrange(MaskID,Visit)%>%fill(BMI)
#============================ Autoantibodies ============================
# KH's table
abtable<-read_excel("TN10_autoabdata_K_Herold_08_2021.xlsx", skip=1)
abtable<-abtable%>%select(MaskID:TGA)%>%filter(!(MaskID%in%c(839473,983719)))%>%arrange(MaskID, `Draw Date`)
str(abtable)
abtable$`Draw Date`<-as.Date(abtable$`Draw Date`)
abtable[,5:9]<-apply(abtable[,5:9], 2, as.numeric)
apply(abtable[,5:9], 2, summary)
# convert autoAB data to binary
abtable<-abtable%>%mutate(GADA = ifelse(GADA > 20, 1, 0), 
                          `IA-2A` = ifelse(`IA-2A` > 5, 1, 0), 
                          mIAA = ifelse(mIAA > .01, 1, 0), 
                          ZnT8A=ifelse(ZnT8A > .02, 1, 0),
                          TGA=ifelse(TGA > .05, 1, 0))
# count a person as 1 if they're seroconverted at any timepoint
# maybe have another way to measure seroconversion: e.g.GADA = mean(GADA) > .3
ab.alltime<-abtable%>%group_by(MaskID)%>%summarise(GADA = max(GADA), 
                                                   `IA-2A` = max(`IA-2A`), 
                                                   mIAA = max(mIAA), 
                                                   ZnT8A=max(ZnT8A))

# trialnet table
abtable2<-tn10meta$Autoantibodies%>%arrange(MaskID,Date_of_Draw)
abtable2$Date_of_Draw<-as.Date(abtable2$Date_of_Draw)
abtable2[,3:9]<-apply(abtable2[,3:9], 2, as.numeric)
str(abtable2)
apply(abtable2[,3:9], 2, summary)
# convert autoAB data to binary
abtable2<-abtable2%>%mutate(GAD65 = ifelse(GAD65 > .032, 1, 0),
                            GAD65H = ifelse(GAD65H > 20, 1, 0), 
                            ICA512=ifelse(ICA512 > 0.049, 1, 0),
                            ICA=ifelse(ICA >= 10, 1, 0),
                            `IA-2H` = ifelse(`IA-2H` > 5, 1, 0),
                            MIAA = ifelse(MIAA > .01, 1, 0),
                            ZNT8=ifelse(ZNT8 > .02, 1, 0))
abtable2$GAD65H2<-abtable2$GAD65H
abtable2$GAD65H2[is.na(abtable2$GAD65H2)]<-abtable2$GAD65[is.na(abtable2$GAD65H2)]


# count a person as 1 if they're seroconverted at any timepoint
ab.alltime2<-abtable2%>%group_by(MaskID)%>%summarise(GAD65H = max(GAD65H2,na.rm = T), 
                                                     `IA-2H` = max(`IA-2H`,na.rm = T), 
                                                     MIAA = max(MIAA,na.rm = T), 
                                                     ZNT8=max(ZNT8,na.rm = T))

ab.summ<-left_join(ab.alltime,ab.alltime2)
# write.csv(ab.summ,"autoAB_alltime_merged.csv",row.names = F)
# data looks largely agree, except:
# 251760,307044,308094,341455,425198,
# 441014,455189,488675,494195,584001,
# 661064,666821,697933,716379,787396,
# 798801,802857,817894,838347,846511,
# 873459,886154,926463,971999

# create a table combining all timepoints from the table
# TN table all have earlier timepoints than KH
ab.all<-abtable%>%select(-Date_Rcvd, -TGA)%>%
  rename(Date_of_Draw=`Draw Date`)
join_by_DOD<-abtable2%>%select(MaskID,Date_of_Draw,GAD65H2,`IA-2H`,MIAA,ZNT8)%>%full_join(ab.all,by = c("MaskID","Date_of_Draw"))%>%arrange(MaskID,Date_of_Draw)  
# there are no overlap between Date_of_Draw in these 2 tables
# TrialNet mostly have timepoints before TN10, KH have timepoints after TN10
# write.csv(join_by_DOD,"autoAB_comp_alltime.csv", row.names = F,na=".")
# people with additional AAb in trial
# 232794,245390,369952,453691,471457,
# 488675,557085,584001,690566,872759
# people with less AAb in trial
# 251760 (+mIAA, -GADA),268041 (-ZnT8A), 
# 437809 (-mIAA),455189 (-IA2A),494195(GADA),817894(-IA2A)


ab.all2<-abtable2%>%select(MaskID,Date_of_Draw,GAD65H2,`IA-2H`,MIAA,ZNT8)%>%
  rename(GADA=GAD65H2,
         `IA-2A`=`IA-2H`,
         mIAA=MIAA,
         ZnT8A=ZNT8)

ab.all2$Visit<-NA
ab.all3<-rbind(ab.all,ab.all2)%>%arrange(MaskID,Date_of_Draw)
# create a column "baseline" containing baseline for each MaskID, matched from table ab.all
# https://stackoverflow.com/questions/25539326/filling-in-columns-with-matching-ids-from-two-dataframes-in-r
ab.all3$baseline<-ab.all[ab.all$Visit=="Baseline",]$Date_of_Draw[match(ab.all3$MaskID,ab.all[ab.all$Visit=="Baseline",]$MaskID)]
# calculate month of visit for each sample based on baseline and Date_of_Draw
ab.all3<-ab.all3%>%mutate(Visit2 = baseline%--%Date_of_Draw%/%months(1))

# calc_months<-function(ref_date,current_date){
#   ref_date%--%current_date/months(1)
# }

ab.all3[is.na(ab.all3$Visit),]$Visit<-paste(ab.all3[is.na(ab.all3$Visit),]$Visit2, "months")
# write.csv(ab.all3,"autoAB_alltime.csv",row.names = F,na=".")
abtableForTN10<-ab.all3%>%filter(Visit %in% c("Baseline","6 months", "12 months","18 months"))%>%select(MaskID,Visit,GADA,`IA-2A`,mIAA,ZnT8A)
abtableForTN10$Visit<-gsub("Baseline",0,abtableForTN10$Visit)
abtableForTN10$Visit<-parse_number(abtableForTN10$Visit)
tn10all.small<-abtableForTN10%>%full_join(tn10all.small, by = c("MaskID","Visit"))%>%arrange(MaskID, Visit)

dim(tn10all.small[tn10all.small$`Treatment Arm`=="Teplizumab" & tn10all.small$Visit == 0 & tn10all.small$GADA == 1,])[1]
dim(tn10all.small[tn10all.small$`Treatment Arm`=="Placebo" & tn10all.small$Visit == 0 & tn10all.small$GADA == 1,])[1]
dim(tn10all.small[tn10all.small$`Treatment Arm`=="Teplizumab" & tn10all.small$Visit == 0 & tn10all.small$mIAA == 1,])[1]
dim(tn10all.small[tn10all.small$`Treatment Arm`=="Placebo" & tn10all.small$Visit == 0 & tn10all.small$mIAA == 1,])[1]
dim(tn10all.small[tn10all.small$`Treatment Arm`=="Teplizumab" & tn10all.small$Visit == 0 & tn10all.small$`IA-2A` == 1,])[1]
dim(tn10all.small[tn10all.small$`Treatment Arm`=="Placebo" & tn10all.small$Visit == 0 & tn10all.small$`IA-2A` == 1,])[1]
dim(tn10all.small[tn10all.small$`Treatment Arm`=="Teplizumab" & tn10all.small$Visit == 0 & tn10all.small$ZnT8A == 1,])[1]
dim(tn10all.small[tn10all.small$`Treatment Arm`=="Placebo" & tn10all.small$Visit == 0 & tn10all.small$ZnT8A == 1,])[1]

# count a person as 1 if they're seroconverted at least 2X 
ab.all4<-ab.all3%>%group_by(MaskID)%>%
  summarise(GADA = mean(GADA,na.rm = T),
            `IA-2A` = mean(`IA-2A`,na.rm = T),
            mIAA = mean(mIAA,na.rm = T),
           ZnT8A=mean(ZnT8A,na.rm = T))
ab.all5<-cbind(ab.all4$MaskID,apply(as.data.frame(ab.all4[,2:5] >.3),2,as.integer))
# write.csv(ab.all5,"autoAB_alltime_merged.csv",row.names = F,na=".")

non_overlap_autoab<-anti_join(abtable,abtable2) 
# compare autoAB data between trialnet's and KH's, KH has 2 extra individuals (839473,983719),
# KH table has Visit which is easier to match with other metadata, 
# KH table and trialnet table very different autoAB data dk what to use

#============================ OGTT related indices and Hba1c ============================
# convert OGTT glucose to dysglycemia/hyperglycemia indicator
ogtt<-tn10meta$OGTT%>%select(-INSM10,-INS0,-INS30,-INS60,-INS90,-INS120)
str(ogtt)
ogtt[,4:28]<-apply(ogtt[,4:28],2,parse_number)
ogtt$Date_of_Draw<-as.Date(ogtt$Date_of_Draw)
ogtt$GLU30to90 <- rowSums(ogtt[,c("GLU30","GLU60","GLU90")]>200)
ogtt<-left_join(ogtt, select(tx.outcome,MaskID,Date.IDDM))
# actually in this cohort there're no followups after T1D, so hyperglycemia should be changed to dysglycemia

# ogtt<-ogtt%>%mutate(glycemia=case_when(
#   GLUM10 > 125 | GLU120 >= 200 ~ "hyperglycemia",
#   between(GLUM10, 110, 125)| between(GLU120, 140, 199)| GLU30to90 >=2 ~ "dysglycemia",
#   TRUE ~ "normal")) # order matters for case_when

# If someone has hyperglycemia would've been diagnosed.. just put dysglycemia for everything then
ogtt<-ogtt%>%mutate(glycemia=case_when(
  Date_of_Draw > Date.IDDM ~ "hyperglycemia",
  GLUM10 >= 110 | GLU120 >= 140 | GLU30to90 >=2 ~ "dysglycemia",
  TRUE ~ "normal")) # order matters for case_when
table(ogtt$glycemia)

# calculate AUC for a bunch of ogtt indices
ogtt$glu_AUC<-apply(ogtt[,c("GLU0","GLU30","GLU60","GLU90","GLU120")],1,function(x)AUC(c(0,30,60,90,120),x)/120)
summary(ogtt$glu_AUC)
ogtt$ins_tot<-apply(ogtt[,c("INST0","INST30","INST60","INST90","INST120")],1,function(x)AUC(c(0,30,60,90,120),x)/120)
summary(ogtt$ins_tot)
ogtt$ins_1h<-apply(ogtt[,c("INST0","INST30","INST60")],1,function(x)AUC(c(0,30,60),x)/60)
summary(ogtt$ins_1h)
ogtt$ins_2h<-apply(ogtt[,c("INST60","INST90","INST120")],1,function(x)AUC(c(60,90,120),x)/60)
summary(ogtt$ins_2h)
ogtt<-ogtt%>%select(MaskID,Visit,glycemia,glu_AUC,mu_cpep,ins_tot,ins_1h,ins_2h,everything())%>%select(-Date.IDDM)
# # confirm how C-pep AUC mean is calculated
# ogttcp<-ogtt%>%select("PEPM10_pmol/mL","PEP0_pmol/mL" ,"PEP30_pmol/mL","PEP60_pmol/mL","PEP90_pmol/mL","PEP120_pmol/mL","mu_cpep")
# ogttcp$cpAUC<-apply(ogttcp[,2:6],1,function(x)AUC(c(0,30,60,90,120),x)/120) # cpAUC should be the same as mu_cpep

# # Baseline date check between sheets --------------------------------------
# # everyone has baseline date from the sheet "Lymphocyte" and "Serum List", maybe use this baseline
# serum.list2<-tn10meta$`Serum List`
# baseserum<-serum.list2%>%filter(Visit == "Baseline")%>%filter(!(Barcode %in% c("TN90GTM","TN069JH")))
# lymphocytes<-tn10meta$Lymphocyte
# baselymph<-lymphocytes[lymphocytes$Visit == "Baseline",]
# baseserumlymph<-baselymph%>%select(MaskID:Date_of_Visit)%>%full_join(baseserum)
# # baseline from "Lymphocyte" and "Serum List" are the same, except 378864 and 700956 (1 day diff)
# 
# tcell<-tn10meta$`T Cell`%>%select(-Comments)
# baset<-tcell[tcell$Visit == "Baseline",]
# baseserumt<-baset%>%full_join(baseserum,by = c("MaskID", "Visit"))
# # baseline from "T cell" and "Serum List" are the same!
# 
# baseab<-ab.all[ab.all$Visit=="Baseline",]
# baseserumab<-baseab%>%full_join(baseserum,by = c("MaskID", "Visit"))
# # baseline from "Autoantibodies" and "Serum List" are the same!

#============================ Hba1c ============================
# Date_of_Draw from "OGTT" and "HbA1c" are the same
hba1c<-tn10meta$HbA1c%>%select(-PID) 
ogtt<-left_join(ogtt,hba1c)

month6<-ogtt[ogtt$Visit=="6 months",] # "6 months" is the time point where everyone has available data
withbase<-ogtt[ogtt$Visit=="Baseline",]
withscreen<-ogtt[ogtt$Visit=="Screening",]
intersect(withbase$MaskID,withscreen$MaskID)
# only one sample has both Screening and Baseline: 494195

# create a column of baseline
# impute baseline date based on 6 months visit, bacause everyone has 6 months visit
ogtt$baseline_6m<-as_date(ogtt[ogtt$Visit=="6 months",]$Date_of_Draw[match(ogtt$MaskID,ogtt[ogtt$Visit=="6 months",]$MaskID)] - months(6))
# for 16 of the individuals they have baseline date so just use baseline date
ogtt$baseline<-as_date(ogtt[ogtt$Visit=="Baseline",]$Date_of_Draw[match(ogtt$MaskID,ogtt[ogtt$Visit=="Baseline",]$MaskID)])
# for those w/o baseline date, use 6-months-imputed-baseline instead
ogtt[is.na(ogtt$baseline),"baseline"]<-ogtt[is.na(ogtt$baseline),"baseline_6m"]
# calculate month of visit for each sample based on baseline and Date_of_Draw
ogtt<-ogtt%>%mutate(Visit2 = baseline%--%Date_of_Draw%/%months(1))

# get baseline date from "Lymphocyte", b/c everyone has baseline in that sheet
lymphocytes<-tn10meta$Lymphocyte
ogtt$baseline_lymph<-as_date(lymphocytes[lymphocytes$Visit=="Baseline",]$Date_of_Visit[match(ogtt$MaskID,lymphocytes[lymphocytes$Visit=="Baseline",]$MaskID)])
# check if lymphocyte baseline are the same as ogtt baseline
testbaseline<-ogtt[ogtt$Visit=="Baseline",]
testbaseline[testbaseline$baseline_lymph != testbaseline$Date_of_Draw,]
# 1 individual(656446) has different baseline from ogtt and Lymphocyte
# for this one use date from Visit=="Baseline"(i.e. ogtt baseline)
ogtt<-ogtt%>%mutate(Visit_lymph = baseline_lymph%--%Date_of_Draw%/%months(1))

# rearrange column order for easier inspection
ogtt%<>%select(baseline,baseline_6m,baseline_lymph,Date_of_Draw,Visit,Visit2,Visit_lymph,everything())

ogtt<-ogtt%>%mutate(diff = Visit2-Visit_lymph)
ogtt[abs(ogtt$diff)>1,] 
# only 1 individual(690566) has 2 month diff between baseline_6m and baseline_lymph
# for this one use date from Visit=="Screening" 

ogtt$baseline_lymph[ogtt$MaskID == 656446]<-ogtt$Date_of_Draw[ogtt$MaskID == 656446 & ogtt$Visit=="Baseline"]
ogtt$baseline_lymph[ogtt$MaskID == 690566]<-ogtt$Date_of_Draw[ogtt$MaskID == 690566 & ogtt$Visit=="Screening"]

ogtt<-ogtt%>%select(-diff,-baseline,-baseline_6m,-Visit_lymph)%>%
  rename(baseline = "baseline_lymph")%>%
  mutate(Visit2 = baseline%--%Date_of_Draw%/%months(1))

# select samples for timepoints available in tn10meta  
# find samples with closest date to baseline
baselineogtt<-ogtt%>%group_by(MaskID)%>%slice_min(abs(Visit2))
# slice(which.min(abs(Visit2))) only gives the first occurrence
# 558054 has 2 samples at time 0: TN01 and Screening
# 610265 has 2 samples at time 0: TN01 and Baseline
# 672212 has 2 samples at time 0: PRN and Baseline
baselineogtt%<>%filter(!(MaskID%in%c(558054,610265,672212) & Visit%in%c("TN01","PRN")))

# make a new column of Visit for TN10 metatable
# set Visits at these Date_of_Draw as Baseline
ogtt$Visit_new<-ogtt$Visit
ogtt$Visit_new[!grepl("month",ogtt$Visit)]<-NA
ogtt[ogtt$Visit2%in%(-1:0) & 
       ogtt$Date_of_Draw%in%baselineogtt$Date_of_Draw & 
       !(ogtt$MaskID%in%c(251760,954069) & ogtt$Visit2 == -1),"Visit_new"]<-"Baseline"

# these samples have usable Visits not marked in "XX months", need to rename
drawdates<-as.Date(c("2019-02-05","2014-02-12","2016-01-12","2014-11-11"))
# maskIDlist<-c(666617,798801,873459)
ogtt[ogtt$Date_of_Draw %in% drawdates,"Visit_new"]<-"18 months"
ogtt[ogtt$Date_of_Draw == "2018-08-14","Visit_new"]<-"12 months"
ogtt<-distinct(ogtt)

ogttForTN10Meta<-ogtt%>%filter(Visit_new %in% c("Baseline","6 months","12 months","18 months"))%>%select(-Visit)%>%rename(Visit="Visit_new")
ogttForTN10Meta%<>%select(MaskID, Visit,Visit2,baseline,Date_of_Draw,everything())
# # remove 3 samples that are not needed 
# ogttForTN10Meta<-ogttForTN10Meta[-which(ogttForTN10Meta$MaskID == "494195" & ogttForTN10Meta$Visit == "Screening"),] # this one already have "Baseline"
# ogttForTN10Meta<-ogttForTN10Meta[-which(ogttForTN10Meta$MaskID == "838347" & ogttForTN10Meta$Visit == "Screening"),] # this one already have "Baseline"
# ogttForTN10Meta<-ogttForTN10Meta[-which(ogttForTN10Meta$MaskID == "420795" & ogttForTN10Meta$Visit2 == "18"),] # this one already have "18 months"

ogttForTN10Meta$Visit<-gsub("Baseline",0,ogttForTN10Meta$Visit)
ogttForTN10Meta$Visit<-parse_number(ogttForTN10Meta$Visit)


# check primary keys
ogttForTN10Meta %>% 
  count(MaskID,Visit) %>% 
  filter(n > 1)
ogttForTN10Meta<-ogttForTN10Meta%>%select(MaskID,Visit,glycemia,glu_AUC,mu_cpep,ins_tot,ins_1h,ins_2h,HbA1c)
tn10all.small<-left_join(tn10all.small,ogttForTN10Meta, by = c("MaskID","Visit"))


summary(tn10all.small[tn10all.small$`Treatment Arm`=="Teplizumab" & tn10all.small$Visit == 0,]$mu_cpep)
summary(tn10all.small[tn10all.small$`Treatment Arm`=="Placebo" & tn10all.small$Visit == 0,]$mu_cpep)

#============================ Lymphocyte ============================
lymphocytes<-tn10meta$Lymphocyte%>%filter(Visit %in% c("Baseline","6 months", "12 months","18 months"))%>%select(-Date_of_Visit)
lymphocytes$Visit<-gsub("Baseline",0,lymphocytes$Visit)
lymphocytes<-data.frame(apply(lymphocytes,2,parse_number))
lymphocytes<-distinct(lymphocytes)
tn10all.small<-left_join(tn10all.small,lymphocytes)
#============================ HLA ============================
hla<-tn10meta$HLA%>%select(MaskID,DR3,DR4)
tn10all.small<-left_join(tn10all.small,hla)
dim(tn10all.small[tn10all.small$`Treatment Arm`=="Teplizumab" & tn10all.small$Visit == 0 & tn10all.small$DR3 == "PRESENT",])[1]
dim(tn10all.small[tn10all.small$`Treatment Arm`=="Teplizumab" & tn10all.small$Visit == 0 & tn10all.small$DR4 == "PRESENT",])[1]
dim(tn10all.small[tn10all.small$`Treatment Arm`=="Placebo" & tn10all.small$Visit == 0 & tn10all.small$DR3 == "PRESENT",])[1]
dim(tn10all.small[tn10all.small$`Treatment Arm`=="Placebo" & tn10all.small$Visit == 0 & tn10all.small$DR4 == "PRESENT",])[1]

#============================ T cells ============================
tcell<-tn10meta$`T Cell`%>%select(-Comments)
# 798801 and 902103 has unscheduled visit
tcell[tcell$Visit=="PRN","Visit"]<-c(18,3)
tcell<-select(tcell,-Date_of_Draw)
tcell$Visit<-gsub("Baseline",0,tcell$Visit)
tcell<-data.frame(apply(tcell,2,parse_number))
tn10all.small<-left_join(tn10all.small,tcell)

write_csv(tn10all.small,"TN10full_smallgate.csv")


#============================ more timepoint, no RI ============================
lymphocytes2<-tn10meta$Lymphocyte%>%filter(!(grepl("Day",Visit)|grepl("weeks",Visit)))%>%rename(Date_of_Draw = "Date_of_Visit")
# 585 rows w/o Day, 522 rows w/o Day and weeks
lymphocytes2[,4:6]<-apply(lymphocytes2[,4:6],2,parse_number)
lymphocytes2$baseline<-as.Date(ogtt[ogtt$Visit_new=="Baseline",]$baseline[match(lymphocytes2$MaskID,ogtt[ogtt$Visit_new=="Baseline",]$MaskID)])
lymphocytes2<-lymphocytes2[!is.na(lymphocytes2$`White Blood Cell Count`) |!is.na(lymphocytes2$`Lymphocytes_%`)|!is.na(lymphocytes2$Lymphocytes_ABS), ]
lymphocytes2%<>%distinct()%>%mutate(Visit2 = baseline%--%Date_of_Draw%/%months(1))

tcell2<-tn10meta$`T Cell`%>%select(-Comments)
tcell2$`%CD8+TIGIT+KLRG1+ of Total CD3`<- parse_number(tcell2$`%CD8+TIGIT+KLRG1+ of Total CD3`)
tcell2$baseline<-as.Date(ogtt[ogtt$Visit_new=="Baseline",]$baseline[match(tcell2$MaskID,ogtt[ogtt$Visit_new=="Baseline",]$MaskID)])
tcell2%<>%drop_na(`%CD8+TIGIT+KLRG1+ of Total CD3`)%>%distinct()%>%mutate(Visit2 = baseline%--%Date_of_Draw%/%months(1))

lymph.t<-full_join(lymphocytes2,tcell2, by = c("MaskID", "Visit", "baseline", "Visit2"))%>%select(-Date_of_Draw.y)%>%rename(Date_of_Draw="Date_of_Draw.x")

ogtt1<-ogtt[!duplicated(ogtt[,c("MaskID","Date_of_Draw","mu_cpep","HbA1c")]),]
# 471457 2013-12-18 seems duplicated but not removed
ogtt1[is.na(ogtt1$Visit_new),"Visit_new"]<-"***"
ogtt.all2<-full_join(ogtt1, lymph.t,by = c("MaskID","baseline", "Visit2","Date_of_Draw"))
ogtt.all2%>%count(Date_of_Draw,MaskID,mu_cpep) %>% filter(n > 1) # 3 samples w/ 2 lines of lymph data for same Date_of_Draw: 802857,369952,425198
dim(ogtt.all2)
ogtt.all2%<>%select(MaskID,Date_of_Draw,baseline,Visit.x,Visit.y,Visit2,`White Blood Cell Count`:`%CD8+TIGIT+KLRG1+ of Total CD3`,everything())
sum(is.na(ogtt.all2$`White Blood Cell Count`))
sum(is.na(ogtt.all2$Visit.y))

lympht_only<-ogtt.all2[is.na(ogtt.all2$Visit_new) ,]# these timepoints dont ogtt for same Date_of_Draw
lympht_only%>%count(Date_of_Draw,MaskID) %>% filter(n > 1) 
# of those not matched with ogtt,  2 samples (802857,369952) 2 lines of different lymph data for same draw date
# merge them to one line, so last command returns <0 rows>
lympht_only[which(lympht_only$MaskID == 802857 & lympht_only$Date_of_Draw == ymd("2012-10-29") & is.na(lympht_only$`Lymphocytes_%`)),]$`Lymphocytes_%`<-lympht_only[which(lympht_only$MaskID == 802857 & lympht_only$Date_of_Draw == ymd("2012-10-29") & !is.na(lympht_only$`Lymphocytes_%`)),]$`Lymphocytes_%`
lympht_only<-lympht_only[-which(lympht_only$MaskID == 802857 & lympht_only$Date_of_Draw == ymd("2012-10-29")& is.na(lympht_only$Lymphocytes_ABS)),]
lympht_only[which(lympht_only$MaskID == 369952 & lympht_only$Date_of_Draw == ymd("2012-12-28") & is.na(lympht_only$`Lymphocytes_%`)),]$`Lymphocytes_%`<-lympht_only[which(lympht_only$MaskID == 369952 & lympht_only$Date_of_Draw == ymd("2012-12-28") & !is.na(lympht_only$`Lymphocytes_%`)),]$`Lymphocytes_%`
lympht_only<-lympht_only[-which(lympht_only$MaskID == 369952 & lympht_only$Date_of_Draw == ymd("2012-12-28")& is.na(lympht_only$Lymphocytes_ABS)),]


lympht_only2<-lympht_only%>%select(-Visit.x,-(`Lymphocytes_%`:Visit_new))%>%rename(Visit = "Visit.y")
test<-left_join(lympht_only2,ogtt1, by =c("MaskID","baseline", "Visit")) 
sum(!is.na(test$Date_of_Draw.y) & abs(test$Visit2.x - test$Visit2.y) <=1 )
# 14 samples can be joined by these 3 columns

# manually merge these 14 samples based only on MaskID, baseline, and Visit
newcols<-c("Visit.y","White Blood Cell Count","Lymphocytes_%","Lymphocytes_ABS" ,"%CD8+TIGIT+KLRG1+ of Total CD3")
for (i in unique(lympht_only$MaskID)){
  for (j in unique(lympht_only[lympht_only$MaskID == i,]$Date_of_Draw)){
    index_lympht = which(lympht_only$MaskID == i & lympht_only$Date_of_Draw == j)
    index_ogtt<-which(ogtt.all2$MaskID==i & ogtt.all2$Visit.x == lympht_only[index_lympht,]$Visit.y)
    if (length(index_ogtt)==1){
      if (is.na(ogtt.all2[index_ogtt,"Visit.y"]) && (abs(ogtt.all2[index_ogtt,]$Visit2- lympht_only[index_lympht,]$Visit2)<=1 ) ){
        ogtt.all2[index_ogtt,newcols]<-lympht_only[index_lympht,newcols]
        ogtt.all2<-ogtt.all2[-which(ogtt.all2$MaskID==i & ogtt.all2$Date_of_Draw == j & is.na(ogtt.all2$Visit.x)),]
      }
    }else if (length(index_ogtt) > 1){
      for (index_ogtt_i in index_ogtt){
        if (is.na(ogtt.all2[index_ogtt_i,"Visit.y"]) && (abs(ogtt.all2[index_ogtt_i,]$Visit2- lympht_only[index_lympht,]$Visit2)<=1 ) ){
          ogtt.all2[index_ogtt_i,newcols]<-lympht_only[index_lympht,newcols]
          ogtt.all2<-ogtt.all2[-which(ogtt.all2$MaskID==i & ogtt.all2$Date_of_Draw == j & is.na(ogtt.all2$Visit.x)),]
        }
      }
    }
  }
}
# only 13 samples got matched
# 656446 on 2015-08-21 did not get matched b/c its closes match for ogtt is on 2015-08-03, already have lympht data


# now match more samples based on different columns
lympht_only3<-ogtt.all2[is.na(ogtt.all2$Visit_new) ,]# these timepoints dont ogtt for same Date_of_Draw
lympht_only3%>%count(Date_of_Draw,MaskID) %>% filter(n > 1) 
# again merge info for 369952 to one line, so last command returns <0 rows>
lympht_only3[which(lympht_only3$MaskID == 369952 & lympht_only3$Date_of_Draw == ymd("2012-12-28") & is.na(lympht_only3$`Lymphocytes_%`)),]$`Lymphocytes_%`<-lympht_only3[which(lympht_only3$MaskID == 369952 & lympht_only3$Date_of_Draw == ymd("2012-12-28") & !is.na(lympht_only3$`Lymphocytes_%`)),]$`Lymphocytes_%`
lympht_only3<-lympht_only3[-which(lympht_only3$MaskID == 369952 & lympht_only3$Date_of_Draw == ymd("2012-12-28")& is.na(lympht_only3$Lymphocytes_ABS)),]

lympht_only4<-lympht_only3%>%select(-Visit.x,-(`Lymphocytes_%`:Visit_new))%>%rename(Visit = "Visit.y")
test2<-left_join(lympht_only4,ogtt1, by =c("MaskID","baseline","Visit2") ) 
dim(ogtt.all2[(ogtt.all2$Date_of_Draw%in%test2$Date_of_Draw.y & ogtt.all2$mu_cpep%in%test2$mu_cpep & is.na(ogtt.all2$Visit.y)), ])
# 8 more samples can be joined by these 3 columns
for (i in unique(lympht_only3$MaskID)){
  for (j in unique(lympht_only3[lympht_only3$MaskID == i,]$Date_of_Draw)){
    index_lympht = which(lympht_only3$MaskID == i & lympht_only3$Date_of_Draw == j)
    index_ogtt<-which(ogtt.all2$MaskID==i & ogtt.all2$Visit2 == lympht_only3[index_lympht,]$Visit2)
    if (length(index_ogtt)==1){
      if (is.na(ogtt.all2[index_ogtt,"Visit.y"])){
        ogtt.all2[index_ogtt,newcols]<-lympht_only3[index_lympht,newcols]
        ogtt.all2<-ogtt.all2[-which(ogtt.all2$MaskID==i & ogtt.all2$Date_of_Draw == j & is.na(ogtt.all2$Visit.x)),]
      }
    }else if (length(index_ogtt) > 1){
      for (index_ogtt_i in index_ogtt){
        if (is.na(ogtt.all2[index_ogtt_i,"Visit.y"])){
          ogtt.all2[index_ogtt_i,newcols]<-lympht_only3[index_lympht,newcols]
          ogtt.all2<-ogtt.all2[-which(ogtt.all2$MaskID==i & ogtt.all2$Date_of_Draw == j & is.na(ogtt.all2$Visit.x)),]
        }
      }
    }
  }
}

ogtt.all2%>% count(MaskID,Date_of_Draw,mu_cpep) %>% filter(n > 1) # of those matched with ogtt, 2 sample (369952,425198) 2 lines of different lymph data for same draw date


ogtt.all2%<>%select(-Visit.y)%>%rename(Visit="Visit.x")

ogtt.all2<-full_join(tx.outcome,ogtt.all2)

write.table(ogtt.all2,"TN10full_ogtt_table.txt", sep="\t", quote = F, row.names = F)
