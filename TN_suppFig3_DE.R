library(tidyverse)
library(magrittr)
library(wesanderson)
library(factoextra)
library(FactoMineR)
library(car) # for ANOVA assumption test
tn10<-read.csv("TN10full_smallgate.csv", header = T, as.is = T)
tn10$Date_of_Draw<-as.Date(tn10$Date_of_Draw)
tn10$Date.IDDM<-as.Date(tn10$Date.IDDM)
tn10$T1D.event<-as.factor(tn10$T1D.event)
tn10$DR3<-factor(tn10$DR3, levels = c("ABSENT","PRESENT"))
tn10$DR4<-factor(tn10$DR4, levels = c("ABSENT","PRESENT"))
tn10$Treatment.Arm<-factor(tn10$Treatment.Arm, levels = c("Placebo","Teplizumab"))
str(tn10)
theme_set(
  theme_classic(base_size = 20)
)

# PCA, baseline -------------------------------------------------------

acab.allIg.base<-tn10%>%filter(Visit == 0)%>%
  filter(!is.na(Bd.IgG2.RI))%>%distinct(MaskID, .keep_all = T)%>%
  mutate(Tfreq = ((X.CD8.TIGIT.KLRG1..of.Total.CD3+1)**(-0.519789386)-1)/-0.519789386)

acab.allIg.base.forPCA<-acab.allIg.base
rownames(acab.allIg.base.forPCA)<-acab.allIg.base.forPCA$MaskID
rownames(acab.allIg.base.forPCA)
acab.allIg.base.forPCA%<>%dplyr::select(contains("Ig"),-X.CD8.TIGIT.KLRG1..of.Total.CD3)
acab.allIg.base.forPCA.res<- prcomp(acab.allIg.base.forPCA, scale = T)
fviz_eig(acab.allIg.base.forPCA.res,addlabels = T)
varexp_toplotBase<-get_eigenvalue(acab.allIg.base.forPCA.res)

# Contributions of variables to PC2
fviz_contrib(acab.allIg.base.forPCA.res, choice = "var", axes = 2, top = 10)
fviz_cos2(acab.allIg.base.forPCA.res, choice = "var", axes = 2)

# Color variable contributions by isotypes
# Split a character vector by a specific character; save the 2nd element in new vector
varlist<-factor(unlist(lapply(strsplit(colnames(acab.allIg.base.forPCA), split='.',fixed = T), `[`, 2)))
fviz_pca_var(acab.allIg.base.forPCA.res,
             
             habillage = varlist,
             # alpha.var = "contrib",
             palette = c("#00AFBB","#FC4E07", "#E7B800"),
             repel = TRUE     # Avoid text overlapping
)+labs(color = "Isotype",title ="")+theme_classic(base_size = 20)

# PCA, 6 months -------------------------------------------------------

acab.allIg.6month<-tn10%>%filter(Visit == 6 & !is.na(Bd.IgG2.RI))
acab.allIg.6month.forPCA<-acab.allIg.6month
rownames(acab.allIg.6month.forPCA)<-acab.allIg.6month.forPCA$MaskID
rownames(acab.allIg.6month.forPCA)
acab.allIg.6month.forPCA%<>%dplyr::select(contains("Ig"),-X.CD8.TIGIT.KLRG1..of.Total.CD3)
acab.allIg.6month.forPCA.res<- prcomp(acab.allIg.6month.forPCA, scale = T)

fviz_eig(acab.allIg.6month.forPCA.res,addlabels = TRUE)
fviz_contrib(acab.allIg.6month.forPCA.res, choice = "var", axes = 2)

# variable contribution plot
varlist<-factor(unlist(lapply(strsplit(colnames(acab.allIg.6month.forPCA), split='.',fixed = T), `[`, 2)))
fviz_pca_var(acab.allIg.6month.forPCA.res,
             habillage = varlist,
             # alpha.var = "contrib",
             palette = c("#00AFBB", "#FC4E07", "#E7B800"),
             repel = TRUE     # Avoid text overlapping
)+labs(color = "Isotype",title="")+theme_classic(base_size = 20)