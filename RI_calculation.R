require(readxl)
require(stringr)
require(comprehenr)
require(dplyr)
require(tidyverse)
require(DescTools)
require(ggplot2)

excelFormat<-function(plateFile1,plateFile2, skipLine=1){
  # read in 2 excel files corresponding to plate 1-7, 8-15
  # return a table w/ samples as rows and RIs as columns
  
  
  # format column names
  row1_plate1<- as.character(read_excel(plateFile1, n_max = 1, col_names = FALSE))
  row1_plate1_split<-str_split(row1_plate1, " ")
  colnames_plate1<-to_vec(for (i in row1_plate1_split) if(length(i)>1) i[2])
  colnames_plate1<-append(colnames_plate1,"sampleID",0)
  
  # read plate files
  plate1<-read_excel(plateFile1, skip = skipLine, col_names = colnames_plate1)
  plate2<-read_excel(plateFile2, skip = skipLine, col_names = colnames_plate1)
  
  
  # format tables
  plate1<-apply(plate1,2,function(x)gsub(',', '.',x))
  plate1<-apply(plate1,2,function(x)str_extract(x, "[:graph:]*"))
  
  plate2<-apply(plate2,2,function(x)gsub(',', '.',x))
  plate2<-apply(plate2,2,function(x)str_extract(x, "[:graph:]*"))
  # if want to use parse_number() from tidyverse:
  # but this method doesnt work for scienticic notation, e.g. 7,33E-3 %
  # plate1<-as.data.frame(plate1)
  # plate2<-as.data.frame(plate2)
  # plate1[,2:7]<-apply(plate1[,2:7],2,function(x)parse_number(x, locale = locale(decimal_mark = ",")))
  
  # combine plates
  plateAll<-rbind(plate1,plate2)
  plateAll<-data.frame(plateAll)
  (head(plateAll))

  # format sampleIDs
  col1_plate_split<-lapply(plateAll[,1],function(x)str_split(x, "_"))
  col1_plate<-to_vec(for (i in col1_plate_split) i[[1]][1])
  plateAll$sampleID<-col1_plate
  
  return(plateAll)
}
  

nonCenteredAUC<-function(plateAll){
  
  plateAll$conc<-rep(c("Red","Blue","Yellow","Green"),dim(plateAll)[1]/4)
  (dilution = log(c(1/1350,1/450,1/150,1/50)))
  colnames_plateAll<-colnames(plateAll)
  
  table1<-plateAll%>%select(sampleID,conc,colnames_plateAll[2])%>%pivot_wider(names_from=conc, values_from=colnames_plateAll[2])
  table1[,2:5]<-apply(table1[,2:5],2,as.numeric)
  
  table1<-table1%>%
    mutate(ri1=apply(table1[,2:5],1,function(y)AUC(x=dilution,y=y, method = "spline")))
  
  riOut<-table1%>%select(sampleID,ri1)
  colnames(riOut)[2]<-paste(colnames_plateAll[2],"RI",sep=".")
  
  for (i in seq(3,dim(plateAll)[2]-1)){
    table2<-plateAll%>%select(sampleID,conc,colnames_plateAll[i])%>%pivot_wider(names_from=conc, values_from=colnames_plateAll[i])
    table2[,2:5]<-apply(table2[,2:5],2,as.numeric)
    table2<-table2%>%mutate(ri2=apply(table2[,2:5],1,function(y)AUC(x=dilution,y=y, method = "spline")))
    riOut<-merge(select(table2,sampleID,ri2),riOut, by = "sampleID" )
    colnames(riOut)[which(colnames(riOut) == "ri2")]<-paste(colnames_plateAll[i],"RI",sep=".")
  }
  return(riOut)
  
}

centerAUCwQC<-function(plateAll, useMedian=F){
  
  # add a column for concentration
  plateAll$conc<-rep(c("Red","Blue","Yellow","Green"),dim(plateAll)[1]/4)
  (dilution = log(c(1/1350,1/450,1/150,1/50)))
  
  # take out each column of isotype*bacteria separately and put to an RI table
  colnames_plateAll<-colnames(plateAll)
  # starts from the 2nd column b/c first column is sampleID
  table1<-plateAll%>%select(sampleID,conc,colnames_plateAll[2])%>%pivot_wider(names_from=conc, values_from=colnames_plateAll[2])
  table1[,2:5]<-apply(table1[,2:5],2,as.numeric)
  table1<-as.data.frame(table1)
  table1$qc<-apply(table1[,2:5],1,function(x)all(x[x>=5] == cummax(x[x>=5])))
  if (useMedian){
    mean_dil_curve_table1<-apply(table1[,2:5],2,median)
  }else{
    mean_dil_curve_table1<-colMeans(table1[,2:5])
  }
  mean_AUC_table1<-AUC(x=dilution, y=mean_dil_curve_table1, method="spline")
  table1<-table1%>%
    mutate(ri1=apply(table1[,2:5],1,function(y)AUC(x=dilution,y=y, method = "spline")-mean_AUC_table1))
  
  
  riOut<-table1%>%select(sampleID,ri1)
  colnames(riOut)[2]<-paste(colnames_plateAll[2],"RI",sep=".")
  
  qcOut<-table1%>%select(sampleID,qc)
  colnames(qcOut)[2]<-paste(colnames_plateAll[2],"qc",sep=".")
  
  for (i in seq(3,dim(plateAll)[2]-1)){
    table2<-plateAll%>%select(sampleID,conc,colnames_plateAll[i])%>%pivot_wider(names_from=conc, values_from=colnames_plateAll[i])
    table2[,2:5]<-apply(table2[,2:5],2,as.numeric)
    table2<-as.data.frame(table2)
    table2$qc<-apply(table2[,2:5],1,function(x)all(x[x>=5] == cummax(x[x>=5])))
    if (useMedian){
      mean_dil_curve_table2<-apply(table2[,2:5],2,median)
    }else{
      mean_dil_curve_table2<-colMeans(table2[,2:5])
    }
    mean_AUC_table2<-AUC(x=dilution, y=mean_dil_curve_table2, method="spline")
    table2<-table2%>%mutate(ri2=apply(table2[,2:5],1,function(y)AUC(x=dilution,y=y, method = "spline")-mean_AUC_table2))
    riOut<-merge(select(table2,sampleID,ri2),riOut, by = "sampleID" )
    colnames(riOut)[which(colnames(riOut) == "ri2")]<-paste(colnames_plateAll[i],"RI",sep=".")
    qcOut<-merge(select(table2,sampleID,qc),qcOut, by = "sampleID" )
    colnames(qcOut)[which(colnames(qcOut) == "qc")]<-paste(colnames_plateAll[i],"qc",sep=".")   
  }
  qcOut$summ<-apply(qcOut[,2:dim(qcOut)[2]],1,function(x)any(x==TRUE))
  return(list(riOut,qcOut))
  
}

lm_eqn <- function(df, x_ind,y_ind){
  y <- df[,y_ind]
  x <- df[,x_ind]
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

smallLargeComp<-function(dataframe){
  for (i in seq(2,dim(dataframe)[2]%/%2+1)){
    ggplot(dataframe,aes(x=dataframe[,i],y=dataframe[,i+dim(dataframe)[2]%/%2]))+geom_point()+
      geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
      geom_text(x=mean(dataframe[,i])+10,y=mean(dataframe[,i+dim(dataframe)[2]%/%2])-10,label = lm_eqn(dataframe,i,i+dim(dataframe)[2]%/%2), parse = TRUE)+
      xlab(paste(colnames(dataframe)[i],"small",sep="_"))+ylab(colnames(dataframe)[i+dim(dataframe)[2]%/%2])
    ggsave(paste(colnames(dataframe)[i], "png", sep = "."))
  }
}

smallLargeCompBA<-function(dataframe){
  for (i in seq(2,dim(dataframe)[2]%/%2+1)){
    toplot<-data.frame(sampleID=dataframe[,1],small=dataframe[,i],large=dataframe[,i+dim(dataframe)[2]%/%2])%>%
      mutate(meanPos = (small+large)/2,
             diffPos =  small-large,
             avg_diff=mean(diffPos,na.rm = T),
             sd_diff=sd(diffPos,na.rm = T),
             lower = avg_diff-1.96*sd_diff,
             upper = avg_diff+1.96*sd_diff)
    ggplot(toplot, aes(x=meanPos, y=diffPos))+geom_point(size = 3)+
      geom_hline(yintercept = toplot$avg_diff[1], color = "grey")+
      geom_hline(yintercept = toplot$lower[1], color = "grey", linetype="dashed")+
      geom_hline(yintercept = toplot$upper[1], color = "grey", linetype="dashed")+
      # ggtitle(paste("Bland-Altman Plot for",unique(.$bug.iso))) +
      ylab("Difference in positive events") +
      xlab(paste(colnames(dataframe)[i], "Mean positive events", sep = "."))+
      geom_text(data=subset(toplot,(diffPos > upper | diffPos < lower)),aes(label=sampleID),vjust = -1.5)
    ggsave(paste(colnames(dataframe)[i], "BA.png", sep = "_"))
  }
}

