#Step 1 statistical analysis
rm(list=ls())
library(stringr)
library(dynamicTreeCut)
library(fastcluster)
library(WGCNA)
library(boot)
library(tidyverse)
library(ggpubr)
library(Hmisc)
library(robustbase)
library(multcomp)
library(survival)
library(survminer)

savepath<-"D:/ZZZ/Manuscripts/Subtyping/Sustain_results/MS_NMO_MOG_sustain"
setwd(savepath)
Data<-read.csv('Sustain_data_Jinyuan_MS_NMO_baseline1.csv',header = TRUE)

Data1<-read.csv('Sustain_data_Jinyuan_MS_NMO_followup.csv',header = TRUE)

rownames(Data)<-Data$PatientID_Image

Data$Sex_F1_M0<-as.factor(Data$Sex_F1_M0);
Data1$Sex_F1_M0<-as.factor(Data1$Sex_F1_M0);

for (i in 15:33)
{fit<-lm(Data[,i]~Age+Sex_F1_M0+EstimatedTotalIntraCranialVol,data=Data)
   fit_value<-predict(fit,newdata=Data)
   Data[,i]<-mean(fit_value)+fit$residuals;
   
   fit_value1<-predict(fit,newdata=Data1)
   Data1[,i+2]<-mean(fit_value1,rm.na=TRUE)+Data1[,i+2]-fit_value1;
   
}

write.csv(Data1,'Sustain_adjusted_data_follow_up.csv')