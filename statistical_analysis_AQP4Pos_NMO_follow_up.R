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
library(pheatmap)

savepath<-"D:/ZZZ/Manuscripts/Subtyping/Sustain_results"
setwd(savepath)
Data0<-read.csv('D:/ZZZ/Manuscripts/Subtyping/Sustain_results/MS_AQP4Pos_NMO/WorkshopOutput_1212_nmo/newAQP4Pos_NMOSDtest.csv',header = TRUE)
rownames(Data0)<-Data0$PatientID_Image


#read MRI metrics 
data_all1<-read.csv('D:/ZZZ/Manuscripts/Subtyping/Sustain_results/New_Clinical_Validation_information1.csv',header = TRUE);

data_all1<-data_all1[!is.na(data_all1$Group)&!is.na(data_all1$Sex_F1_M0)&!is.na(data_all1$Age),]

data_all2<-read.csv('D:/ZZZ/Manuscripts/Subtyping/Sustain_results/tem_clinical_data_Volume_20221119_complete.csv',header = TRUE);
#select data
data_all2<-data_all2[data_all2$baseline=='baseline'&is.na(data_all2$Repeat_data)&!is.na(data_all2$Group)&
                       !is.na(data_all2$Sex_F1_M0)&!is.na(data_all2$Age),]

Nor_var<-c('PatientID_Image','Group','Age','Sex_F1_M0','Education_years','Duration_Month','relapses',
           'WMH_volume','WM.hypointensities.1',
           'CortexVol','SubCortGrayVol','CerebralWhiteMatterVol',
           'SUS_Brainstem','SUS_CerebellumGM','SUS_CerebellumWM','SUS_MUCCA123','EstimatedTotalIntraCranialVol',
           'choroid_volume','MMSE','MOCA','CVLT','BVLT','PASAT','EDSS')
m2<-dim(data_all2);m1<-dim(data_all1)


data_all3<-rbind(data_all2[,Nor_var],data_all1[,Nor_var]);

add_var<-c('SDMT','COWAT','Beton',
           'Followup_relapse','Followup_EDSS','Followup_time_years',
           'EDSS_Progress','SPMS_conversion',
           
           'T1T2Ratio_cerebral_GM','T1T2Ratio_cerebral_WM','T1T2Ratio_cerebellum_GM','T1T2Ratio_cerebellum_WM',
           'T1T2Ratio_brainstem','FA_cerebral_GM','FA_cerebral_WM','FA_cerebellum_GM','FA_cerebellum_WM','FA_brainstem',
           'MD_cerebral_GM','MD_cerebral_WM','MD_cerebellum_GM','MD_cerebellum_WM','MD_brainstem',
           'DC_Cerebral_GM','DC_Cerebellum_GM',
           'fALFF_Cerebral_GM','fALFF_Cerebellum_GM',
           'NFL','GFAP','AB40','AB42','TAU','TREM2',
           'SUS_MUCCA123','CortexVol','SubCortGrayVol','CerebralWhiteMatterVol',
           'SUS_CerebellumGM','SUS_CerebellumWM','SUS_Brainstem','EstimatedTotalIntraCranialVol')

data_all3[1:m2[1],add_var]<-data_all2[,add_var]

  
Data_clinical<-data_all3
  
#Data_clinical<-read.csv('tem_clinical_data_Volume_20221119_complete.csv',header = TRUE)

#Data_clinical<-Data_clinical[is.na(Data_clinical$Repeat_data),]

#Data_clinical<-Data_clinical[!is.na(Data_clinical$Group)&!is.na(Data_clinical$Sex_F1_M0)&!is.na(Data_clinical$Age),]

#Data_clinical<-Data_clinical[Data_clinical$Group=='MS'|Data_clinical$Group=='AQP4Pos'|Data_clinical$Group=='AQP4Neg'|Data_clinical$Group=='MOGAD',]
rownames(Data_clinical)<-Data_clinical$PatientID_Image

index<-intersect(rownames(Data0),rownames(Data_clinical))

Data0<-Data0[index,]
Data_clinical<-Data_clinical[index,]

Data<-Data_clinical;

Data[,'ml_subtype']<-Data0$ml_subtype
Data[,'ml_stage']<-Data0$ml_stage

Data[,'prob_ml_subtype']<-Data0$prob_ml_subtype
Data[,'prob_ml_stage']<-Data0$prob_ml_stage


setwd(paste0(savepath,"/Results/AQP4_NMOSD"))

Data<-Data[!(is.na(Data$ml_subtype)|Data$ml_subtype==""),]
#Data<-Data[Data$ml_stage>0,] #remove stage 0

Data[Data$Group=='AQP4Pos_NMOSD'&Data$ml_stage==0,"Group_subtype"]='B0';
Data[Data$Group=='AQP4Pos_NMOSD'&Data$ml_subtype==0&Data$ml_stage>0,"Group_subtype"]='B1';
Data[Data$Group=='AQP4Pos_NMOSD'&Data$ml_subtype==1&Data$ml_stage>0,"Group_subtype"]='B2';
Data[Data$Group=='AQP4Pos_NMOSD'&Data$ml_subtype==2&Data$ml_stage>0,"Group_subtype"]='B3';

Data<-Data[Data$Group=='AQP4Pos_NMOSD',]

Data$Group_subtype<-factor(Data$Group_subtype,levels=c('B0','B1','B2','B3'))


Data_fl<-read.csv('D:/ZZZ/Manuscripts/Subtyping/Sustain_results/Results/follow_up/newAQP4Pos_NMOSDeval.csv',header = TRUE)


rownames(Data_fl)<-Data_fl$PatientID_Image

# index1<-intersect(rownames(Data),rownames(Data_fl))
# Data<-Data[index,]
# Data_fl<-Data_fl[index,]

#Data_new<-

Data_fl$baseline_for_follow_up

Data_fl$follow_up_time_days

index2<-NULL;
for (j in rownames(Data_fl))
  
{
  if (Data_fl[j,'follow_up_time_days']>60)
  {
    index2<-c(index2,j)
  }
  
  if (Data_fl[j,'follow_up_time_days']<360)
  {
    Data_fl[j,'fl_years']<-1;
  }
    
    else if (Data_fl[j,'follow_up_time_days']<720)
    {
    Data_fl[j,'fl_years']<-2; 
    }
  
    else if (Data_fl[j,'follow_up_time_days']<1080)
    {
    Data_fl[j,'fl_years']<-3; 
    }
    

}

Data_fl<-Data_fl[index2,]


n=0;
for(i in rownames(Data))
{
  n1=0;
  for (j in rownames(Data_fl))

    {
     if (i==Data_fl[j,'baseline_for_follow_up'])
     {     
         n1=n1+1;
         Data[i,'individual']<-i
         Data[i,'follow_up_time_days']<-0
         Data_fl[j,'individual']<-i
         Data[i,'fl_point']<-0
         Data_fl[j,'fl_point']<-n1
         Data[i,'fl_years']<-0
         
     }
  
    }

}


Data_fl[Data_fl$Group=='AQP4Pos_NMOSD'&Data_fl$ml_stage==0,"Group_subtype"]='B0';
Data_fl[Data_fl$Group=='AQP4Pos_NMOSD'&Data_fl$ml_subtype==0&Data_fl$ml_stage>0,"Group_subtype"]='B1';
Data_fl[Data_fl$Group=='AQP4Pos_NMOSD'&Data_fl$ml_subtype==1&Data_fl$ml_stage>0,"Group_subtype"]='B2';
Data_fl[Data_fl$Group=='AQP4Pos_NMOSD'&Data_fl$ml_subtype==2&Data_fl$ml_stage>0,"Group_subtype"]='B3';

col_index<-intersect(colnames(Data_fl),colnames(Data))

Data<-Data[,col_index]
Data_fl<-Data_fl[,col_index]



Data_new<-rbind(Data,Data_fl);

Data_new<-Data_new[!is.na(Data_new$individual),]

Data_new[Data_new$Group_subtype=='B0', 'Group_Index']=1;
Data_new[Data_new$Group_subtype=='B1', 'Group_Index']=2;
Data_new[Data_new$Group_subtype=='B2', 'Group_Index']=3;
Data_new[Data_new$Group_subtype=='B3', 'Group_Index']=4;


tem_name1<-Data_new[Data_new$fl_years==0,'individual']
tem_name2<-Data_new$individual


Subtype_matrix<-data.frame(matrix(NA,length(tem_name1),max(Data_new$fl_years)+1));
Subtype_pro_matrix<-data.frame(matrix(NA,length(tem_name1),max(Data_new$fl_years)+1));
Stage_matrix<-data.frame(matrix(NA,length(tem_name1),max(Data_new$fl_years)+1));
Stage_pro_matrix<-data.frame(matrix(NA,length(tem_name1),max(Data_new$fl_years)+1));


for(i in 1:length(tem_name1))
{
  for (j in 1:length(tem_name2))
    
  {
    if (tem_name1[i]==tem_name2[j])
    {     
      
      Subtype_matrix[i,Data_new$fl_years[j]+1]<-Data_new$Group_Index[j];
      
      Subtype_pro_matrix[i,Data_new$fl_years[j]+1]<-Data_new$prob_ml_subtype[j];
      
      Stage_matrix[i,Data_new$fl_years[j]+1]<-Data_new$ml_stage[j];
      
      Stage_pro_matrix[i,Data_new$fl_years[j]+1]<-Data_new$prob_ml_stage[j];
      
      
    }
    
  }
  
}


#identify the stable matrix
del_index<-NULL
for (i in 1:dim(Subtype_matrix)[1])
{
  n=0;
  for (j in 1:dim(Subtype_matrix)[2])
  {
    if (Subtype_matrix[i,j]!=Subtype_matrix[i,1]&&!is.na(Subtype_matrix[i,j]))
        {
          n=1;
        }
    
  }
  
  if (n==1)
  {
    del_index<-c(del_index,i);
  }
}


pdf(file=paste0('Subtype_matrix.pdf'),width=5,height=10)
p<-pheatmap(Subtype_matrix,display_numbers =Stage_matrix,number_color = 'black',fontsize_number = 10, 
            cluster_rows=F,cluster_cols=F,color=colorRampPalette(c('lightblue','white','firebrick3'))(4),
                 legend_breaks=c(1, 2, 3, 4),legend_labels=c(1,2,3,4),na_col = 'white')
print(p)
dev.off()

#pheatmap(Stage_matrix,display_numbers = T, cluster_rows=F,cluster_cols=F,color=colorRampPalette(c('navy','white','firebrick3'))(40))


library(ROCR)
#inlude normal appearing
conMatrix1<-table(Subtype_matrix[!is.na(Subtype_matrix[,2]),1],Subtype_matrix[!is.na(Subtype_matrix[,2]),2])

conMatrix2<-table(Subtype_matrix[!is.na(Subtype_matrix[,3]),1],Subtype_matrix[!is.na(Subtype_matrix[,3]),3])

conMatrix3<-table(Subtype_matrix[!is.na(Subtype_matrix[,4]),1],Subtype_matrix[!is.na(Subtype_matrix[,4]),4])

write.csv(conMatrix1,'follow_up_subtype_1y.csv')
write.csv(conMatrix2,'follow_up_subtype_2y.csv')
write.csv(conMatrix3,'follow_up_subtype_3y.csv')


#exlude normal appearing
conMatrix1<-table(Subtype_matrix[!is.na(Subtype_matrix[,2])&Subtype_matrix[,1]>1,1],Subtype_matrix[!is.na(Subtype_matrix[,2])&Subtype_matrix[,1]>1,2])

conMatrix2<-table(Subtype_matrix[!is.na(Subtype_matrix[,3])&Subtype_matrix[,1]>1,1],Subtype_matrix[!is.na(Subtype_matrix[,3])&Subtype_matrix[,1]>1,3])

conMatrix3<-table(Subtype_matrix[!is.na(Subtype_matrix[,4])&Subtype_matrix[,1]>1,1],Subtype_matrix[!is.na(Subtype_matrix[,4])&Subtype_matrix[,1]>1,4])


write.csv(conMatrix1,'follow_up_subtype_1y_del_NA.csv')
write.csv(conMatrix2,'follow_up_subtype_2y_del_NA.csv')
write.csv(conMatrix3,'follow_up_subtype_3y_del_NA.csv')


Subtype_matrix_stable<-Subtype_matrix[-del_index,]
Subtype_pro_matrix_stable<-Subtype_pro_matrix[-del_index,]
Stage_matrix_stable<-Stage_matrix[-del_index,]
Stage_pro_matrix_stable<-Stage_pro_matrix[-del_index,]

Subtype_matrix_unstable<-Subtype_matrix[-del_index,]
Subtype_pro_matrix_unstable<-Subtype_pro_matrix[-del_index,]
Stage_matrix_unstable<-Stage_matrix[-del_index,]
Stage_pro_matrix_unstable<-Stage_pro_matrix[-del_index,]


pdf(file=paste0('Subtype_matrix_stable.pdf'),width=5,height=10)
p<-pheatmap(Subtype_matrix_stable,display_numbers =Stage_matrix_stable,number_color = 'white',fontsize_number = 10, cluster_rows=F,cluster_cols=F,color=colorRampPalette(c('navy','white','firebrick3'))(4),
         legend_breaks=c(1, 2, 3, 4),legend_labels=c(1,2,3,4),na_col = 'white')
print(p)
dev.off()

# pdf(file=paste0('(Stage_matrix_stable.pdf'),width=10,height=20)
# pheatmap(Stage_matrix_stable,display_numbers = T,number_color = 'white',fontsize_number = 10, cluster_rows=F,cluster_cols=F,color=colorRampPalette(c('navy','firebrick3'))(40),
#          legend_breaks=c(0, 10, 20, 30,40),legend_labels=c(0,10,20,30,40),na_col = 'white')
# dev.off()

# # library(caret)
# # confusionMatrix(Subtype_matrix$x1,Subtype_matrix$x2)
# 
# colnames(Data)
# 
# my_comparisons<-list(c('B0','B1'),
#                    c('B0','B2'),
#                    c('B0','B3'),
#                    
#                    c('B1','B2'),
#                    c('B1','B3'),
#                    
#                    c('B2','B3'))
# 
# my_comparisons1<-list(
#                      c('B1','B2'),
#                      c('B1','B3'),
#                      
#                      c('B2','B3'))
# 
# 
# corrected_method='fdr';
# #ml_stage
# Data1<-Data
# Data1<-Data1[!(is.na(Data1$ml_stage)|Data1$ml_stage==''),]
# Group_subtype1<-factor(Data1$Group_subtype)
# Fc<-lm(Data1$ml_stage~Group_subtype1)
# ml_stage_results<-summary(glht(Fc, linfct=mcp(Group_subtype1="Tukey")))
#
# pdf(file='ml_stage.pdf')
# p1 = ggboxplot(Data1, x = "Group_subtype", y = "ml_stage",
#                color = "Group_subtype", palette = "jco",add = "jitter",legend="none")+
#   
#   stat_compare_means(comparisons = my_comparisons1,
#                      
#                      p.adjust.method = corrected_method,
#                      method='wilcox.test',
#                      hide.ns = F,
#                      #label="p.adj",
#                      # p.adjust.method = "Tukey",
#                      # hide.ns=TRUE
#                      
#                      label = "p.fromat"
#   )+
#   
#   stat_compare_means(method='kruskal.test')
# print(p1)
# dev.off()
# 
# 
# #Age
# Data1<-Data
# Data1<-Data1[!(is.na(Data1$Age)|Data1$Age==''),]
# Group_subtype1<-factor(Data1$Group_subtype)
# Fc<-lm(Data1$Age~Group_subtype1)
# Age_results<-summary(glht(Fc, linfct=mcp(Group_subtype1="Tukey")))
# 
# pdf(file='age.pdf')
# p1 = ggboxplot(Data1, x = "Group_subtype", y = "Age",
#                color = "Group_subtype", palette = "jco",add = "jitter",legend="none")+
#   
#   stat_compare_means(comparisons = my_comparisons,
#                      
#                      p.adjust.method = corrected_method,
#                      method='wilcox.test',
#                      hide.ns = F,
#                      #label="p.adj",
#                      # p.adjust.method = "Tukey",
#                      # hide.ns=TRUE
#                      
#                      label = "p.fromat"
#   )+
#   
#   stat_compare_means(method='kruskal.test')
# print(p1)
# dev.off()
# 
# #sex
# Data1<-Data
# Data1<-Data1[!(is.na(Data1$Sex_F1_M0)|Data1$Sex_F1_M0=="[]"|Data1$Sex_F1_M0==""|Data1$Sex_F1_M0=="NaN"),]
# Group_subtype1<-factor(Data1$Group_subtype)
# Data1$Sex_F1_M0<-as.numeric(Data1$Sex_F1_M0)
# 
# tab12<-table(Data1$Group_subtype,Data1$Sex_F1_M0)
# 
# 
# # df <- as.data.frame(tab12)
# # 
# # colnames(df)<-c('Subtype','Sex','Freq')
# # 
# # ggbarplot(df, x = "Subtype", y = "Freq",
# #           color = "Sex", position = position_dodge(),
# #           #palette = c("brown", "blue", "black", "red")
# #           )
# pdf(file='Sex.pdf')
# p1<-ggplot(Data1, aes(x=as.character(Group_subtype),fill=as.factor(Sex_F1_M0)))+
#   geom_bar(stat="count",width=0.5,position='dodge')+
#   #scale_fill_manual(values=c('#304156','#e69f00'))+
#   geom_text(stat='count',aes(label=..count..), color="black", 
#             size=5,position=position_dodge(0.5),vjust=-0.5 )+
#   ggtitle("堆砌柱状图") +
#   xlab("Subtype")+
#   labs(fill="Sex")+
#   theme_minimal()
# 
# print(p1)
# dev.off()
# 
# #Education_year
# Data1<-Data
# Data1<-Data1[!(is.na(Data1$Education_years)|Data1$Education_year==""|Data1$Education_years=="NaN"),]
# Group_subtype1<-factor(Data1$Group_subtype)
# Data1$Education_years<-as.numeric(Data1$Education_years)
# Fc<-lm(Data1$Education_years ~Group_subtype1)
# Education_years_results<-summary(glht(Fc, linfct=mcp(Group_subtype1="Tukey")))
# 
# pdf(file='Education_years.pdf')
# p1 = ggboxplot(Data1, x = "Group_subtype", y = "Education_years",
#                color = "Group_subtype", palette = "jco",add = "jitter")+
#   
#   stat_compare_means(comparisons = my_comparisons,
#                      
#                      p.adjust.method =corrected_method,
#                      method='wilcox.test',
#                      hide.ns = F,
#                      #label="p.adj",
#                      # p.adjust.method = "Tukey",
#                      # hide.ns=TRUE
#                      
#                      label = "p.fromat"
#   )+
#   stat_compare_means(method='kruskal.test')
# print(p1)
# dev.off()
# 
# 
# 
# #Disease duration
# Data1<-Data
# Data1<-Data1[!(is.na(Data1$Duration_Month)|Data1$Duration_Month==""|Data1$Duration_Month=="NaN"),]
# Group_subtype1<-factor(Data1$Group_subtype)
# Data1$Duration_Month<-as.numeric(Data1$Duration_Month)
# Fc<-lm(Data1$Duration_Month ~Group_subtype1)
# Education_years_results<-summary(glht(Fc, linfct=mcp(Group_subtype1="Tukey")))
# 
# pdf(file='Duration_Month.pdf')
# p1 = ggboxplot(Data1, x = "Group_subtype", y = "Duration_Month",
#                color = "Group_subtype", palette = "jco",add = "jitter")+
#   
#   stat_compare_means(comparisons = my_comparisons,
#                      
#                      p.adjust.method = corrected_method,
#                      method='wilcox.test',
#                      hide.ns = F,
#                      #label="p.adj",
#                      # p.adjust.method = "Tukey",
#                      # hide.ns=TRUE
#                      
#                      label = "p.fromat"
#   )+
#   stat_compare_means(method='kruskal.test')
# print(p1)
# dev.off()
# 
# 
# 
# var<-c('relapses','EDSS','Followup_relapse','Followup_EDSS')
# 
# library(lme4)
# 
# for (i in var)
# {Data1<-Data
# Data1<-Data1[!(is.na(Data1[,i])|Data1[,i]==""|Data1[,i]=="NaN"),]
# Data1$Group_subtype<-factor(Data1$Group_subtype)
# 
# Data1$relapses<-as.numeric(Data1$relapses)
# 
# Fc1<-lm(Data1[,i] ~Data1$Group_subtype+Data1$ml_stage+Data1$Age+Data1$Sex_F1_M0)
# print(i)
# print(summary(Fc1))
# 
# 
# 
# pdf(file=paste0(i,'.pdf'))
# p1 = ggboxplot(Data1, x = "Group_subtype", y = i,
#                color = "Group_subtype", palette = "jco",add = "jitter")+
#   
#   stat_compare_means(comparisons = my_comparisons,
#                      
#                      p.adjust.method = corrected_method,
#                      method='wilcox.test',
#                      hide.ns = F,
#                      #label="p.adj",
#                      # p.adjust.method = "Tukey",
#                      # hide.ns=TRUE
#                      
#                      label = "p.fromat"
#   )+
#   stat_compare_means(method='kruskal.test')
# print(p1)
# dev.off()
# 
# 
# 
# }
# 
# 
# 
# var<-c('MOCA','MMSE','SDMT','CVLT','BVLT','COWAT','Beton','PASAT')
# 
# for (i in var)
# {Data1<-Data
# Data1<-Data1[!(is.na(Data1[,i])|Data1[,i]==""|Data1[,i]=="NaN"),]
# Data1$Group_subtype<-factor(Data1$Group_subtype)
# 
# Data1$relapses<-as.numeric(Data1$relapses)
# 
# Fc1<-lm(Data1[,i] ~Data1$Group_subtype+Data1$ml_stage+Data1$Age+Data1$Sex_F1_M0)
# print(i)
# print(summary(Fc1))
# 
# pdf(file=paste0(i,'.pdf'))
# p1 = ggboxplot(Data1, x = "Group_subtype", y = i,
#                color = "Group_subtype", palette = "jco",add = "jitter")+
#   
#   stat_compare_means(comparisons = my_comparisons,
#                      
#                      p.adjust.method = corrected_method,
#                      method='wilcox.test',
#                      hide.ns = F,
#                      #label="p.adj",
#                      # p.adjust.method = "Tukey",
#                      # hide.ns=TRUE
#                      
#                      label = "p.fromat"
#   )+
#   stat_compare_means(method='kruskal.test')
# print(p1)
# dev.off()
# }
# 
# 
# # var<-c('T1T2Ratio_cerebral_GM','T1T2Ratio_cerebral_WM','T1T2Ratio_cerebellum_GM','T1T2Ratio_cerebellum_WM',
# #       'T1T2Ratio_brainstem','FA_cerebral_GM','FA_cerebral_WM','FA_cerebellum_GM','FA_cerebellum_WM','FA_brainstem',
# #       'MD_cerebral_GM','MD_cerebral_WM','MD_cerebellum_GM','MD_cerebellum_WM','MD_brainstem',
# #       'DC_Cerebral_GM','DC_Cerebellum_GM',
# #       'fALFF_Cerebral_GM','fALFF_Cerebellum_GM',
# #       'NFL','GFAP','AB40','AB42','TAU','TREM2',
# #       'SUS_MUCCA123','CortexVol','SubCortGrayVol','CerebralWhiteMatterVol',
# #       'SUS_CerebellumGM','SUS_CerebellumWM','SUS_Brainstem','EstimatedTotalIntraCranialVol')
# 
# var<-c('WMH_volume','WM.hypointensities.1',
#        'CortexVol','SubCortGrayVol','CerebralWhiteMatterVol',
#        'SUS_Brainstem','SUS_CerebellumGM','SUS_CerebellumWM','SUS_MUCCA123','EstimatedTotalIntraCranialVol',
#        'choroid_volume',
#        'FA_cerebral_WM','FA_cerebellum_WM','FA_brainstem',
#        'MD_cerebral_WM','MD_cerebellum_WM','MD_brainstem',
#        'DC_Cerebral_GM','DC_Cerebellum_GM',
#        'fALFF_Cerebral_GM','fALFF_Cerebellum_GM')
# 
# #Data[,var]<-as.numeric(Data[,var])
# 
# for (i in var)
# {Data1<-Data
# Data1<-Data1[!(is.na(Data1[,i])|Data1[,i]==""|Data1[,i]=="NaN"|is.na(Data1[,i])|Data1[,i]=='[]'),]
# Data1$Group_subtype<-factor(Data1$Group_subtype)
# 
# 
# Fc1<-lm(Data1[,i] ~Data1$Group_subtype+Data1$ml_stage+Data1$Age+Data1$Sex_F1_M0)
# print(i)
# print(summary(Fc1))
# 
# pdf(file=paste0(i,'.pdf'))
# p1 = ggboxplot(Data1, x = "Group_subtype", y = i,
#                color = "Group_subtype", palette = "jco",add = "jitter")+
# 
#   stat_compare_means(comparisons = my_comparisons,
# 
#                      p.adjust.method =corrected_method,
#                      method='wilcox.test',
#                      hide.ns = F,
#                      #label="p.adj",
#                      # p.adjust.method = "Tukey",
#                      # hide.ns=TRUE
# 
#                      label = "p.fromat"
#   )+
#   stat_compare_means(method='kruskal.test')
# print(p1)
# dev.off()
# }
# 
# 
# var<-c('NFL','GFAP','AB40','AB42','TAU','TREM2')
# 
# for (i in var)
# {Data1<-Data
# Data1<-Data1[!(is.na(Data1[,i])|Data1[,i]==""|Data1[,i]=="NaN"|is.na(Data1[,i])|Data1[,i]=='[]'),]
# Data1$Group_subtype<-factor(Data1$Group_subtype)
# 
# 
# Fc1<-lm(Data1[,i] ~Data1$Group_subtype+Data1$ml_stage+Data1$Age+Data1$Sex_F1_M0)
# print(i)
# print(summary(Fc1))
# 
# pdf(file=paste0(i,'.pdf'))
# p1 = ggboxplot(Data1, x = "Group_subtype", y = i,
#                color = "Group_subtype", palette = "jco",add = "jitter")+
#   
#   stat_compare_means(comparisons = my_comparisons,
#                      
#                      p.adjust.method = corrected_method,
#                      method='wilcox.test',
#                      hide.ns = F,
#                      #label="p.adj",
#                      # p.adjust.method = "Tukey",
#                      # hide.ns=TRUE
#                      
#                      label = "p.fromat"
#   )+
#   stat_compare_means(method='kruskal.test')
# print(p1)
# dev.off()
# }
# 
# 
# 
# library(ggpmisc)
# 
# LETTERS<-c('B0','B1','B2','B3')
# NG<-length(LETTERS)
# Ass_var<-c('Age','Education_years','Duration_Month','relapses','EDSS',
#            'WMH_volume','WM.hypointensities.1','choroid_volume',
#            'MMSE','MOCA','CVLT','BVLT','PASAT','SDMT','COWAT','Beton',
#            'NFL','GFAP','AB40','AB42','TAU','TREM2')
# 
# Res_corr<-data.frame(matrix(0,length(Ass_var),2*NG))
# 
# for (i in 1:length(Ass_var))
# {
#   
#   Data[,'tem_ass']<-Data[,Ass_var[i]]
#   Data1<-Data
#   Data1<-Data1[!is.na(Data1$tem_ass),]
#   
#   
#   pdf(file=paste0('Stage_assocaited_with_',Ass_var[i],'.pdf'),width=20,height=3)
#   p1<-ggplot(data = Data1, aes(x =ml_stage , y = tem_ass, color=Group_subtype),palette = "jco") +
#     geom_smooth(method = "lm",size=2, se=FALSE) +
#     #scale_color_manual(values=colors_2)+
#     geom_point(size=1,
#                alpha=0.3)+
#     facet_wrap(~Group_subtype,scales="free",nrow=1)+
#     
#     # stat_fit_glance(method = 'lm',
#     #                 
#     #                 
#     #                 method.args = list(formula=y~x),
#     #                 
#     #                 mapping = aes(label = sprintf('R^2~"="~%.3f~~italic(P)~"="~%.2g', stat(r.squared), stat(p.value))),
#     #                 
#     #                 parse = TRUE,label.x = 0.95,label.y = 0.95)+
#     stat_correlation(mapping=aes(label=paste(after_stat(cor.label),
#                                              after_stat(p.value.label),
#                                              after_stat(n.label),
#                                              sep='*"; "*')))+
#     scale_x_continuous(name ='Disease Stage')+
#     
#     scale_y_continuous(#expand = c(0, 0),#设定x轴和y轴的交叉点
#       
#       name =Ass_var[i],#设定y轴标题
#       
#       #breaks=seq(0,50,10),#设定Y轴的数据间隔
#       
#       #limits = c(10,50) #设定Y轴的数据上下限
#     )+
#     
#     theme_bw()+
#     theme(panel.grid = element_blank())
#   
#   print(p1)
#   dev.off()
#   
#   
#   
#   library(ppcor)
#   
#   for (j in 1:NG)
#     
#   {
#     
#     tem_data<-Data1[Data1$Group_subtype==LETTERS[j],c('ml_stage','tem_ass')]
#     
#     if (dim(tem_data)[1]>3)
#     {  
#       cor_res<-pcor(tem_data,method='kendall')
#       
#       Res_corr[i,1+2*(j-1)]<-cor_res$estimate[2,1]
#       Res_corr[i,2+2*(j-1)]<-cor_res$p.value[2,1]
#       colnames(Res_corr)[c(1+2*(j-1),2+2*(j-1))]<-paste0(LETTERS[j],c('_R','_P'))
#       rownames(Res_corr)[i]<-Ass_var[i]
#     }
#   }
# }
# 
# 
# write.csv(Res_corr,"All_var_correlation.csv")
# #Followup_time_years
# Data1<-Data
# Data1<-Data1[!(is.na(Data1$Followup_time_years)|Data1$Followup_time_years=="[]"|Data1$Followup_time_years==""|Data1$Followup_time_years=="NaN"),]
# Group_subtype1<-factor(Data1$Group_subtype)
# Data1$Followup_time_years<-as.numeric(Data1$Followup_time_years)
# Fc1<-lm(Data1$Followup_time_years ~Data1$Age+Data1$Sex_F1_M0)
# Data1$Followup_time_years<-Fc1$residuals+mean(Data1$Followup_time_years)
# 
# Fc<-lm(Data1$Followup_time_years ~Group_subtype1)
# Beton_results<-summary(glht(Fc, linfct=mcp(Group_subtype1="Tukey")))
# 
# pdf(file='Followup_time_years.pdf')
# p1 = ggboxplot(Data1, x = "Group_subtype", y = "Followup_time_years",
#                color = "Group_subtype", palette = "jco",add = "jitter")+
#   
#   stat_compare_means(comparisons = my_comparisons,
#                      
#                      p.adjust.method =corrected_method,
#                      method='t.test',
#                      hide.ns = TRUE,
#                      label = "p.signif"
#                      #label.y=c()
#   )+
#   stat_compare_means(method='kruskal.test',label.y=0)
# print(p1)
# dev.off()
# 
# 
# 
# #EDSS_Progress
# Data1<-Data
# Data1<-Data1[!(is.na(Data1$EDSS_Progress)|Data1$EDSS_Progress=="[]"|Data1$EDSS_Progress==""|Data1$EDSS_Progress=="NaN"),]
# Group_subtype1<-factor(Data1$Group_subtype)
# Data1$EDSS_Progress<-as.numeric(Data1$EDSS_Progress)
# 
# tab12<-table(Data1$Group_subtype,Data1$EDSS_Progress)
# 
# 
# # df <- as.data.frame(tab12)
# # 
# # colnames(df)<-c('Subtype','Sex','Freq')
# # 
# # ggbarplot(df, x = "Subtype", y = "Freq",
# #           color = "Sex", position = position_dodge(),
# #           #palette = c("brown", "blue", "black", "red")
# #           )
# pdf(file='EDSS_Progress.pdf')
# p1<-ggplot(Data1, aes(x=as.character(Group_subtype),fill=as.factor(EDSS_Progress)))+
#   geom_bar(stat="count",width=0.5,position='dodge')+
#   #scale_fill_manual(values=c('#304156','#e69f00'))+
#   geom_text(stat='count',aes(label=..count..), color="black", 
#             size=5,position=position_dodge(0.5),vjust=-0.5 )+
#   ggtitle("堆砌柱状图") +
#   xlab("Subtype")+
#   labs(fill="EDSS_Progress")+
#   theme_minimal()
# 
# print(p1)
# dev.off()
# 
# 
# # #SPMS_conversion
# # Data1<-Data
# # Data1<-Data1[!(is.na(Data1$SPMS_conversion)|Data1$SPMS_conversion=="[]"|Data1$SPMS_conversion==""|Data1$SPMS_conversion=="NaN"),]
# # Group_subtype1<-factor(Data1$Group_subtype)
# # Data1$SPMS_conversion<-as.numeric(Data1$SPMS_conversion)
# # 
# # tab12<-table(Data1$Group_subtype,Data1$SPMS_conversion)
# # 
# # 
# # # df <- as.data.frame(tab12)
# # # 
# # # colnames(df)<-c('Subtype','Sex','Freq')
# # # 
# # # ggbarplot(df, x = "Subtype", y = "Freq",
# # #           color = "Sex", position = position_dodge(),
# # #           #palette = c("brown", "blue", "black", "red")
# # #           )
# # pdf(file='SPMS_conversion.pdf')
# # p1<-ggplot(Data1, aes(x=as.character(Group_subtype),fill=as.factor(SPMS_conversion)))+
# #   geom_bar(stat="count",width=0.5,position='dodge')+
# #   #scale_fill_manual(values=c('#304156','#e69f00'))+
# #   geom_text(stat='count',aes(label=..count..), color="black", 
# #             size=5,position=position_dodge(0.5),vjust=-0.5 )+
# #   ggtitle("堆砌柱状图") +
# #   xlab("Subtype")+
# #   labs(fill="SPMS_conversion")+
# #   theme_minimal()
# # 
# # print(p1)
# # dev.off()
# # 
# 
# #生存曲线分析
# Data1<-Data
# Data1<-Data1[!(is.na(Data1$EDSS_Progress)|Data1$EDSS_Progress=="[]"|Data1$EDSS_Progress==""|Data1$EDSS_Progress=="NaN"),]
# #Data1$ml_stage<-as.factor(Data1$ml_stage)
# Data1$Sex_F1_M0<-as.factor(Data1$Sex_F1_M0)
# 
# res<-pairwise_survdiff(Surv(Followup_time_years,EDSS_Progress)~Group_subtype,data = Data1)
# 
# print(res)
# 
# fit<-survfit(Surv(Data1$Followup_time_years,Data1$EDSS_Progress)~Data1$Group_subtype,data=Data1)
# 
# pdf(file='EDSS_Progress_Surrival_curve.pdf')
# p1<-ggsurvplot(fit,pval=TRUE,conf.int = FALSE,risk.table = TRUE,risk.table.y.text.col=TRUE,risk.table.col="strata",)
# print(p1)
# 
# dev.off()
# 
# pdf(file='EDSS_Progress_Cumulative_hazard.pdf')
# p1<-ggsurvplot(fit,pval=TRUE,conf.int = FALSE,risk.table = TRUE,risk.table.y.text.col=TRUE,risk.table.col="strata",fun="cumhaz")
# print(p1)
# dev.off()
# 
# 
# #res.cox<-coxph(formula = Surv(Followup_time_years,EDSS_Progress) ~ Group_subtype+Age+Sex_F1_M0+Duration_Month, data = Data1)
# 
# res.cox<-coxph(formula = Surv(Followup_time_years,EDSS_Progress) ~ Group_subtype+ml_stage+Age+Sex_F1_M0, data = Data1)
# #res.cox<-coxph(formula = Surv(Followup_time_years,EDSS_Progress) ~ Group_subtype+Age, data = Data1)
# 
# summary(res.cox)
# 
# pdf(file='EDSS_Progress_ggforest_cox.pdf')
# ggforest(res.cox, data = Data1,  #数据集
#          main = 'Hazard ratio',  #标题
#          #cpositions = c(0.05, 0.15, 0.35),  #前三列距离
#          fontsize = 1, #字体大小
#          noDigits = 3 #保留HR值以及95%CI的小数位数
# )
# dev.off()
# 
# # ggsurvplot(data=Data1,survfit(res.cox),color="#2E9FDF",ggtheme=theme_minimal())
# # 
# # subtype_df<-within(Data1,{Group_subtype<-factor(Group_subtype,labels=c('A0','A1','A2'))})
# # 
# # res.cox<-coxph(formula = Surv(Followup_time_years,EDSS_Progress) ~ Group_subtype+Age+Sex._F1_M0, data = subtype_df)
# # 
# # 
# # ggforest(res.cox, data = Data1,  #数据集
# #          main = 'Hazard ratio',  #标题
# #          #cpositions = c(0.05, 0.15, 0.35),  #前三列距离
# #          fontsize = 1, #字体大小
# #          noDigits = 3 #保留HR值以及95%CI的小数位数
# # )
# # 
# #ggsurvplot(res.cox,subtype_df,conf.int=TRUE,legend.labs=c("Subtype=1","Subtype=2","Subtype=3"),ggtheme=theme_minimal())
# 
# 
# #生存曲线分析
# Data1<-Data
# Data1<-Data1[!(is.na(Data1$Followup_relapse)|Data1$Followup_relapse=="[]"|Data1$Followup_relapse==""|Data1$Followup_relapse=="NaN"),]
# 
# Data1[Data1$Followup_relapse==0,'Followup_relapse']=0
# Data1[Data1$Followup_relapse>0,'Followup_relapse']=1
# 
# #Data1$ml_stage<-as.factor(Data1$ml_stage)
# Data1$Sex_F1_M0<-as.factor(Data1$Sex_F1_M0)
# 
# res<-pairwise_survdiff(Surv(Followup_time_years,Followup_relapse)~Group_subtype,data = Data1)
# 
# print(res)
# 
# fit<-survfit(Surv(Data1$Followup_time_years,Data1$Followup_relapse)~Data1$Group_subtype,data=Data1)
# 
# pdf(file='Followup_relapse.pdf')
# p1<-ggsurvplot(fit,pval=TRUE,conf.int = FALSE,risk.table = TRUE,risk.table.y.text.col=TRUE,risk.table.col="strata",)
# print(p1)
# 
# dev.off()
# 
# pdf(file='Followup_relapse_Cumulative_hazard.pdf')
# p1<-ggsurvplot(fit,pval=TRUE,conf.int = FALSE,risk.table = TRUE,risk.table.y.text.col=TRUE,risk.table.col="strata",fun="cumhaz")
# print(p1)
# dev.off()
# 
# 
# res.cox<-coxph(formula = Surv(Followup_time_years,Followup_relapse) ~ Group_subtype+ml_stage+Sex_F1_M0, data = Data1)
# 
# 
# #res.cox<-coxph(formula = Surv(Followup_time_years,EDSS_Progress) ~ Group_subtype+Age, data = Data1)
# 
# summary(res.cox)
# 
# pdf(file='Followup_relapse_ggforest_cox.pdf')
# ggforest(res.cox, data = Data1,  #数据集
#          main = 'Hazard ratio',  #标题
#          #cpositions = c(0.05, 0.15, 0.35),  #前三列距离
#          fontsize = 1, #字体大小
#          noDigits = 3 #保留HR值以及95%CI的小数位数
# )
# dev.off()
# 
# 
# 
# 
# Data1<-Data
# Data1<-Data1[!(is.na(Data1$SPMS_conversion)|Data1$SPMS_conversion=="[]"|Data1$SPMS_conversion==""|Data1$SPMS_conversion=="NaN"),]
# 
# #Data1$ml_stage<-as.factor(Data1$ml_stage)
# Data1$Sex_F1_M0<-as.factor(Data1$Sex_F1_M0)
# 
# res<-pairwise_survdiff(Surv(Followup_time_years,SPMS_conversion)~Group_subtype,data = Data1)
# 
# print(res)
# 
# fit<-survfit(Surv(Data1$Followup_time_years,Data1$SPMS_conversion)~Data1$Group_subtype,data=Data1)
# 
# pdf(file='SPMS_conversion_Surrival_curve.pdf')
# p1<-ggsurvplot(fit,pval=TRUE,conf.int = FALSE,risk.table = TRUE,risk.table.y.text.col=TRUE,risk.table.col="strata",)
# print(p1)
# 
# dev.off()
# 
# pdf(file='SPMS_conversion_Cumulative_hazard.pdf')
# p1<-ggsurvplot(fit,pval=TRUE,conf.int = FALSE,risk.table = TRUE,risk.table.y.text.col=TRUE,risk.table.col="strata",fun="cumhaz")
# print(p1)
# dev.off()
# 
# 
# #res.cox<-coxph(formula = Surv(Followup_time_years,SPMS_conversion) ~ Group_subtype+Age+Sex_F1_M0+Duration_Month+ml_stage, data = Data1)
# 
# 
# res.cox<-coxph(formula = Surv(Followup_time_years,SPMS_conversion) ~ Group_subtype+ml_stage+Age+Sex_F1_M0, data = Data1)
# 
# #res.cox<-coxph(formula = Surv(Followup_time_years,EDSS_Progress) ~ Group_subtype+Age, data = Data1)
# 
# summary(res.cox)
# 
# pdf(file='SPMS_conversion_ggforest_cox.pdf')
# ggforest(res.cox, data = Data1,  #数据集
#          main = 'Hazard ratio',  #标题
#          #cpositions = c(0.05, 0.15, 0.35),  #前三列距离
#          fontsize = 1, #字体大小
#          noDigits = 3 #保留HR值以及95%CI的小数位数
# )
# dev.off()
# 
# 
