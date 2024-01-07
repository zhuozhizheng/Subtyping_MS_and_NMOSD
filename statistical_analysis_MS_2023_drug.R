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
library(eoffice)


savepath<-"D:/ZZZ/Manuscripts/Subtyping/Sustain_results/Results/MS_0411"
setwd(savepath)
Data0<-read.csv('D:/ZZZ/Manuscripts/Subtyping/Sustain_results/MS_AQP4Pos_NMO/WorkshopOutput_1212_ms/newMStest.csv',header = TRUE)
rownames(Data0)<-Data0$PatientID_Image


#read MRI metrics 
data_all1<-read.csv('D:/ZZZ/Manuscripts/Subtyping/Sustain_results/New_Clinical_Validation_information1.csv',header = TRUE);

data_all1<-data_all1[!is.na(data_all1$Group)&!is.na(data_all1$Sex_F1_M0)&!is.na(data_all1$Age),]

data_all2<-read.csv('D:/ZZZ/Manuscripts/Subtyping/Sustain_results/tem_clinical_data_Volume_20221119_complete_drug_follow_up.csv',header = TRUE);
#select data
data_all2<-data_all2[data_all2$baseline=='baseline'&is.na(data_all2$Repeat_data)&!is.na(data_all2$Group)&
                       !is.na(data_all2$Sex_F1_M0)&!is.na(data_all2$Age),]

Nor_var<-c('PatientID_Image','Group','Age','Sex_F1_M0','Education_years','Duration_Month','relapses',
           'WMH_volume','WM.hypointensities.1',
           'CortexVol','SubCortGrayVol','CerebralWhiteMatterVol',
           'SUS_Brainstem','SUS_CerebellumGM','SUS_CerebellumWM','SUS_MUCCA123','EstimatedTotalIntraCranialVol',
           'choroid_volume','MMSE','MOCA','CVLT','BVLT','PASAT','EDSS','Dataset')
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
           'SUS_CerebellumGM','SUS_CerebellumWM','SUS_Brainstem','EstimatedTotalIntraCranialVol','sensitive_drug_DMT',
           'Followup_EDSS_drug','Followup_relapses_drug')

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
  

setwd(paste0(savepath,"/Results/MS"))
Data<-Data[!(is.na(Data$ml_subtype)|Data$ml_subtype==""),]
#Data<-Data[Data$ml_stage>0,] #remove stage 0

Data[Data$Group=='MS'&Data$ml_stage==0,"Group_subtype"]='MS-NA';
Data[Data$Group=='MS'&Data$ml_subtype==0&Data$ml_stage>0,"Group_subtype"]='MS-SC';
Data[Data$Group=='MS'&Data$ml_subtype==1&Data$ml_stage>0,"Group_subtype"]='MS-DGM';
Data[Data$Group=='MS'&Data$ml_subtype==2&Data$ml_stage>0,"Group_subtype"]='MS-C';



Data$Group_subtype<-factor(Data$Group_subtype,levels=c('MS-NA','MS-C','MS-SC','MS-DGM'))


Data<-Data[Data$Group=='MS',]


Data_relapse_drug<-Data[!is.na(Data$Followup_relapses_drug),]

tem_data<-Data_relapse_drug[,c("Group_subtype",'sensitive_drug_DMT','relapses','Followup_relapses_drug')];

n1=nrow(tem_data)
tem_data<-rbind(tem_data,tem_data);
n2=nrow(tem_data)

tem_data[(n1+1):n2,'relapses']=tem_data[(n1+1):n2,'Followup_relapses_drug']
tem_data[(n1+1):n2,'T']<-c('1year');
tem_data[1:n1,'T']<-c('Baseline');

tem_data[1:n1,'Sub']<-paste0('Sub',c(1:n1));
tem_data[(n1+1):n2,'Sub']<-paste0('Sub',c(1:n1));

tem_data$sensitive_drug_DMT<-factor(tem_data$sensitive_drug_DMT,
                                 levels=c('None','IFN_β','TF','RTX','Fingolimod','siponimod'))

tem_data$T<-factor(tem_data$T,levels=c('Baseline','1year'))

ggplot(tem_data, aes(x = T, y = relapses)) +
  facet_wrap("~ Group_subtype",nrow = 1, ncol = 4)+
  geom_point(size = 5) + 
  geom_point(size = 5,shape=21) +  #绘制散点
  geom_line(aes(group = Sub,color=sensitive_drug_DMT), lwd = 1)+
  theme_bw()
  


Data_EDSS_drug<-Data[!is.na(Data$Followup_EDSS_drug),]

tem_data<-Data_EDSS_drug[,c("Group_subtype",'sensitive_drug_DMT','EDSS','Followup_EDSS_drug')];

n1=nrow(tem_data)
tem_data<-rbind(tem_data,tem_data);
n2=nrow(tem_data)

tem_data[(n1+1):n2,'EDSS']=tem_data[(n1+1):n2,'Followup_EDSS_drug']
tem_data[(n1+1):n2,'T']<-c('1year');
tem_data[1:n1,'T']<-c('Baseline');

tem_data[1:n1,'Sub']<-paste0('Sub',c(1:n1));
tem_data[(n1+1):n2,'Sub']<-paste0('Sub',c(1:n1));

tem_data$sensitive_drug_DMT<-factor(tem_data$sensitive_drug_DMT,
                                    levels=c('None','IFN_β','TF','RTX','Fingolimod','siponimod'))

tem_data$T<-factor(tem_data$T,levels=c('Baseline','1year'))

ggplot(tem_data, aes(x = T, y = EDSS)) +
  facet_wrap("~ Group_subtype",nrow = 1, ncol = 4)+
  geom_point(size = 5) + 
  geom_point(size = 5,shape=21) +  #绘制散点
  geom_line(aes(group = Sub,color=sensitive_drug_DMT), lwd = 1)+
  theme_bw()


# 
# 
# 
# Data$sensitive_drug_DMT<-as.factor(Data$sensitive_drug_DMT)
# 
# Data1<-Data;
# Data1<-Data1[!is.na(Data1$sensitive_drug_DMT)&Data1$sensitive_drug_DMT!='',]
# 
# tem_var<-unique(Data1$sensitive_drug_DMT)
# 
# Data1$sensitive_drug_DMT<-factor(Data1$sensitive_drug_DMT,
#                                  levels=c('None','IFN_β','TF','RTX','Fingolimod','siponimod'))
# 
# 
# ggplot(Data1, aes(y = Group_subtype)) +
#   geom_bar(aes(fill = sensitive_drug_DMT), position = position_stack(reverse = TRUE)) +
#   theme(legend.position = "top")
# 
# pdf('Drugs.pdf',width=5,height=3)
# ggplot(Data1, aes(x = Group_subtype)) +
#   geom_bar(aes(fill = sensitive_drug_DMT), position = 'fill',width=0.5) +
#   theme(legend.position = "top")+theme_bw()
# dev.off()
# 
# library(data.table)
# ka<-xtabs(~Data1$Group_subtype+Data1$sensitive_drug_DMT,data=Data1)
# 
# ka<-xtabs(~Data1$sensitive_drug_DMT+Data1$Group_subtype,data=Data1)
# 
# chisq.test(ka)
# 
# 
# Data$Sex_F1_M0<-factor(Data$Sex_F1_M0,levels=c(0,1));
# 
# add_var<-c('SDMT','COWAT','Beton',
#            'Followup_relapse','Followup_EDSS','Followup_time_years',
#            'EDSS_Progress','SPMS_conversion',
#            
#            'T1T2Ratio_cerebral_GM','T1T2Ratio_cerebral_WM','T1T2Ratio_cerebellum_GM','T1T2Ratio_cerebellum_WM',
#            'T1T2Ratio_brainstem','FA_cerebral_GM','FA_cerebral_WM','FA_cerebellum_GM','FA_cerebellum_WM','FA_brainstem',
#            'MD_cerebral_GM','MD_cerebral_WM','MD_cerebellum_GM','MD_cerebellum_WM','MD_brainstem',
#            'DC_Cerebral_GM','DC_Cerebellum_GM',
#            'fALFF_Cerebral_GM','fALFF_Cerebellum_GM',
#            'NFL','GFAP','AB40','AB42','TAU','TREM2',
#            'SUS_MUCCA123','CortexVol','SubCortGrayVol','CerebralWhiteMatterVol',
#            'SUS_CerebellumGM','SUS_CerebellumWM','SUS_Brainstem','EstimatedTotalIntraCranialVol')
# library(table1)
# library(boot)
# 
# table1(~ Age+factor(Sex_F1_M0)+Education_years+Duration_Month+relapses+EDSS+Followup_EDSS+Followup_time_years+Followup_relapse+
#          
#          factor(EDSS_Progress)+factor(SPMS_conversion)+
#          
#          MMSE+MOCA+SDMT+CVLT+BVLT+PASAT+COWAT+Beton+
#          
#          #T1T2Ratio_cerebral_GM+T1T2Ratio_cerebral_WM+T1T2Ratio_cerebellum_GM+T1T2Ratio_cerebellum_WM+T1T2Ratio_brainstem+
#          
#          
#          FA_cerebral_GM+FA_cerebral_WM+FA_cerebellum_GM+FA_cerebellum_WM+FA_brainstem+
#          
#          MD_cerebral_GM+MD_cerebral_WM+MD_cerebellum_GM+MD_cerebellum_WM+MD_brainstem+
#          
#          DC_Cerebral_GM+DC_Cerebellum_GM+fALFF_Cerebral_GM+fALFF_Cerebellum_GM+
#          
#          NFL+GFAP+AB40+AB42+TAU+TREM2+
#          
#          choroid_volume+WMH_volume+WM.hypointensities.1+EstimatedTotalIntraCranialVol+
#          
#          SUS_MUCCA123+CortexVol+SubCortGrayVol+CerebralWhiteMatterVol+SUS_CerebellumGM+SUS_CerebellumWM+SUS_Brainstem|Group_subtype,
#          
#          data=Data,topclass='Rtable1-grid Rtable1-zebra')
# 
# 
# 
# 
# 
# setwd(savepath)
# str<-c('Age','Education_years','Duration_Month','relapses','EDSS',
#        'Followup_time_years','Followup_relapse','Followup_EDSS',
#       
#       'MMSE','MOCA','SDMT','CVLT','BVLT','PASAT','COWAT','Beton',
#       
#       'WMH_volume','WM.hypointensities.1','choroid_volume','EstimatedTotalIntraCranialVol',
#       'CortexVol','SubCortGrayVol','CerebralWhiteMatterVol','SUS_CerebellumGM',
#       'SUS_CerebellumWM','SUS_Brainstem','SUS_MUCCA123',
#       
#       'FA_cerebral_WM','FA_cerebellum_WM','FA_brainstem',
#       'MD_cerebral_WM','MD_cerebellum_WM','MD_brainstem',
#       'DC_Cerebral_GM','DC_Cerebellum_GM','fALFF_Cerebral_GM','fALFF_Cerebellum_GM')
# 
# Results<-data.frame(matrix(NA,length(str),6))
# 
# for(i in 1:length(str))
# 
# { 
# 
#   if(i>3)   
# {tem<-Data[,c(str[i],'Group_subtype','Age','Sex_F1_M0','Dataset')];
# 
# colnames(tem)<-c('tem_var','Group_subtype','Age','Sex_F1_M0','Dataset');
# 
# tem<-tem[!is.na(tem$tem_var),];
# 
# #fit<-lm(tem_var~Age+Sex_F1_M0+Dataset,data=tem)
# 
# #tem$tem_var<-fit$residuals+mean(fit$fitted.values)
# 
# res1<-kruskal.test(tem_var~Group_subtype,data=tem)
# res2<-with(data=tem,
#      pairwise.wilcox.test(x=tem_var,g=Group_subtype,p.adjust.method = "fdr")
# )
# 
# Results[i,1]<-res1$statistic;
# Results[i,2]<-res1$p.value;
# Results[i,3]<-res2$p.value[1,1];
# Results[i,4]<-res2$p.value[2,1];
# Results[i,5]<-res2$p.value[3,1];
# Results[i,6]<-str[i];
# Results[i,7]<-res2$p.value[2,2];
# Results[i,8]<-res2$p.value[3,2];
# Results[i,9]<-res2$p.value[3,3];
# 
# 
#   }
# 
#   if(i<=3)   
#   {tem<-Data[,c(str[i],'Group_subtype')];
#   
#   colnames(tem)<-c('tem_var','Group_subtype');
#   
#   res1<-kruskal.test(tem_var~Group_subtype,data=tem)
#   res2<-with(data=tem,
#        pairwise.wilcox.test(x=tem_var,g=Group_subtype,p.adjust.method = "fdr")
#   )
#   
#   Results[i,1]<-res1$statistic;
#   Results[i,2]<-res1$p.value;
#   Results[i,3]<-res2$p.value[1,1];
#   Results[i,4]<-res2$p.value[2,1];
#   Results[i,5]<-res2$p.value[3,1];
#   Results[i,6]<-str[i];
#   Results[i,7]<-res2$p.value[2,2];
#   Results[i,8]<-res2$p.value[3,2];
#   Results[i,9]<-res2$p.value[3,3];
#   
#   }
#   
# }
# 
# 
# colnames(Results)<-c('Chi-square','p','pS1','ps2','ps3','var','ps2-s1','ps3-s1','ps3-s2')
# #write.csv(Results,'With_NA_as_reference.csv')
# 
# write.csv(Results,'With_NA_as_reference_fdr.csv')
# 
# #one vs all others
# 
# G=c('MS-C','MS-SC','MS-DGM')
# 
# 
# for(gg in G)
# {
# tem_Data<-Data[Data$Group_subtype!='MS-NA',];
# 
# tem_Data$Group_subtype<-as.character(tem_Data$Group_subtype)
# 
# tem_Data[tem_Data$Group_subtype!=gg,'Group_subtype']<-c('others')
# 
# tem_Data$Group_subtype<-factor(tem_Data$Group_subtype,levels=c('others',gg))
# 
# Results<-data.frame(matrix(NA,length(str),6))
# 
# for(i in 1:length(str))
#   
# { 
#   
#   if(i>3)   
#   {tem<-tem_Data[,c(str[i],'Group_subtype','Age','Sex_F1_M0','Dataset')];
#   
#   colnames(tem)<-c('tem_var','Group_subtype','Age','Sex_F1_M0','Dataset');
#   
#   tem<-tem[!is.na(tem$tem_var),];
#   
#   #fit<-lm(tem_var~Age+Sex_F1_M0+Dataset,data=tem)
#   
#   #tem$tem_var<-fit$residuals+mean(fit$fitted.values)
#   
#   
#   res1<-wilcox.test(tem[tem$Group_subtype==gg,'tem_var'],tem[tem$Group_subtype=='others','tem_var'])
#   
#   
#   Results[i,1]<-res1$statistic
#   Results[i,2]<-res1$p.value;
#   Results[i,3]<-str[i];
#   
#   
#   }
#   
#   if(i<=3)   
#   {tem<-tem_Data[,c(str[i],'Group_subtype')];
#   
#   colnames(tem)<-c('tem_var','Group_subtype');
#   
#   res1<-wilcox.test(tem[tem$Group_subtype==gg,'tem_var'],tem[tem$Group_subtype=='others','tem_var'])
#   
#   
#   Results[i,1]<-res1$statistic
#   Results[i,2]<-res1$p.value;
#   Results[i,3]<-str[i];
#   
#   }
#   
# }
# 
# write.csv(Results,paste0(gg,'_With_others_as_reference.csv'))
# 
# }
# 
# 
# 
# my_comparisons<-list(c('MS-NA','MS-C'),
#                      c('MS-NA','MS-SC'),
#                      c('MS-NA','MS-DGM'),
#                    
#                      c('MS-C','MS-SC'),
#                      c('MS-C','MS-DGM'),
#                    
#                      c('MS-SC','MS-DGM'))
# 
# 
# my_comparisons1<-list(c('MS-C','MS-SC'),
#                       c('MS-C','MS-DGM'),
#                      
#                       c('MS-C','MS-DGM'))
# 
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
#   stat_compare_means(method='kruskal.test',label.y=0)
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
#   stat_compare_means(method='kruskal.test',label.y=0)
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
#   stat_compare_means(method='kruskal.test',label.y=-5)
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
#   stat_compare_means(method='kruskal.test',label.y=-5)
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
#   stat_compare_means(method='kruskal.test',label.y=-5)
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
#   stat_compare_means(method='kruskal.test',label.y=-5)
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
#   stat_compare_means(method='kruskal.test',label.y=-5)
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
#   stat_compare_means(method='kruskal.test',label.y=-5)
# print(p1)
# dev.off()
# }
# 
# 
# 
# library(ggpmisc)
# 
# LETTERS<-c('MS-NA','MS-C','MS-SC','MS-DGM')
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
# ressult<-summary(res.cox)
# 
# pdf(file='EDSS_Progress_ggforest_cox.pdf')
# p<-ggforest(res.cox, data = Data1,  #数据集
#          main = 'Hazard ratio',  #标题
#          #cpositions = c(0.05, 0.15, 0.35),  #前三列距离
#          fontsize = 1, #字体大小
#          noDigits = 3 #保留HR值以及95%CI的小数位数
# )
# print(p)
# dev.off()
# 
# write.csv(ressult$coefficients,'EDSS_Progress_ggforest_cox.csv')
# 
# topptx(figure = p,filename = "EDSS_Progress_ggforest_cox.pptx")
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
# res.cox<-coxph(formula = Surv(Followup_time_years,Followup_relapse) ~ Group_subtype+ml_stage+Age+Sex_F1_M0, data = Data1)
# 
# 
# #res.cox<-coxph(formula = Surv(Followup_time_years,EDSS_Progress) ~ Group_subtype+Age, data = Data1)
# 
# ressult<-summary(res.cox)
# 
# pdf(file='Followup_relapse_ggforest_cox.pdf')
# p<-ggforest(res.cox, data = Data1,  #数据集
#          main = 'Hazard ratio',  #标题
#          #cpositions = c(0.05, 0.15, 0.35),  #前三列距离
#          fontsize = 1, #字体大小
#          noDigits = 3 #保留HR值以及95%CI的小数位数
# )
# print(p)
# dev.off()
# 
# write.csv(ressult$coefficients,'Followup_relapse_ggforest_cox.csv')
# topptx(figure = p,filename = "Followup_relapse_ggforest_cox.pptx")
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
# ressult<-summary(res.cox)
# 
# pdf(file='SPMS_conversion_ggforest_cox.pdf')
# p<-ggforest(res.cox, data = Data1,  #数据集
#          main = 'Hazard ratio',  #标题
#          #cpositions = c(0.05, 0.15, 0.35),  #前三列距离
#          fontsize = 1, #字体大小
#          noDigits = 3 #保留HR值以及95%CI的小数位数
# )
# print(p)
# dev.off()
# 
# write.csv(ressult$coefficients,'SPMS_conversion_ggforest_cox.csv')
# 
# topptx(figure = p,filename = "SPMS_conversion_ggforest_cox.pptx")