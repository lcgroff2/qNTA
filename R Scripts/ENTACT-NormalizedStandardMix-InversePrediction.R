#lgroff
library(ggplot2)
library(readxl)
library(xlsx)
library(investr)
library(data.table)
library(gridExtra)
library(bestNormalize)
source('C:/Users/lgroff/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/R/R Scripts/bestNormalize-CustomFunctions.R')

datadir <- 'C:/Users/lgroff/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/Data/ENTACT Semi Quant/Standard Mixtures/JRS Tracer Analysis'
figdir <- paste0(datadir,'/Figures')
setwd(datadir)

#Positive Mode Analysis:
posdata_1 <- read_xlsx('Supp_Table1_Pos.xlsx')
#Filter out values that were dropped for regression:
drop_pos <- unique(posdata_1[!is.na(posdata_1$Drop_for_Reg),'Preferred_Name'])
#posfilter_1 <- posdata_1[is.na(posdata_1$Drop_for_Reg),]
#Generate response factors from ratio of I/C:
posdata_1$RF <- posdata_1$Normalized_Intensity/posdata_1$Concentration
chemunique_pos <- unique(posdata_1$New_ID) 
#Prediction Interval level:
alpha = 0.01 #0.01 for 99% pred. interval, 0.05 for 95%
level = 1-alpha #One-sided distribution
#load in data frame with JRS RFs, slopes, and intercepts:
posdata_2 <- read_xlsx('Supp_Table2_Pos.xlsx')
#Perform individual cal. curve calibrations given a y0 value from their data subset (random data point choice):
posdata_2$y0 <- NA
posdata_2$Known_Conc <- NA
posdata_2$Conc_CC <- NA
posdata_2$Conc_CC_Upper <- NA
lmlist_pos <- vector('list',length(chemunique_pos))
plotlist_pos <- lmlist_pos
predlist_pos <- lmlist_pos
setwd(figdir)
set.seed(16384)
tmpidx_pos<-sample(1:100,nrow(posdata_2),size=nrow(posdata_2))/100
for (i in 1:length(chemunique_pos)){
  tmpdf <- posdata_1[posdata_1$New_ID==chemunique_pos[[i]],]
  tmpdf <- tmpdf[complete.cases(tmpdf$Normalized_Intensity),]
  posdata_2$y0[[i]] <- tmpdf[[ceiling(tmpidx_pos[[i]]*nrow(tmpdf)),'Normalized_Intensity']]
  posdata_2$Known_Conc[[i]] <- tmpdf[[ceiling(tmpidx_pos[[i]]*nrow(tmpdf)),'Concentration']]
  if(posdata_2$Slope_Full[[i]]==0){
    posdata_2$Intercept_Full[[i]] <- NA
    posdata_2$RF[[i]] <- posdata_2$y0[[i]]/posdata_2$Known_Conc[[i]]
  } 
  #store regression model objects in a list for later use:
  lmlist_pos[[i]] <- lm(tmpdf,formula = log10(Normalized_Intensity)~log10(Concentration))
  
  #store concentration predictions for individual chemical cal curves:
  posdata_2$Conc_CC[[i]] <- 10^((log10(posdata_2$y0[[i]])-lmlist_pos[[i]]$coefficients[1])/lmlist_pos[[i]]$coefficients[2])
  
  #store 99% prediction interval data for each chemical cal curve:
  predlist_pos[[i]] <- predict(lmlist_pos[[i]],
                           level=0.99,
                           interval='prediction')
  #run calibrate on n > 3 intensities to get 99% upper bound:
  tmpgg <- ggplot(tmpdf,aes(x=log10(Concentration),y=log10(Normalized_Intensity)))+
    geom_point()+
    geom_smooth(method='lm',formula=y~x,se=FALSE)+
  labs(x='Log10 Concentration',
       y= 'Log10 Normalized Intensity',
       title = unique(tmpdf$Preferred_Name))
  if (length(tmpdf$Intensity) > 3){
    tmpcal<-try(calibrate(lmlist_pos[[i]],
                      y0=log10(posdata_2$y0[[i]]),
                      level=0.99,
                      interval='inversion'),silent=T)
    #store individual regression plots in a list:
    tmpgg <- tmpgg+
      geom_ribbon(aes(ymin=predlist_pos[[i]][,'lwr'],
                      ymax=predlist_pos[[i]][,'upr']),
                  alpha=0.25)
    if (class(tmpcal) != 'try-error'){
      posdata_2$Conc_CC_Upper[[i]] <- 10^tmpcal$upper
      # posdata_2$Conc_CC[[i]]<-10^tmpcal$estimate
    }
  }
  plotlist_pos[[i]] <- tmpgg

  # ggsave(filename = paste0(unique(tmpdf$Preferred_Name),'_fit.png'),
  #        device='png',
  #        plot=tmpgg,
  #        dpi=300)
}
rm(tmpdf,tmpgg)

#Deal with zero slope cases, using RF = y0/injected conc instead:
for(i in 1:length(posdata_2$Slope_Full)){
  if(posdata_2$Slope_Full[[i]]==0){
    posdata_2$Intercept_Full[[i]] <- NA
    posdata_2$RF[[i]] <- posdata_2$y0[[i]]/posdata_2$Known_Conc[[i]]
  } 
}

# Conc_CC_plot_pos<-ggplot(posdata_2, 
#                          aes(x=Conc_CC_Upper))+
#   geom_histogram(fill='blue',color='black')+
#   scale_x_log10()+
#   labs(x='Log10Conc_CC_Upper',
#        y='Frequency',
#        title='ESI+, n = 74')

#Store Error Metrics:
posdata_2$Err_CCUppervsCCEst <- posdata_2$Conc_CC_Upper/posdata_2$Conc_CC

Conc_CC_Err_plot_pos<-ggplot(posdata_2,
                         aes(x=Err_CCUppervsCCEst))+
  geom_histogram(fill='blue',
                 color='black',
                 alpha=0.5,
                 na.rm=T)+
  scale_x_log10()+
  labs(x='log10 Conc_CC_0.99/Conc_CC',
       y='Frequency',
       title=paste0('Median Error = ',
                    round(median(posdata_2$Conc_CC_Upper/posdata_2$Conc_CC,na.rm=T),
                          digits=4)))

# Conc_CC_Err_plot<-grid.arrange(Conc_CC_plot_pos,Conc_CC_Err_plot_pos,nrow=1)

#Transform RF distribution and determine Pred. Interval bounds from transformed distribution:
posdata_2$RF_6thRoot <- posdata_2$RF^(1/6) #6th root transformation
#Upper prediction interval bound
pos_RF6th_upper<- mean(posdata_2$RF_6thRoot,na.rm=T)+
  qt(level,length(posdata_2$RF_6thRoot)-1)*sd(posdata_2$RF_6thRoot,na.rm=T)*sqrt(1+1/length(posdata_2$RF_6thRoot))
pos_RF6th_lower<- mean(posdata_2$RF_6thRoot,na.rm=T)-
  qt(1-0.05,length(posdata_2$RF_6thRoot)-1)*sd(posdata_2$RF_6thRoot,na.rm=T)*sqrt(1+1/length(posdata_2$RF_6thRoot))
pos_RF_plot<- ggplot(posdata_2,
                     aes(x=RF_6thRoot))+
  geom_histogram(fill='darkgreen',
                 color='black',alpha=0.5)+
  geom_vline(xintercept=pos_RF6th_lower,
             color='black', lwd=1.5, lty=2)+
  geom_vline(xintercept=pos_RF6th_upper,
             color='black',lwd=1.5,lty=2)+
  labs(x='Response Factor^(1/6)',
       y='Frequency',
       title='ESI+, n = 551')

#Calculate RF-based upper bound concentration estimates:
posdata_2$Conc_DefaultRF_Upper <- posdata_2$y0/pos_RF6th_lower^6
Conc_RF_plot_pos<-ggplot(posdata_2,aes(x=Conc_DefaultRF_Upper))+
  geom_histogram(fill='slateblue',
                 color='black',
                 alpha=0.5)+
  scale_x_log10()+
  labs(x='log10 Conc_RF_0.99',
       y='Frequency',
       title='ESI+, n = 551')

#Store Error Metrics:
posdata_2$Err_RFUppervsCCEst<- posdata_2$Conc_DefaultRF_Upper/posdata_2$Conc_CC
posdata_2$Err_RFUppervsCCUpper<- posdata_2$Conc_DefaultRF_Upper/posdata_2$Conc_CC_Upper


Conc_RF_Conc_CC_pos <- ggplot(posdata_2,aes(x=Err_RFUppervsCCEst))+
  geom_histogram(fill='slateblue',
                 color='black')+
  scale_x_log10()+
  labs(x='Upper Conc_RF/Conc_CC')+
  geom_vline(xintercept=1,
             lwd=1.5,lty=2,color='black')

Conc_RF_Conc_CC_Upper_pos <- ggplot(posdata_2,aes(x=Err_RFUppervsCCUpper))+
  geom_histogram(fill='slateblue',
                 color='black',
                 alpha=0.5)+
  scale_x_log10()+
  labs(x='Upper Conc_RF/Upper Conc_CC')+
  geom_vline(xintercept=1,
             lwd=1.5,lty=2,color='black')

# Conc_RF_Conc_CC_plot_pos<- grid.arrange(Conc_RF_Conc_CC_pos,Conc_RF_Conc_CC_Upper_pos,nrow=1)
#Calculate Default RF method concentrations with I/C RFs:
pos_RF6th_lower_IC <- mean(posdata_1$RF^(1/6),na.rm=T)-qt(level,df=length(posdata_1$RF)-1)*sd(posdata_1$RF^(1/6),na.rm=T)*sqrt(1+1/length(posdata_1$RF))
posdata_2$Conc_DefaultRF_Upper_IC <- posdata_2$y0/pos_RF6th_lower_IC^6

#Calculate log(Ionization Efficiency) vs. log(RF) RF bounds from Kruve methods using investr:
datadir2<-'C:/Users/lgroff/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/Data/ENTACT Semi Quant/Standard Mixtures/Anneli Updates'
posdata_anneli <- read_xlsx(paste0(datadir2,'/ENTACT_pos_predictions_Update_2021_03_15.xlsx'))
posdata2_smiles <- read_xlsx(paste0(datadir2,'/posdata2_DTXSID_MSReadySMILES.xlsx'))
posdata_2$MSReady_SMILES <- posdata2_smiles$MS_READY_SMILES
posdata_2$logIE_Pred <- NA
posdata_2$SMILES_Match <- NA
#Match on MS_Ready SMILES String and store logIE_pred:
for (i in 1:length(posdata_2$MSReady_SMILES)){
  if (sum(posdata_anneli$SMILES == unlist(strsplit(posdata_2$MSReady_SMILES[[i]],'[,]')))>0){
    posdata_2$logIE_Pred[[i]] <- as.numeric(unique(posdata_anneli[posdata_anneli$SMILES == unlist(strsplit(posdata_2$MSReady_SMILES[[i]],'[,]')),'logIE_pred'])[[1]][[1]])
    posdata_2$SMILES_Match[[i]] <- unique(posdata_anneli[posdata_anneli$SMILES == unlist(strsplit(posdata_2$MSReady_SMILES[[i]],'[,]')),'SMILES'])[[1]][[1]]
  }
  #for some salts with complex counterions (e.g. doxylamine succinate, 
  #where both pieces have a SMILES string comma-separated),
  #one segment of SMILES string will match if exists in Kruve predictions.
  #Split string, and find the matching piece. Store its logIE_pred:
  if (length(unlist(strsplit(posdata_2$MSReady_SMILES[[i]],'[,]'))) > 1){
    for (j in 1:length(unlist(strsplit(posdata_2$MSReady_SMILES[[i]],'[,]')))){
      if(sum(posdata_anneli$SMILES == unlist(strsplit(posdata_2$MSReady_SMILES[[i]],'[,]'))[[j]])>0){
        posdata_2$logIE_Pred[[i]] <- as.numeric(unique(posdata_anneli[posdata_anneli$SMILES == unlist(strsplit(posdata_2$MSReady_SMILES[[i]],'[,]'))[[j]],'logIE_pred'])[[1]][[1]])
        posdata_2$SMILES_Match[[i]] <- unique(posdata_anneli[posdata_anneli$SMILES == unlist(strsplit(posdata_2$MSReady_SMILES[[i]],'[,]'))[[j]],'SMILES'])[[1]][[1]]
      }
    }
  }
}
#Test other transforms for RF and IE:
posdata_2$IE_Pred <- 10^posdata_2$logIE_Pred
best_RF<-bestNormalize(posdata_2$RF,
                       standardize=F,
                       allow_orderNorm=F,
                       new_transforms= custom_transform,
                       loo=T)
posdata_2$bestN_RF<-best_RF$x.t
#Test for IE best normalization transformation:
best_IE<-bestNormalize(posdata_2$IE_Pred,
                       standardize=F,
                       allow_orderNorm=F,
                       new_transforms= custom_transform,
                       loo=T)
posdata_2$bestN_IE<-best_IE$x.t

# write.csv(posdata_2,file = 'ENTACTMixes_Kruve_logIEMatchesESI+.csv')
posdata_anneli[posdata_anneli$SMILES == unlist(strsplit(posdata_2[posdata_2$Preferred_Name == 'GW473178E methyl benzene sulphonic acid',]$MSReady_SMILES,'[,]'))[[2]],'SMILES']
anneli_lm_pos<-lm(posdata_2,
                  formula=posdata_2$RF_6thRoot~posdata_2$logIE_Pred)
anneli_predint_pos<-predict(anneli_lm_pos,
                            newdata=posdata_2,
                            interval='prediction',
                            level=0.99)
anneli_fig_pos<-ggplot(data=posdata_2,
                       aes(y=RF_6thRoot,x=logIE_Pred))+
  geom_point(na.rm=T)+
  geom_smooth(method='lm',
              formula=y~x,
              se=FALSE,
              na.rm=T)+
  geom_ribbon(aes(ymin=anneli_predint_pos[,'lwr'],
                  ymax=anneli_predint_pos[,'upr']),
              fill='darkgray',
              alpha=0.5)+
  labs(y=expression(Response~Factor^{(1/6)}),
       x=expression(log[10]~(Predicted~IE)))+
  # geom_vline(xintercept=pos_RF6th_lower,lwd=1.5,lty=2,color='red')+
  # geom_hline(yintercept=3.25,lwd=1.5,lty=2,color='red')+
  theme_classic()+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
    panel.border = element_rect(color='black',fill=NA,size=1.5),
    axis.ticks = element_line(size = .5),
    axis.ticks.length = unit(.2, "cm"))+
  # annotation_logticks(sides='l',outside=T, scaled=T)+
  coord_cartesian(clip='off')+
  annotation_logticks(base=2,sides='b',outside=T,
                      mid=unit(0.1,'cm'),long=unit(0.1,'cm'))

# ggsave(anneli_fig_pos,filename='RF6thRoot_vs_Log10IE.jpeg',device='jpeg',dpi=1200)

#determine threshold y0 above which log(IE) vs RF^(1/6) Inverse Estimation
#improves the prediction accuracy. Where pos_RF6th_lower intersects upper PI:
# tmplm<-lm(formula=anneli_predint_pos[,'upr']~posdata_2$RF_6thRoot)
# y0_threshold_pos<-tmplm$coefficients[2]*pos_RF6th_lower+tmplm$coefficients[1]

# ggsave('SI_Fig_log10IEvsRF6th_HiRes.jpeg',plot=anneli_fig_pos,dpi=1200,device='jpeg')
#save figure:
# ggsave(plot=anneli_fig_pos,
#        filename='lg2-logIEvs4thRF_99PredInt.png',
#        dpi=300,
#        device='png')
# qqnorm(anneli_lm$residuals,main='Log(IE) vs. RF^(1/6) Residuals\nNormal Q-Q Plot')
# qqline(anneli_lm$residuals,lwd=3,col='red')
#check regression residual linearity
# IERF_lm_pos <- lm(posdata_2,formula = logIE_Pred ~ RF_4thRoot)
# qqnorm(IERF_lm_pos$residuals,main = 'log IE vs log RF Residual\n QQ Plot')
# qqline(IERF_lm_pos$residuals,col='red',lwd=3)

posdata_2$logIE_RF_0.99 <- NA
posdata_2$logIE_RF_0.01 <- NA
posdata_2$logIE_RF_fit <- NA
#Run inverse predictions using the log IE and RF for each 
#chemical with the regression
lm_RFvsIE<-lm(posdata_2,formula = RF_6thRoot~logIE_Pred)
for (i in 1:length(posdata_2$RF_6thRoot)){
  if (!is.na(posdata_2$logIE_Pred[[i]])){
    pred_rf_pos<-predict.lm(lm_RFvsIE,
                       newdata=data.frame(logIE_Pred=posdata_2$logIE_Pred[[i]]),
                       interval='prediction',
                       level=0.99)
    posdata_2$logIE_RF_0.01[[i]]<-pred_rf_pos[,'lwr']
    posdata_2$logIE_RF_0.99[[i]]<-pred_rf_pos[,'upr']
    posdata_2$logIE_RF_fit[[i]]<-pred_rf_pos[,'fit']
  }
}
rm(tmpcal)
# write.xlsx(x=posdata_2,file='posdata_2_RFvsIE_SQ.xlsx')
#histogram of 1st percentile RF predictions from calibration prior to imputation:
RF_p0.01_pos_fig_pre<-ggplot(posdata_2,aes(x=logIE_RF_0.01),na.rm=T)+
  geom_histogram(fill='salmon',alpha=0.5,color='black',na.rm=T)+
  # scale_x_log10()+
  scale_y_log10()+
  geom_vline(xintercept=pos_RF6th_lower,lwd=1.5,lty=2)+
  # annotate(geom='text',x=55,y=45,label='tRF_0.01')+
  labs(x='1st Percentile RF^(1/6)',
       y='Frequency',
       title='ESI+, n = 551')

posdata_2$logIE_RF_0.01_IC <- posdata_2$logIE_RF_0.01  
#Impute 1st percentile RF_6thRoot for logIE_RF_0.01 values that are negative:
for (i in 1:length(posdata_2$logIE_RF_0.01)){
  if (!is.na(posdata_2$logIE_RF_0.01[[i]]) &
      posdata_2$logIE_RF_0.01[[i]] < pos_RF6th_lower){
    posdata_2$logIE_RF_0.01[[i]]<-pos_RF6th_lower
  }
  if (!is.na(posdata_2$logIE_RF_0.01_IC[[i]]) &
      posdata_2$logIE_RF_0.01_IC[[i]] < pos_RF6th_lower_IC){
    posdata_2$logIE_RF_0.01_IC[[i]]<-pos_RF6th_lower_IC
  }
}

#histogram of 1st percentile RF predictions from calibration post-imputation:
RF_p0.01_pos_fig_post<-ggplot(posdata_2,aes(x=logIE_RF_0.01),na.rm=T)+
  geom_histogram(fill='salmon',alpha=0.5,color='black',na.rm=T)+
  scale_x_log10()+
  scale_y_log10()

#examine normality of lower bound RF distribution
qqnorm(posdata_2$logIE_RF_0.01,main='QQ Plot\n1st percentile RF^(1/6) from calibration')
qqline(posdata_2$logIE_RF_0.01,col='red',lwd=3)
qqline(pos_RF6th_lower,col='blue',lwd=2,lty=2)
text(-2,11,'RF_0.01',col='blue')
#Calculate Kruve RF-based concentration estimates:
posdata_2$logIE_Conc_Upper <- posdata_2$y0/posdata_2$logIE_RF_0.01^6
posdata_2$logIE_Conc_Upper_IC <- posdata_2$y0/posdata_2$logIE_RF_0.01_IC^6

#Plot concentration histogram:
logIERF_conc_fig_pos <- ggplot(posdata_2,aes(x=logIE_Conc_Upper))+
  geom_histogram(fill='darkorange',alpha=0.5,color='black')+
  scale_x_log10()+
  labs(x='log10 Conc_0.99_IE',
       y='Frequency',
       title='ESI+, n = 551')

#Store Error Metrics:
posdata_2$Err_IEUppervsCCEst <- posdata_2$logIE_Conc_Upper/posdata_2$Conc_CC
posdata_2$Err_IEUppervsCCUpper <- posdata_2$logIE_Conc_Upper/posdata_2$Conc_CC_Upper
posdata_2$Err_RFUppervsCCEst_IC <- posdata_2$Conc_DefaultRF_Upper_IC/posdata_2$Conc_CC
posdata_2$Err_IEUppervsCCEst_IC <- posdata_2$logIE_Conc_Upper_IC/posdata_2$Conc_CC

#Compare to Pure RF distribution-based Concentration Predictions:
posdata_2$Err_RFUppervsIEUpper <- posdata_2$Conc_DefaultRF_Upper/posdata_2$logIE_Conc_Upper
RF_vs_IE_plot_pos<-ggplot(posdata_2,
                          aes(x=Err_RFUppervsIEUpper))+
  geom_histogram(fill='royalblue',alpha=0.5,color='black',na.rm=T)+
  scale_x_log10()+
  scale_y_log10()+
  labs(y='Frequency',
       x='Log10 Upper RF Conc/Upper IE Conc')+
  geom_vline(xintercept=1,
             lwd=1.5, lty=2)

upper_err_fig_pos<-ggplot(posdata_2)+
  geom_histogram(aes(x=PredError_DefaultRF_Upper,fill='Default RF PI'),
                 color='black',alpha=0.5)+
  geom_histogram(aes(x=PredError_logIERF_Upper,fill='Log(IE) vs Log(RF) PI'),
                 color='black',alpha=0.5)+
  scale_x_log10()+
  labs(x='Upper Bound Estimate / Upper Bound CC Estimate',
       y='Frequency',
       fill='Method:',
       title= 'ESI+ n = 74')+
  geom_vline(xintercept=1,
             lwd=1.5,
             lty=2)
  # theme(legend.position = 'top')

#Error Quotient Outlier identification:
# outlier_err_fig_pos<-anneli_fig_pos+
#   geom_point(data=posdata_2[posdata_2$Err_IEUppervsCCEst < 1,],
#              aes(x=logIE_Pred,y=RF_6thRoot),
#             color='red',size=2,na.rm=T)+
#   geom_hline(aes(yintercept=pos_RF6th_lower))
write.csv(posdata_2[,c("Err_RFUppervsCCEst","Err_IEUppervsCCEst",'Err_RFUppervsIEUpper')],'Pos_RFUppervsIEUpper.csv')
#save updated posdata_2
# write.xlsx(posdata_2,file='posdata_2_RFvsIE_SQ.xlsx')

#Add MS-Ready DTXCID and InChiKeys:
setwd('C:/Users/lgroff/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/Data/ENTACT Semi Quant/Standard Mixtures/Anneli Updates')
entact_merge_ids<-read.xlsx('ENTACT_Mixture_Key_Merge_Updated_020420.xlsx',sheetIndex = 1)
posdata_2$MS_Ready_DTXCID <- as.character(NA)
posdata_2$MS_Ready_InChiKey <- as.character(NA) 
for (i in 1:length(posdata_2$DTXSID)){
  posdata_2$MS_Ready_DTXCID[[i]] <- unique(entact_merge_ids[entact_merge_ids$DTXSID %in% posdata_2$DTXSID[[i]],"MS_Ready_DTXCID"])
  posdata_2$MS_Ready_InChiKey[[i]] <- unique(entact_merge_ids[entact_merge_ids$DTXSID %in% posdata_2$DTXSID[[i]],"MS_Ready_Inchikey"])
}

posdata_1$MS_Ready_DTXCID <- as.character(NA)
posdata_1$MS_Ready_InChiKey <- as.character(NA)
for (i in 1:length(posdata_1$DTXSID)){
  posdata_1$MS_Ready_DTXCID[[i]] <- unique(entact_merge_ids[entact_merge_ids$DTXSID %in% posdata_1$DTXSID[[i]],"MS_Ready_DTXCID"])
  posdata_1$MS_Ready_InChiKey[[i]] <- unique(entact_merge_ids[entact_merge_ids$DTXSID %in% posdata_1$DTXSID[[i]],"MS_Ready_Inchikey"])
}
# write.csv(posdata_1, file = 'posdata1updated.csv')
# write.csv(posdata_2, file = 'posdata2updated.csv')

pos_RFComparison<- ggplot()+
  geom_histogram(aes(x=posdata_1$RF^(1/6),
                     fill='I/C',
                     y=..density..),
                 color='black',
                 alpha=0.5)+
  geom_histogram(aes(x=posdata_2$RF^(1/6),
                     fill='Intercept',
                     y=..density..),
                 color='black',
                 alpha=0.5)+
  theme_classic()+
  geom_vline(aes(xintercept=pos_RF6th_lower_IC),
             color='salmon',lwd=1.5,lty=2)+
  geom_vline(aes(xintercept=pos_RF6th_lower),
             color='turquoise',lwd=1.5,lty=2)+
  labs(x=expression(Response~Factor^(1/6)),
       title='ESI+')

# ggsave(filename='ENTACT_ESI+_RFComparison.jpg',plot=pos_RFComparison,device='jpeg',dpi=1200)

#save posdata_2 dataframe to .csv:
#write.csv(posdata_2,file='Supp_Table_2_Pos_ConcPredUpdate.csv') 

#Save Figure:
# ggsave(plot=RF_vs_IE_plot_pos,
#        filename='lg2-Histogram_RFConcIEConcRatio.png',
#        device='png',
#        dpi=300)

#Compare to CNL Bayesian rjags calcs:
cnl_bayesint<-read.xlsx(paste0(datadir2,'/CNL_Bayes_PI_RF6thvsIE.xlsx'),
                       sheetIndex=1)
CNL_VS_LCG_PIs<- ggplot()+
  geom_point(data=posdata_2,
             aes(x=logIE_Pred,
                 y=RF_6thRoot))+
  geom_ribbon(data=cnl_bayesint,
              aes(x=IE_new,
                  ymin=uncertain_lower,
                  ymax=uncertain_upper,
                  fill='Bayesian'),
              alpha=0.5)+
  geom_ribbon(data=posdata_2,
              aes(x=logIE_Pred,
                  ymin=anneli_predint_pos[,'lwr'],
                  ymax=anneli_predint_pos[,'upr'],
                  fill='Frequentist'),
              alpha=0.5)+
  theme_classic()+
  labs(x=expression(log[10]~Ionization~Efficiency),
       y=expression(Resposne~Factor^(1/6)))
# ggsave(plot=CNL_VS_LCG_PIs,filename = 'Bayes_vs_Freq_PIs.jpeg',dpi=1200,device='jpeg')

#Repeat for Negative mode:
setwd(datadir)

#Negitive Mode Analysis:
negdata_1 <- read_xlsx('Supp_Table1_Neg.xlsx')
#Filter out values that were dropped for regression:
drop_neg <- unique(negdata_1[!is.na(negdata_1$Drop_for_Reg),'Preferred_Name'])
chemunique_neg <- unique(negdata_1$New_ID) 
#Prediction Interval level:
alpha = 0.01 #0.01 for 99% pred. interval, 0.05 for 95%
level = 1-alpha/2
#load in data frame with JRS RFs, slopes, and intercepts:
negdata_2 <- read_xlsx('Supp_Table2_Neg.xlsx')

#Perform individual cal. curve calibrations given a y0 value from their data subset (random data point choice):
negdata_2$y0 <- NA
negdata_2$Known_Conc <- NA
negdata_2$Conc_CC <- NA
negdata_2$Conc_CC_Upper <- NA
negdata_1$RF_IC <- negdata_1$Normalized_Intensity/negdata_1$Concentration
lmlist_neg <- vector('list',length(chemunique_neg))
plotlist_neg <- lmlist_neg
predlist_neg <- lmlist_neg
setwd(figdir)
set.seed(16384)
tmpidx_neg<-sample(1:100,nrow(negdata_2),size=nrow(negdata_2))/100
for (i in 1:length(chemunique_neg)){
  tmpdf <- negdata_1[negdata_1$New_ID==chemunique_neg[[i]],]
  tmpdf <- tmpdf[complete.cases(tmpdf$Normalized_Intensity),]
  negdata_2$y0[[i]] <- tmpdf[[ceiling(tmpidx_neg[[i]]*nrow(tmpdf)),'Normalized_Intensity']]
  negdata_2$Known_Conc[[i]] <- tmpdf[[ceiling(tmpidx_neg[[i]]*nrow(tmpdf)),'Concentration']]

  #store regression model objects in a list for later use:
  lmlist_neg[[i]] <- lm(tmpdf,formula = log10(Normalized_Intensity)~log10(Concentration))
  negdata_2$Slope_Full[[i]] <- lmlist_neg[[i]]$coefficients[[2]]
  negdata_2$Intercept_Full[[i]] <- lmlist_neg[[i]]$coefficients[[1]]
  #store concentration predictions for individual chemical cal curves:
  negdata_2$Conc_CC[[i]] <- 10^((log10(negdata_2$y0[[i]])-lmlist_neg[[i]]$coefficients[1])/lmlist_neg[[i]]$coefficients[2])
  
  #store 99% prediction interval data for each chemical cal curve:
  predlist_neg[[i]] <- predict(lmlist_neg[[i]],
                               level=0.99,
                               interval='prediction')
  
  #run calibrate on n > 3 intensities to get 99% upper bound:
  if (length(tmpdf$Intensity) > 3){
    tmpcal<-try(calibrate(lmlist_neg[[i]],
                          y0=log10(negdata_2$y0[[i]]),
                          level=0.99,
                          interval='inversion'),silent=T)
    if (class(tmpcal) != 'try-error'){
      negdata_2$Conc_CC_Upper[[i]] <- 10^tmpcal$upper
      # negdata_2$Conc_CC[[i]]<- 10^tmpcal$estimate
    }
  }
  #store individual regression plots in a list:
  tmpgg <- ggplot(tmpdf,aes(x=log10(Concentration),y=log10(Normalized_Intensity)))+
    geom_point()+
    geom_smooth(method='lm',formula=y~x,se=FALSE)+
    labs(x='Log10 Concentration',
         y= 'Log10 Normalized Intensity',
         title = unique(tmpdf$Preferred_Name))
  if (!is.na(predlist_neg[[i]][1,'lwr'])){
    tmpgg<- tmpgg+geom_ribbon(aes(ymin=predlist_neg[[i]][,'lwr'],
                                  ymax=predlist_neg[[i]][,'upr']),
                              alpha=0.25)
  }
  plotlist_neg[[i]] <- tmpgg
  #Deal with zero slope cases, using RF = y0/injected conc instead:
  if(nrow(tmpdf) == 1){
    negdata_2$Intercept_Full[[i]] <- NA
    negdata_2$RF[[i]] <- negdata_2$y0[[i]]/negdata_2$Known_Conc[[i]]
  }else{
   negdata_2$RF[[i]] <- 10^negdata_2$Intercept_Full[[i]]  
  }
}
  
  # ggsave(filename = paste0(unique(tmpdf$Preferred_Name),'fit.png'),
  #        device='png',
  #        plot=tmpgg,
  #        dpi=300)


Conc_CC_plot_neg<-ggplot(negdata_2, 
                         aes(x=Conc_CC_Upper))+
  geom_histogram(fill='blue',color='black')+
  scale_x_log10()+
  labs(x='Log10 Conc_CC_0.95',
       y='Frequency',
       title='ESI-, n = 5')

Conc_CC_Err_plot_neg<-ggplot(negdata_2,
                             aes(x=Conc_CC_Upper/Conc_CC))+
  geom_histogram(fill='blue',color='black',alpha=0.5,na.rm=T)+
  scale_x_log10()+
  labs(x='log10 Conc_CC_0.95/Conc_CC',
       y='Frequency',
       title=paste0('Median Error = ',round(median(negdata_2$Conc_CC_Upper/negdata_2$Conc_CC,na.rm=T),digits=4)))
Upper_Conc_CC_Err_plot_neg<-grid.arrange(Conc_CC_plot_neg,Conc_CC_Err_plot_neg,nrow=1)

#Transform RF distribution and determine Pred. Interval bounds from transformed distribution:
negdata_2$logRF <- log10(negdata_2$RF) #log transformation Normalizes ESI-
#Upper prediction interval bound
neg_logRF_upper<- mean(negdata_2$logRF)+
  qt(level,length(negdata_2$logRF)-1)*sd(negdata_2$logRF)*sqrt(1+1/length(negdata_2$logRF))
neg_logRF_lower<- mean(negdata_2$logRF)-
  qt(level,length(negdata_2$logRF)-1)*sd(negdata_2$logRF)*sqrt(1+1/length(negdata_2$logRF))
neg_logRF_IC_lower<-mean(log10(negdata_1$RF_IC),na.rm=T)-
  qt(level,length(negdata_1$RF_IC)-1)*sd(log10(negdata_1$RF_IC),na.rm=T)*sqrt(1+1/length(negdata_1$RF_IC))

neg_RF_plot<-ggplot(negdata_2,
       aes(x=logRF))+
  geom_histogram(fill='darkgreen',
                 color='black',alpha=0.5)+
  geom_vline(xintercept=neg_logRF_lower,
             color='black', lwd=1.5, lty=2)+
  geom_vline(xintercept=neg_logRF_upper,
             color='black',lwd=1.5,lty=2)+
  labs(x='log10 Response Factor',
       y='Frequency',
       title='ESI-, n = 269')
  # lims(x=c(0,180))

#Compare Intercept RFs in negdata_2 to I/C RFs in negdata_1:
RFComparison <- ggplot()+
  geom_histogram(aes(x=log10(negdata_1$RF_IC),
                     fill='I/C',
                     y=..density..),
                 color='black',
                 alpha=0.5)+
  geom_histogram(aes(x=log10(negdata_2$RF),
                     fill='Intercept',
                     y=..density..),
                 color='black',
                 alpha=0.5)+
  labs(x=expression(Log[10]~Response~Factor),
       title='ESI-')+
  theme_classic()+
  geom_vline(aes(xintercept=neg_logRF_IC_lower),
             color='salmon',lwd=1.5,lty=2)+
  geom_vline(aes(xintercept=neg_logRF_lower),
             color='turquoise',lwd=1.5,lty=2)

# ggsave(filename='ENTACT_ESI-_RFComparison.jpg',plot=RFComparison,device='jpeg',dpi=1200)

#Calculate RF-based upper bound concentration estimates:
negdata_2$Conc_DefaultRF_Upper <- negdata_2$y0/10^neg_logRF_lower
#Store Error Metrics:
negdata_2$Err_RFUppervsCCEst<- negdata_2$Conc_DefaultRF_Upper/negdata_2$Conc_CC
negdata_2$Err_RFUppervsCCUpper<- negdata_2$Conc_DefaultRF_Upper/negdata_2$Conc_CC_Upper

Conc_RF_Conc_CC_neg <- ggplot(negdata_2,aes(x=Err_RFUppervsCCEst))+
  geom_histogram(fill='slateblue',
                 color='black')+
  scale_x_log10()+
  labs(x='Upper Conc_RF/Conc_CC',
       title='ESI- n = 269')+
  geom_vline(xintercept=1,
             lwd=1.5,lty=2,color='black')

Conc_RF_Conc_CC_Upper_neg <- ggplot(negdata_2,aes(x=Err_RFUppervsCCUpper))+
  geom_histogram(fill='slateblue',
                 color='black',
                 alpha=0.5)+
  scale_x_log10()+
  labs(x='Upper Conc_RF/Upper Conc_CC',
       title='ESI-, n = 5')+
  geom_vline(xintercept=1,
             lwd=1.5,lty=2,color='black')+
  lims(x=c(1,300))

#Calculate log(Ionization Efficiency) vs. log(RF) RF bounds from Kruve methods using investr:
setwd('C:/Users/lgroff/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/Data/ENTACT Semi Quant/Standard Mixtures/Anneli Updates')
negdata_anneli <- read_xlsx('All_NegData_for_Evaluation.xlsx')
negdata_2$logIE_Pred <- NA
negdata_anneli <- negdata_anneli[negdata_anneli$DTXSID %in% negdata_2$DTXSID,]
for (i in 1:length(negdata_2$DTXSID)){
  if (sum(negdata_anneli$DTXSID == negdata_2$DTXSID[[i]])>0){
    negdata_2$logIE_Pred[[i]] <- as.numeric(unique(negdata_anneli[negdata_anneli$DTXSID == negdata_2$DTXSID[[i]],"Log_Pred_Ionization_Eff"])[[1]][[1]])
  }
}
anneli_lm_neg<-lm(negdata_2,
                  formula=negdata_2$logRF~negdata_2$logIE_Pred)
anneli_fig_neg<-ggplot(data=negdata_2,
                       aes(y=logRF,x=logIE_Pred))+
  geom_point(na.rm=T)+
  geom_smooth(method='lm',
              formula=y~x,
              se=FALSE,
              na.rm=T)
anneli_predint_neg<-predict(anneli_lm_neg,
                            newdata=negdata_2,
                            interval='prediction',
                            level=0.99)
anneli_fig_neg<-ggplot(data=negdata_2,
                       aes(y=logRF,x=logIE_Pred))+
  geom_point(na.rm=T)+
  geom_smooth(method='lm',
              formula=y~x,
              se=FALSE,
              na.rm=T)+
  geom_ribbon(aes(ymin=anneli_predint_neg[,'lwr'],
                  ymax=anneli_predint_neg[,'upr']),
              fill='darkgray',
              alpha=0.5)+
  labs(y=expression(log[10]~Response~Factor),
       x=expression(log[10]~(Predicted~IE)))+
  theme_classic()+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        panel.border = element_rect(color='black',fill=NA,size=1.5),
        axis.ticks = element_line(size = .5),
        axis.ticks.length = unit(.2, "cm"))+
  # annotation_logticks(sides='l',outside=T, scaled=T)+
  coord_cartesian(clip='off')+
  annotation_logticks(base=2,sides='b',outside=T,
                      mid=unit(0.1,'cm'),long=unit(0.1,'cm'))

# ggsave('SI_Fig_log10RFvslog10IE_ESI-_HiRes.jpeg',plot=anneli_fig_neg,dpi=1200,device='jpeg')
#check regression residual linearity
# IERF_lm_neg <- lm(negdata_2,formula = logIE_Pred ~ RF_4thRoot)
# qqnorm(IERF_lm_neg$residuals,main = 'log IE vs log RF Residual\n QQ Plot')
# qqline(IERF_lm_neg$residuals,col='red',lwd=3)

negdata_2$logIE_RF_0.99 <- NA
negdata_2$logIE_RF_0.01 <- NA
negdata_2$logIE_RF_fit <- NA
#Run inverse predictions using the log IE and RF for each 
#chemical with the regression
lm_RFvsIE_neg<-lm(negdata_2,formula = logRF~logIE_Pred)
for (i in 1:length(negdata_2$logRF)){
  if (!is.na(negdata_2$logIE_Pred[[i]])){
    pred_rf_neg<-predict.lm(lm_RFvsIE_neg,
                            newdata=data.frame(logIE_Pred=negdata_2$logIE_Pred[[i]]),
                            interval='prediction',
                            level=0.99)
    negdata_2$logIE_RF_0.01[[i]]<-pred_rf_neg[,'lwr']
    negdata_2$logIE_RF_0.99[[i]]<-pred_rf_neg[,'upr']
    negdata_2$logIE_RF_fit[[i]]<-pred_rf_neg[,'fit']
  }
}
rm(tmpcal)

#histogram of 1st percentile RF predictions from calibration prior to imputation:
RF_p0.01_neg_fig_pre<-ggplot(negdata_2,aes(x=logIE_RF_0.01),na.rm=T)+
  geom_histogram(fill='salmon',alpha=0.5,color='black',na.rm=T)+
  geom_vline(xintercept=neg_logRF_lower,lwd=1.5,lty=2)+
  annotate(geom='text',label='tRF_0.01',x=-15,y=45)+
  lims(y=c(0,45))+
  labs(x='5th Percentile log(RF)',
       y='Frequency',
       title='ESI-, n = 269')

#Impute 1st percentile RF_4thRoot for logIE_RF_0.01 values that are negative:
for (i in 1:length(negdata_2$logIE_RF_0.01)){
  if (!is.na(negdata_2$logIE_RF_0.01[[i]]) &
      negdata_2$logIE_RF_0.01[[i]] < neg_logRF_lower){
    negdata_2$logIE_RF_0.01[[i]]<-neg_logRF_lower
  }
}

#histogram of 1st percentile RF predictions from calibration negt-imputation:
RF_p0.01_neg_fig<-ggplot(negdata_2,aes(x=logIE_RF_0.01),na.rm=T)+
  geom_histogram(fill='salmon',alpha=0.5,color='black',na.rm=T)+
  scale_y_log10()

#examine normality of lower bound RF distribution
qqnorm(negdata_2$logIE_RF_0.01,main='QQ Plot\n1st percentile RF^(1/4) from calibration')
qqline(negdata_2$logIE_RF_0.01,col='red',lwd=3)

#Calculate Kruve RF-based concentration estimates:
negdata_2$logIE_Conc_Upper <-NA
negdata_2$logIE_Conc_Upper <- negdata_2$y0/10^negdata_2$logIE_RF_0.01

#Plot concentration histogram:
logIERF_conc_fig_neg <- ggplot(negdata_2,aes(x=logIE_Conc_Upper))+
  geom_histogram(fill='darkorange',alpha=0.5,color='black')+
  scale_x_log10()+
  labs(x='log10 Conc_0.95_IE',
       y='Frequency',
       title='ESI-, n = 269')
#Store Error Metrics:
negdata_2$Err_IEUppervsCCEst<- negdata_2$logIE_Conc_Upper/negdata_2$Conc_CC
negdata_2$Err_IEUppervsCCUpper<- negdata_2$logIE_Conc_Upper/negdata_2$Conc_CC_Upper

#Plot prediction error histogram:
Conc_IE_Conc_CC_neg <- ggplot(negdata_2,aes(Err_IEUppervsCCEst))+
  geom_histogram(fill='red',
                 color='black')+
  geom_vline(xintercept=1,
             lty=2, lwd=1.5, color='black')+
  scale_x_log10()+
  labs(x='Upper Conc_IE/Conc_CC',
       title='ESI-, n = 269')

Conc_IE_Conc_CC_Upper_neg <- ggplot(negdata_2,aes(Err_IEUppervsCCUpper))+
  geom_histogram(fill='red',
                 color='black',
                 alpha=0.5)+
  geom_vline(xintercept=1,
             lty=2, lwd=1.5, color='black')+
  scale_x_log10()+
  labs(x='Upper Conc_IE/Upper Conc_CC',
       title='ESI-, n = 5')

#Compare to Pure RF distribution-based Concentration Predictions:
negdata_2$RFConc_vs_IEConc <- negdata_2$Conc_DefaultRF_Upper/negdata_2$logIE_Conc_Upper
RF_vs_IE_plot_neg<-ggplot(negdata_2,
                          aes(x=RFErr_vs_IEErr))+
  geom_histogram(fill='royalblue',alpha=0.5,color='black',na.rm=T)+
  scale_x_log10()+
  scale_y_log10()+
  labs(y='Frequency',
       x='Log10 Upper RF Conc/Upper IE Conc')+
  geom_vline(xintercept=1,
             lwd=1.5, lty=2)
write.xlsx(negdata_2,file = 'negdata_2_RF_IE_flip.xlsx')
upper_err_fig_neg<-ggplot(negdata_2)+
  geom_histogram(aes(x=PredError_DefaultRF_Upper,fill='Default RF PI'),
                 color='black',alpha=0.5)+
  geom_histogram(aes(x=PredError_logIERF_Upper,fill='Log(IE) vs Log(RF) PI'),
                 color='black',alpha=0.5)+
  scale_x_log10()+
  labs(x='Upper Bound Estimate / Upper Bound CC Estimate',
       y='Frequency',
       fill='Method:',
       title= 'ESI- n = 5')+
  geom_vline(xintercept=1,
             lwd=1.5,
             lty=2)
# theme(legend.position = 'top')

negdata_2$MS_Ready_DTXCID <- as.character(NA)
negdata_2$MS_Ready_InChiKey <- as.character(NA) 
for (i in 1:length(negdata_2$DTXSID)){
  negdata_2$MS_Ready_DTXCID[[i]] <- unique(entact_merge_ids[entact_merge_ids$DTXSID %in% negdata_2$DTXSID[[i]],"MS_Ready_DTXCID"])
  negdata_2$MS_Ready_InChiKey[[i]] <- unique(entact_merge_ids[entact_merge_ids$DTXSID %in% negdata_2$DTXSID[[i]],"MS_Ready_Inchikey"])
}

negdata_1$MS_Ready_DTXCID <- as.character(NA)
negdata_1$MS_Ready_InChiKey <- as.character(NA)
for (i in 1:length(negdata_1$DTXSID)){
  negdata_1$MS_Ready_DTXCID[[i]] <- unique(entact_merge_ids[entact_merge_ids$DTXSID %in% negdata_1$DTXSID[[i]],"MS_Ready_DTXCID"])
  negdata_1$MS_Ready_InChiKey[[i]] <- unique(entact_merge_ids[entact_merge_ids$DTXSID %in% negdata_1$DTXSID[[i]],"MS_Ready_Inchikey"])
}
negdata_2$MS_Ready_DTXCID <- as.character(NA)
negdata_2$MS_Ready_InChiKey <- as.character(NA) 
for (i in 1:length(negdata_2$DTXSID)){
  negdata_2$MS_Ready_DTXCID[[i]] <- unique(entact_merge_ids[entact_merge_ids$DTXSID %in% negdata_2$DTXSID[[i]],"MS_Ready_DTXCID"])
  negdata_2$MS_Ready_InChiKey[[i]] <- unique(entact_merge_ids[entact_merge_ids$DTXSID %in% negdata_2$DTXSID[[i]],"MS_Ready_Inchikey"])
}
# #Error Quotient Outlier identification:
# outlier_err_fig_neg<-anneli_fig_neg+
#   geom_point(data=posdata_2[(posdata_2$logIE_Conc_Upper / posdata_2$Known_Conc) < 1,]$Preferred_Name,
#              aes(y=logIE_Pred,x=logRF),
#              color='red',size=2)

#Percentage of concentration estimates that yield improvement 
#in accuracy versus the default RF method:
#sum(!is.na(negdata_2$RFConc_vs_IEConc) & negdata_2$RFConc_vs_IEConc > 1)/sum(!is.na(negdata_2$RFConc_vs_IEConc))

#save negdata_2 dataframe to .csv:
#write.csv(negdata_2,file='Supp_Table_2_Neg_ConcPredUpdate.csv') 

#Save Figure:
# ggsave(plot=RF_vs_IE_plot_neg,
#        filename='lg2-Histogram_RFConcIEConcRatio.png',
#        device='png',
#        dpi=300)

negdata_2$Err_CCUppervsCCEst <- negdata_2$Conc_CC_Upper/negdata_2$Conc_CC
negdata_2$CASRN <- as.character(negdata_2$CASRN)
negdata_1$MS_Ready_DTXCID <- as.character(NA)
negdata_1$MS_Ready_InChiKey <- as.character(NA)
for (i in 1:length(negdata_1$DTXSID)){
  negdata_1$MS_Ready_DTXCID[[i]] <- unique(entact_merge_ids[entact_merge_ids$DTXSID %in% negdata_1$DTXSID[[i]],"MS_Ready_DTXCID"])
  negdata_1$MS_Ready_InChiKey[[i]] <- unique(entact_merge_ids[entact_merge_ids$DTXSID %in% negdata_1$DTXSID[[i]],"MS_Ready_Inchikey"])
}

write.xlsx(negdata_2,file = 'negdata_2_RFvsIE_SQ.xlsx')

#aggregate figures for ESI+/-:
Conc_CC_Upper_plot<-grid.arrange(Conc_CC_plot_pos,Conc_CC_plot_neg,nrow=1)
tRF_plot<-grid.arrange(pos_RF_plot,neg_RF_plot,nrow=1)
tRF_Conc_plot<-grid.arrange(Conc_RF_plot_pos,Conc_RF_plot_neg,nrow=1)
anneli_plot<-grid.arrange(anneli_fig_pos+labs(title='ESI+, 99% PI'),anneli_fig_neg,nrow=1)
logIERF_pre<-grid.arrange(RF_p0.01_pos_fig_pre,RF_p0.01_neg_fig_pre,nrow=1)
logIERF_conc_fig<-grid.arrange(logIERF_conc_fig_pos,logIERF_conc_fig_neg,nrow=1)
Conc_Method_CC_fig<-grid.arrange(upper_err_fig_pos,upper_err_fig_neg,nrow=1)
Upper_Conc_CC_Err_plot_pos<-grid.arrange(Conc_CC_plot_pos,Conc_CC_Err_plot_pos,nrow=1)
Upper_Conc_CC_Err_plot_neg<-grid.arrange(Conc_CC_plot_neg,Conc_CC_Err_plot_neg,nrow=1)

ViolinSummary_df<- melt(setDT(posdata_2[,c("Preferred_Name","Err_CCUppervsCCEst","Err_IEUppervsCCEst","Err_IEUppervsCCUpper",
                                           "Err_RFUppervsCCEst","Err_RFUppervsCCUpper","Err_RFUppervsIEUpper")]),
                        id.vars = c('Preferred_Name'),
                        measure.vars = c("Err_CCUppervsCCEst","Err_IEUppervsCCEst","Err_IEUppervsCCUpper",
                                         "Err_RFUppervsCCEst","Err_RFUppervsCCUpper","Err_RFUppervsIEUpper"),
                        variable.name = 'Method',
                        value.name = 'Error')
ViolinSummary_plot<- ggplot(ViolinSummary_df,
                            aes(x=Error,y=Method),na.rm=T)+
  # geom_boxplot(width=0.25)+
  geom_violin(fill=NA,scale='width',)+
  scale_x_log10()+
  labs(y='Method Comparison',x=' Concentration Estimate Ratio\nMethod 1 / Method 2')+
  scale_y_discrete(labels=c('CC_Upper / CC_Est', 'IE_Upper / CC_Est', 'IE_Upper / CC_Upper',
                            'RF_Upper / CC_Est','RF_Upper / CC_Upper','RF_Upper / IE_Upper'))

ViolinSummaryPlot2 <- ggplot(ViolinSummary_df[ViolinSummary_df$Method==c("Err_CCUppervsCCEst",
                                                                         "Err_IEUppervsCCEst",
                                                                         "Err_RFUppervsCCEst")],
                             aes(x=Error,y=Method),na.rm=T)+
  # geom_boxplot(width=0.25)+
  geom_violin(fill=NA,lwd=0.5,scale='width')+
  scale_x_log10()+
  labs(x='Error Ratio (Method 1 / Method 2)',y='')+
  scale_y_discrete(labels=c('CC_Upper / CC_Est','IE_Upper / CC_Est','RF_Upper / CC_Est'))+
  geom_jitter(alpha=0.5,color='royalblue')+
  geom_boxplot(width=0.1,lwd=1)+
  theme_classic()

anneli_fig_pos+
  geom_point(data=posdata_2[posdata_2$Err_CCUppervsCCEst > 6,],
             aes(y=logIE_Pred,x=RF_6thRoot),
             color='red',size=2)
pos_extrema<-posdata_2[posdata_2$Err_CCUppervsCCEst > 100,]