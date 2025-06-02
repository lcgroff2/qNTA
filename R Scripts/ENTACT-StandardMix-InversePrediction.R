library(ggplot2)
library(readxl)
library(xlsx)
library(investr)
library(data.table)
library(gridExtra)

datadir<-'C:/Users/lgroff/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/Data/ENTACT Semi Quant/Standard Mixtures'
figdir<-paste0(datadir,'/Figures/Individual Linear Fits R')
setwd(datadir)
posdata <- read_xlsx('All_Neg_Data_Clean.xlsx')
posdata[(complete.cases(posdata$Questionable_Hit)==FALSE),'Questionable_Hit']<-0
#rows that are NA on DTXSID are NA on everything, remove those, questionable hits, and chemicals
#that don't have low conc data.
posdata <- posdata[(posdata$Questionable_Hit != 1) &
                   complete.cases(posdata$DTXSID) &
                   complete.cases(posdata$BlankSub_Med_Low),]
#Tack additional naming onto isomers:
for (i in 1:length(posdata$Isomer_Number)){
  if (posdata$DTXSID_Isomer[i] == 1){
    if (posdata$Isomer_Order[i] == 'First'){
      posdata$Preferred_Name[i] <- paste0(posdata$Preferred_Name[i],'_1')
      posdata$DTXSID[i] <- paste0(posdata$DTXSID[i],'_1')
    }else if (posdata$Isomer_Order[i] == 'Second'){
        posdata$Preferred_Name[i] <- paste0(posdata$Preferred_Name[i],'_2')
        posdata$DTXSID[i] <- paste0(posdata$DTXSID[i],'_2')
    }else if (posdata$Isomer_Order[i] == 'Third'){
          posdata$Preferred_Name[i] <- paste0(posdata$Preferred_Name[i],'_3')
          posdata$DTXSID[i] <- paste0(posdata$DTXSID[i],'_3')
    }
  }
}
#make long format concentration and intensity tables with columns for the levels, and concatenate them:
poslongconc <- melt(setDT(posdata[,c('DTXSID','Preferred_Name','n_Repeated_Hits',
                                     'Corrected_LowConc','Corrected_MidConc','Corrected_HighConc')]),
                    id.vars = c('DTXSID','Preferred_Name','n_Repeated_Hits'),
                    measure.vars = c('Corrected_LowConc','Corrected_MidConc','Corrected_HighConc'),
                    variable.name = 'Concentration_Level')
names(poslongconc)[length(poslongconc)] <- 'Concentration'
poslongresponse <- melt(setDT(posdata[,c('DTXSID','Preferred_Name',
                                         'BlankSub_Med_Low','BlankSub_Med_Mid','BlankSub_Med_High')]),
                        id.vars = c('DTXSID','Preferred_Name'),
                        measure.vars = c('BlankSub_Med_Low','BlankSub_Med_Mid','BlankSub_Med_High'),
                        variable.name = 'Intensity_Level')
names(poslongresponse)[length(poslongresponse)] <- 'Intensity'
poslong <- cbind(poslongconc[,c('DTXSID','Preferred_Name','n_Repeated_Hits',
                                'Concentration_Level','Concentration')],
                 poslongresponse[,c('Intensity_Level','Intensity')])
#Filter out individual cases of NA on Mid or High Intensity, since it screws stuff up later with plotting:
poslong <- poslong[complete.cases(poslong$Intensity),]


#Plot data with 95% confidence interval, 95% prediction interval, and calibration interval given an intensity:
chemunique <- unique(poslong$Preferred_Name)
#Global Cal Curve Analysis (Low and Mid only):
pos_low_mid <- poslong[poslong$Concentration_Level != 'Corrected_HighConc',]

# ggsave(plot=slopehist,filename='Individual Curves_NoRepeatHitFilter_Histogram.png')
##save individual calibration curve data in an excel file:
# write.xlsx(pos_regress_df,file='PosMode_CalCurves_InverseCals.xlsx',
#            sheetName='Individual Linear Regressions',append=FALSE)
# write.xlsx(predcal, file='PosMode_CalCurves_InverseCals.xlsx',
#            sheetName='Prediction Interval Derived Calibrations (95%)', append=TRUE)
# write.xlsx(confcal, file='PosMode_CalCurves_InverseCals.xlsx',
#            sheetName='Confidence Interval Derived Calibrations (95%)', append=TRUE)
# ggplot(pos_regress_df_lowmid,aes(x=slope))+
#   geom_histogram(aes(y=..density..),
#                  color='black',fill='orange4')+
#   stat_function(fun=dnorm,args=list( mean = mean(pos_regress_df_lowmid$slope,na.rm=T),
#                                       sd=sd(pos_regress_df_lowmid$slope,na.rm=T)),
#                 lwd=1.5,color='red')+
#   labs(title=paste0('Slope distribution of 4th root Intensity_norm vs. 4th root Concentration\nµ = ',as.character(round(mean(pos_regress_df_lowmid$slope, na.rm=T),2)),
#                     ', ?? = ',as.character(round(sd(pos_regress_df_lowmid$slope, na.rm=T),2))))
# qqnorm(pos_regress_df_lowmid$slope)
# qqline(pos_regress_df_lowmid$slope,lwd=3,col='orange2')

pos_regress_df_lowmid <- data.frame(chemical = character(length(chemunique)),
                                    slope = double(length(chemunique)),
                                    intercept = double(length(chemunique)),
                                    RF = double(length(chemunique)))

predcal_lowmid <- data.frame(chemical = character(length(chemunique)),
                             concentration = double(length(chemunique)),
                             lbound = double(length(chemunique)),
                             est = double(length(chemunique)),
                             ubound = double(length(chemunique)),
                             y0=double(length(chemunique)),
                             glob_ubound = double(length(chemunique)),
                             err_msg = character(length(chemunique)))

pos_low_mid$Conc_4th<-pos_low_mid$Concentration^(1/4)
pos_low_mid$Intens_4th<-pos_low_mid$Intensity^(1/4)
calex_allregress<-ggplot(data=poslong,aes(x=log10(Concentration),y=log10(Intensity)))
#Redo slopes with just low/mid data:
for (i in 1:length(chemunique)){
  xydat2<-pos_low_mid[pos_low_mid$Preferred_Name == chemunique[i],
                      c('DTXSID','Preferred_Name','Concentration','Intensity', 'Conc_4th',
                        'Concentration_Level')]
  #comment out if statement line and lower bracket if you want regressions for everything, not just the
  #61 chems we get calibrates for:
  # if (length(xydat2$Intensity) > 2){ #Comment this line only and second-to-last bracket
    regress_lowmid<-lm(data=xydat2,
                       formula=log10(Intensity)~log10(Concentration))
    pos_regress_df_lowmid$chemical[i] <- chemunique[[i]]
    pos_regress_df_lowmid$slope[i] <- regress_lowmid$coefficients[[2]]
    pos_regress_df_lowmid$intercept[i] <- regress_lowmid$coefficients[[1]]
    
    if (#pos_regress_df_lowmid$slope[[i]] > 0.7 & pos_regress_df_lowmid$slope[[i]] < 1.3 &
        !is.na(pos_regress_df_lowmid$slope[[i]])){
      #Since medians don't always equate to intensities we have known concentrations for,
      #Force it to pick the closest known y0 to the median on the low end.
      predcal_lowmid$y0[[i]] <- max(xydat2[xydat2$Intensity <= 
                                             median(xydat2$Intensity),]$Intensity)
      predcal_lowmid$concentration[[i]] <- xydat2[xydat2$Intensity == 
                                                    predcal_lowmid$y0[[i]],]$Concentration
      errtest_p2 <- try(cal3<-calibrate(regress_lowmid,y0=log10(predcal_lowmid$y0[[i]]),
                                        interval='inversion',level=0.99),
                        silent=TRUE)
      if (class(errtest_p2) == 'try-error') {
        predcal_lowmid$chemical[i] <- chemunique[i]
        predcal_lowmid$lbound[i] <- NA
        predcal_lowmid$ubound[i] <- NA
        predcal_lowmid$est[i] <- NA
        predcal_lowmid$err_msg[i] <- errtest_p2[[1]]
      } else {
        predcal_lowmid$chemical[i] <-chemunique[i]
        predcal_lowmid$lbound[i] <- cal3$lower
        predcal_lowmid$est[i] <- cal3$estimate
        predcal_lowmid$ubound[i] <- cal3$upper
        predcal_lowmid$err_msg[i] <-NA
      }
      #perform 99% calibration interval off of regression 95% confidence interval and store in dataframe per chemical
      predcal_lowmid$err_msg[predcal_lowmid$ubound==Inf] <- "The calibration line is not well determined."
      predint_lowmid<-predict(regress_lowmid,interval='prediction',level=0.99)
      xydat2<-cbind(xydat2,predint_lowmid)
      calex<-ggplot(data=xydat2,aes(y=log10(Intensity), x=log10(Concentration)))+
        geom_point(size=3)+
        labs(title=xydat2$Preferred_Name,
             x='Log10(Concentration)',y='Log10(Intensity)')+
        geom_smooth(method='lm',formula = y~x,aes(color='royalblue'),se=FALSE,
                    lwd=1.5,alpha=0.5)
      calex2<-calex+
        geom_ribbon(data=xydat2,aes(ymin = lwr, ymax = upr, fill='99% Prediction'),alpha=0.25)+
        geom_hline(aes(yintercept=log10(predcal_lowmid$y0[[i]]),color='black'),lwd=1.5,lty=2)+
        geom_vline(aes(xintercept=cal3$estimate,color='purple'), lwd=1.5,lty=2)+
        geom_segment(aes(x=cal3$lower, y=log10(predcal_lowmid$y0[[i]]), xend=cal3$lower,
                         yend=-Inf,color='red'),lwd=1.5,lty=2)+
        geom_segment(aes(x=cal3$upper, y=log10(predcal_lowmid$y0[[i]]), xend=cal3$upper,
                         yend=-Inf,color='red'),lwd=1.5,lty=2)+
        scale_fill_manual("Y-Regression Intervals",
                          values = c('darkslategray','darkgreen'),
                          guide = guide_legend(override.aes= list(
                            linetype=c(0,0)
                          )))+
        scale_color_manual('Regression Statistics',
                           values=c('black',
                                    'purple',
                                    'red',
                                    'royalblue'),
                           labels=c('Y0 (New Observed Intensity)',
                                    'X0 (Calibration Estimate)',
                                    '99% Prediction Interval',
                                    'Regression Line'),
                           guide = guide_legend(override.aes = list(
                             fill = NA)))
      calex_allregress<-calex_allregress+
        geom_point(data=xydat2,size=3)+
        labs(x='Log10(Concentration)',y='Log10(Intensity)')+
        geom_smooth(data=xydat2,method='lm',color='red',formula = y~x,se=FALSE,
                    lwd=0.5,alpha=0.25)
      #save figures:
      # ggsave(plot=calex2,filename=paste0(unique(xydat2$DTXSID),'_lowmid.svg'),device='svg',
      #        height=4, width=6, units= 'in')
      # while (!is.null(dev.list()))   dev.off()
      # ggsave(plot=calex2,filename=paste0(unique(xydat2$DTXSID),'_lowmid.png'),device='png',
      #        dpi=300,height=4,width=6, units = 'in')
      # while (!is.null(dev.list()))   dev.off()
      # print(paste0('chemical ',i,' of ',length(chemunique), ' (',
      #              as.character(round(100*i/length(chemunique),2)),'% Complete) ',
      #              unique(xydat2$DTXSID),'_lowmid.svg',' saved successfully.'))
    }
  # }
}
qqnorm(pos_regress_df_lowmid[pos_regress_df_lowmid$slope < 1.36 &
                             pos_regress_df_lowmid$slope > 0.55,]$slope, 
       main=paste0('Individual Regression Slope\nNormal Q-Q Plot\n',),
       ylab='Slope Quantiles')
qqline(pos_regress_df_lowmid[pos_regress_df_lowmid$slope > 1.3 &
                               pos_regress_df_lowmid$slope < 0.7,],lwd=3,lty=2,col='red')
ggplot(pos_regress_df_lowmid)+
  geom_histogram(aes(x=slope),fill='red',color='black',alpha=0.5)

pos_regress_df_lowmid$RF <- 10^pos_regress_df_lowmid$intercept
ggplot(pos_regress_df_lowmid)+
  geom_histogram(aes(x=RF^(1/4),y=..density..),bins=18,
                 color='black',fill='orange2',alpha=0.5)+
  labs(x='Local Response Factors^(1/4)',
       y='Frequency')

qqnorm(pos_regress_df_lowmid$intercept,main='Log Regression Intercept\nNormal Q-Q Plot',
       ylab ='Intercept Quantiles')
qqline(pos_regress_df_lowmid$intercept,lwd=3,lty=2,col='green4')

predcal_lowmid$glob_est<-0
pos_regress_glob<-lm(data=pos_low_mid,formula=log10(Intensity) ~ log10(Concentration))
slope_glob<-pos_regress_glob$coefficients[2]
intercept_glob<-pos_regress_glob$coefficients[1]
#Try a test intensity with log(I) ~ 6 given issues with other chemicals
for (i in 1:length(predcal_lowmid$y0)){
  if (predcal_lowmid$y0[[i]] != 0){
    cal_pred_glob<-calibrate(pos_regress_glob,
                             y0=log10(predcal_lowmid$y0[[i]]),
                             interval='inversion',
                             level=0.99)
    predcal_lowmid$glob_ubound[[i]]<-cal_pred_glob$upper
    predcal_lowmid$glob_est[[i]]<-cal_pred_glob$estimate
  } else {
    cal_pred_glob<-calibrate(pos_regress_glob,
                             y0=log10(median(predcal_lowmid$y0)),
                             interval='inversion',
                             level=0.99)
    predcal_lowmid$glob_ubound[[i]]<-cal_pred_glob$upper
    predcal_lowmid$glob_est[[i]]<-cal_pred_glob$estimate
  }
}
predint_glob<-predict(pos_regress_glob,
                      newdata=pos_low_mid,
                      interval='prediction',
                      level=0.99)
pos_low_mid<-cbind(pos_low_mid,
                   predint_glob)
#run if you need to re-do prediction interval
#to remove prior interval columns
# pos_low_mid<-subset(pos_low_mid,select = -c(lwr,upr,fit)) 

# #Plot Global Cal. Residuals vs. 4th-root Conc.:
# pos_low_mid<-cbind(pos_low_mid,pos_regress_glob$residuals)
# names(pos_low_mid)[length(pos_low_mid)]<-'glob_residuals'
# ggplot(pos_low_mid,aes(x=Conc_4th,y=glob_residuals))+
#   geom_point(size=3)+
#   labs(x='mM Concentration^(1/4)',
#        y='Global Calibration Curve\nResiduals')

# Plot results of calibration interval analyses:
globalcal<-ggplot(data=pos_low_mid,
                  aes(x=Conc_4th,y=Intens_4th))+
  geom_point(size = 2)+
  geom_smooth(method='lm',
              formula= y~x,
              aes(color='Y Regression'),
              lwd=1.5,lty=1,
              se=FALSE)+
  geom_ribbon(aes(ymin=lwr,ymax=upr,fill='99% Prediction'),
              alpha=0.25)+
  geom_hline(aes(yintercept=median(predcal_lowmid$y0)^(1/4),
                 color='y0'),lwd=1.5,lty='dashed')+
                      geom_segment(aes(x=cal_pred_glob$estimate,y=median(predcal_lowmid$y0)^(1/4),
                   color='Concentration Estimate'),
               xend=cal_pred_glob$estimate, yend=-Inf,lwd=1.5,lty='dashed')+
  geom_segment(aes(x=cal_pred_glob$lower,y=median(predcal_lowmid$y0)^(1/4),
                   color='95% Calibration\nInterval (Pred.)'),
               xend=cal_pred_glob$lower,yend=-Inf,lwd=1.5,lty='dashed')+
  geom_segment(aes(x=cal_pred_glob$upper,y=median(predcal_lowmid$y0)^(1/4),
                   color='95% Calibration\nInterval (Pred.)'),
               xend=cal_pred_glob$upper,yend=-Inf,lwd=1.5,lty='dashed')+
  labs(title=paste0('Global Calibration Curve Results\nLow & Mid Concentrations Only\nIntensity^(1/4) = ',
                    as.character(round(pos_regress_glob$coefficients[[2]],4)),'*Concentration^(1/4) + ',
                    as.character(round(pos_regress_glob$coefficients[[1]],4))),
       x='mM Concentration^(1/4)',
       y='Intensity^(1/4)')+
  scale_color_manual('',aesthetics='color',
                     values=c('red','forestgreen','purple',
                              'blue','royalblue'),
                     labels=c('99% Prediction Interval', 'X0 (Estimated Concentration)',
                              'Regression Line','Y0 (New Observed Intensity)'),
                     guide=guide_legend(override.aes = list(
                       linetype = c(1,1,1,1))))+
  scale_fill_manual('',aesthetics = 'fill',
                    values=c('darkred'),
                    labels=c('Regression Prediction Interval'),
                    guide=guide_legend(override.aes = list(
                      linetype=c('blank'))))
setwd(figdir)

##Q-Q plot of global cal curve residuals
# qqnorm(pos_regress_glob$residuals,main='Global Calibration Curve (4th-root Space)\nQ-Q Plot of Residuals')
# qqline(pos_regress_glob$residuals,lwd=3,lty=2,col='purple')
# ggsave(plot=globalcal,filename='GlobalCalCurve.svg',device='svg',
#        height=4,width=6,units='in')
#Plot upper bound concentrations (95% Prediction Interval) 
#compared to the upper bound of the global cal curve
predcal_lowmid_filtered <- predcal_lowmid[predcal_lowmid$ubound != 0 &
                                            predcal_lowmid$ubound != Inf &
                                            complete.cases(predcal_lowmid$ubound),]

glob_df<- data.frame(chemical = character(length(chemunique)),
                     intensity = double(length(chemunique)),
                     concentration = double(length(chemunique)),
                     ubound = double(length(chemunique)))

for (i in 1:length(chemunique)){
  glob_df$chemical[i]<-chemunique[i]
  glob_df$intensity[i]<-log10(pos_low_mid[pos_low_mid$Preferred_Name == chemunique[i],]$Intensity[[1]])
  glob_df$concentration[i]<-log10(pos_low_mid[pos_low_mid$Preferred_Name == chemunique[i],]$Concentration[[1]])
  glob_cal<-calibrate(pos_regress_glob, y0=glob_df$intensity[[i]], interval='inversion', level=0.99)
  glob_df$ubound[i]<-glob_cal$upper
}
glob_df$ubound_individual<-predcal_lowmid$ubound
glob_df$EQ<-glob_df$ubound^4/glob_df$concentration^4

ggplot(predcal_lowmid_filtered)+
  geom_histogram(aes(x=glob_ubound^4/concentration,fill='ubound preds'),
                 binwidth=0.05,alpha=0.5,color='black')+
  labs(x='Predicted Conc. / Known Conc.',
       y='Frequency',
       title='Global Calibration Curve Error Quotients\nUpper Bound Concentration Estimates')+
  scale_x_log10()+
  geom_vline(aes(xintercept=1),color='black',lwd=1.5,lty=2)
#Filter with trimmed slope range 0.7 < x < 1.3 and combine pos_regress_df_lowmid and predcal_lowmid_filtered columns:
pos_regress_df_lowmid$known_conc <- NA
pos_regress_df_lowmid[pos_regress_df_lowmid$chemical %in% predcal_lowmid_filtered$chemical,]$known_conc<-
  predcal_lowmid_filtered[predcal_lowmid_filtered$chemical %in% pos_regress_df_lowmid$chemical,]$concentration
pos_regress_df_lowmid$glob_ubound <- NA
pos_regress_df_lowmid[pos_regress_df_lowmid$chemical %in% predcal_lowmid_filtered$chemical,]$glob_ubound<-
  predcal_lowmid_filtered[predcal_lowmid_filtered$chemical %in% pos_regress_df_lowmid$chemical,]$glob_ubound
pos_regress_df_lowmid$glob_est <- NA
pos_regress_df_lowmid[pos_regress_df_lowmid$chemical %in% predcal_lowmid_filtered$chemical,]$glob_est<-
  predcal_lowmid_filtered[predcal_lowmid_filtered$chemical %in% pos_regress_df_lowmid$chemical,]$glob_est
pos_regress_df_lowmid$indiv_ubound <- NA
pos_regress_df_lowmid[pos_regress_df_lowmid$chemical %in% predcal_lowmid_filtered$chemical,]$indiv_ubound<-
  predcal_lowmid_filtered[predcal_lowmid_filtered$chemical %in% pos_regress_df_lowmid$chemical,]$ubound
pos_regress_df_lowmid$indiv_est <- NA
pos_regress_df_lowmid[pos_regress_df_lowmid$chemical %in% predcal_lowmid_filtered$chemical,]$indiv_est<-
  predcal_lowmid_filtered[predcal_lowmid_filtered$chemical %in% pos_regress_df_lowmid$chemical,]$est


EQhist_glob<-ggplot(pos_regress_df_lowmid,
                    aes(x=glob_ubound^4/known_conc))+
  geom_histogram(aes(fill='darkred'),
                 alpha=0.5,color='black',na.rm=TRUE)+
  labs(x='log10(Global Upper Bound Conc. / Known Conc.)',
       y='Frequency')+
  geom_vline(aes(xintercept=1),lty=2,lwd=1.5)+
  scale_x_log10()+
  scale_fill_manual(values='darkred',labels='Global EQ')
pos_regres
uboundhist_glob_indiv<-ggplot(pos_regress_df_lowmid,aes(x=glob_ubound^4/10^indiv_ubound))+
  geom_histogram(aes(fill='Inverse Regression\nUpper Bounds'),bins = 18,
                 alpha=0.5,color='black',na.rm=TRUE)+
  labs(x='Global Upper Bound^4 / 10^Local Upper Bound',
       y='Frequency')+
  geom_vline(aes(xintercept=1),lty=2,lwd=1.5)+
  scale_x_log10()+
  scale_fill_manual(values='maroon')
# ggsave(plot=uboundhist,filename='ENTACTStandardMixtures_UBoundConcEstimates_GlobVSIndividual.png',dpi=300)


uboundhist_glob_est<-ggplot(pos_regress_df_lowmid,
                            aes(x=glob_ubound^4/indiv_est^4))+
  geom_histogram(bins=18,
                 alpha=0.5,
                 fill='forestgreen',
                 color='black',
                 na.rm=TRUE)+
  labs(x='Global Upper Bound / Local Point Estimate',
       y='Frequency',
       title='Global Calibration Curve')+
  scale_x_log10()
pos_regress_df_lowmid$glob_vs_indiv_ubound<-0
pos_regress_df_lowmid$glob_vs_indiv_ubound<-
  pos_regress_df_lowmid$glob_ubound^4/ 10^pos_regress_df_lowmid$indiv_ubound

uboundhist_indiv_est<-ggplot(pos_regress_df_lowmid,
                             aes(x=10^indiv_ubound/10^indiv_est))+
  geom_histogram(bins = 18,
                 alpha=0.7,
                 fill='darkgray',
                 color='black',
                 na.rm=TRUE)+
  labs(x='Local Upper Bound / Local Point Estimate',
       y='Frequency',
       title='Individual Calibration Curves')+
  scale_x_log10()
median(10^pos_regress_df_lowmid$indiv_ubound/10^pos_regress_df_lowmid$indiv_est,na.rm=T)
#histogram RFs = 10^(intercept)
ggplot(pos_regress_df_lowmid,
       aes(x=RF^(1/4)))+
  geom_histogram(fill='blue',color='black',bins=18)+
  labs(x='Response Factor^(1/4)',
       y='Frequency',
       title='Individual Calibration Curves\n RF = [10^(Indiv. Intercept)]^1/4 Distribution')

qqnorm(pos_regress_df_lowmid[pos_regress_df_lowmid$chemical %in% predcal_lowmid_filtered$chemical,]$RF^(1/4))
qqline(pos_regress_df_lowmid[pos_regress_df_lowmid$chemical %in% predcal_lowmid_filtered$chemical,]$RF^(1/4),
       lwd=3,col='blue')
qqnorm(pos_regress_df_lowmid$RF^(1/4),main='All Response Factors\nNormal Q-Q Plot')
qqline(pos_regress_df_lowmid$RF^(1/4),lwd=3,col='blue')

#Show bounding for RFs and for global cal for single chemical, ethionamide:
Ethionamide_x0 <- predcal_lowmid_filtered[predcal_lowmid_filtered$chemical == 'Ethionamide',]$concentration
Ethionamide_y0 <- predcal_lowmid_filtered[predcal_lowmid_filtered$chemical == 'Ethionamide',]$y0
pos_regress_df_lowmid$RFconc <- (Ethionamide_y0^(1/4)/pos_regress_df_lowmid$RF^(1/4))^4
Ethionamide_glob<-calibrate(pos_regress_glob,
                            interval='inversion',
                            y0=Ethionamide_y0^(1/4),
                            level=0.99)
#for stat intervals:
t_crit_99p5_n523<-qt(0.995,df=length(pos_regress_df_lowmid$RF)-1)
RF_mean_4throot <- mean(pos_regress_df_lowmid$RF^(1/4))
RF_sigma_4throot <- sd(pos_regress_df_lowmid$RF^(1/4))
RF_predint_99p5 <- RF_mean_4throot+t_crit_99p5_n523*
  RF_sigma_4throot*sqrt(1+(1/length(pos_regress_df_lowmid$RF)))
RF_predint_0p5 <- RF_mean_4throot-t_crit_99p5_n523*
  RF_sigma_4throot*sqrt(1+(1/length(pos_regress_df_lowmid$RF)))
eth_df<- data.frame(bounding = character(2),
                    upper = double(2),
                    lower = double(2))
eth_df$bounding <- c('Global Calibration','Response Factor')
eth_df$upper <- c(Ethionamide_glob$upper^4,(Ethionamide_y0^(1/4)/RF_predint_0p5)^4)
eth_df$lower <- c(Ethionamide_glob$lower^4,(Ethionamide_y0^(1/4)/RF_predint_99p5)^4)
#plot upper and lower bounds for each 
ggplot(eth_df,aes(x=bounding))+
  geom_point(aes(y=upper,color='Upper Bound'),size=5)+
  geom_point(aes(y=lower,color='Lower Bound'),size=5)+
  labs(x='Bounding Method',y='Predicted Concentration\nBounds (mM)',title='Ethionamide')+
  geom_hline(aes(yintercept=Ethionamide_x0),color='forestgreen', alpha=0.5, lwd=1.5,lty=2)+
  scale_color_manual('',values=c('royalblue','darkred'),labels=c('Lower','Upper'))+
  scale_y_log10()
subset(pos_regress_df_lowmid, select = -c(RFconc, RFupper))
pos_regress_df_lowmid$y0 <-0
pos_regress_df_lowmid$y0 <- predcal_lowmid$y0
pos_regress_df_lowmid$RF_conc_upper <- (pos_regress_df_lowmid$y0^(1/4)/RF_predint_0p5)^4
pos_regress_df_lowmid$RF_conc_lower <- (pos_regress_df_lowmid$y0^(1/4)/RF_predint_99p5)^4
pos_regress_df_lowmid$RF_upper_glob_upper <- pos_regress_df_lowmid$RF_conc_upper/pos_regress_df_lowmid$glob_ubound^4
pos_regress_df_lowmid$RF_upper_indiv_upper <- pos_regress_df_lowmid$RF_conc_upper/10^pos_regress_df_lowmid$indiv_ubound
pos_regress_df_lowmid$RF_upper_indiv_est <- pos_regress_df_lowmid$RF_conc_upper/10^pos_regress_df_lowmid$indiv_est

ggplot(pos_regress_df_lowmid)+
  geom_histogram(aes(x=RF_upper_glob_upper,fill='RF Upper / Global Upper'),bins=18,na.rm=T,alpha=0.5,color='black')+
  scale_x_log10(limits=c(1,2000))+
  labs(x='Concentration Prediction Ratios',
       y='Frequency')
ggplot(pos_regress_df_lowmid)+
  geom_histogram(aes(x=RF_upper_indiv_upper,fill='RF Upper / Individual Upper'),bins=18,na.rm=T,alpha=0.5,color='black')+
  scale_x_log10(limits=c(1,2000))+
  labs(x='Concentration Prediction Ratios',
       y='Frequency')
  ggplot(pos_regress_df_lowmid)+  
  geom_histogram(aes(x=RF_upper_indiv_est,fill='RF Upper / Individual Estimate'),bins=18,na.rm=T,alpha=0.5,color='black')+
  scale_x_log10(limits=c(1,2000))+
  labs(x='Concentration Prediction Ratios',
       y='Frequency')

ggplot(pos_regress_df_lowmid)+
  geom_histogram(aes(x=(Ethionamide_y0^(1/4)/RF^(1/4))^4),
                 bins=18,fill='purple',color='black')+
  geom_vline(aes(xintercept=Ethionamide_x0,color='2'),lwd=1.5,lty=2)+
  geom_vline(aes(color='1', xintercept=
                   pos_regress_df_lowmid[pos_regress_df_lowmid$chemical == 'Ethionamide',]$glob_ubound^4),
             lwd=1.5,lty=2)+
  geom_vline(aes(color='3', xintercept=
                   pos_regress_df_lowmid[pos_regress_df_lowmid$chemical == 'Ethionamide',]$RF_conc_upper),
             lwd=1.5,lty=2)+
  scale_x_log10()+
  labs(title='Ethionamide',
       x='RF-Predicted mM Concentrations',
       y='Frequency')+
  scale_color_manual('',values=c('red','black','forestgreen'),
                     labels=c('Global Upper Bound Estimate\nfrom 99% Prediction Interval',
                              'Known Concentration',
                              'RF Upper Bound Estimate\nfrom 99% Prediction Interval'))
pos_regress_df_lowmid$y0<-NA
pos_regress_df_lowmid[pos_regress_df_lowmid$chemical %in% predcal_lowmid$chemical,]$y0<-
  predcal_lowmid[predcal_lowmid$chemical %in% pos_regress_df_lowmid$chemical,]$y0
pos_regress_df_lowmid$glob_ubound<-NA
pos_regress_df_lowmid[pos_regress_df_lowmid$chemical %in% predcal_lowmid$chemical,]$glob_ubound<-
  predcal_lowmid[predcal_lowmid$chemical %in% pos_regress_df_lowmid$chemical,]$glob_ubound

intercept_lower<- mean(pos_regress_df_lowmid$intercept)-qt(0.995,df=length(pos_regress_df_lowmid$intercept)-1)*
  sd(pos_regress_df_lowmid$intercept)*sqrt(1+(1/length(pos_regress_df_lowmid$intercept)))
pos_regress_df_lowmid$RFconc<-NA
pos_regress_df_lowmid$RFconc_upper<-log10(pos_regress_df_lowmid$y0)/intercept_lower
pos_regress_df_lowmid$RF_glob_ratio<-NA
pos_regress_df_lowmid$RF_glob_ratio<-10^pos_regress_df_lowmid$RFconc_upper/10^pos_regress_df_lowmid$glob_ubound

ggplot(pos_regress_df_lowmid)+
  geom_histogram(aes(x=RF_glob_ratio),fill='red',color='black',alpha=0.5,bins=28)+
  geom_vline(aes(xintercept=1),lwd=1.5,lty=2,color='black')+
  labs(x='Upper Bound Concentration Ratios',
       y='Frequency',
       title='ESI-, log(RF) upper bound vs. global cal. upper bound')
