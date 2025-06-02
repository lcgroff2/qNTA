#Inverse calibration analysis for Brita GC-HRMS (Louis Groff 12 Feb. 2021)
library(ggplot2)
library(data.table)
library(investr)
library(xlsx)
library(readxl)
library(bestNormalize)
library(lme4)

#Change this to the directory that holds the Brita Data:
gcdatadir <- 'C:/Users/lgroff/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/Data/LiberatoreH GC-MS Data/Brita'
setwd(gcdatadir) #Change to your working directory if R isn't already pointing at it.
gcdata <- read_xlsx('BritaGCDataFrame.xlsx') # load in data (from readxl package)
#Melt TIC data from wide column format to long column format.
TIC_melt<-melt(setDT(gcdata[,c('Component Name','RT',
                               "TIC1 0.4 ppm","TIC2 0.4 ppm","TIC3 0.4 ppm",
                               "TIC1 2 ppm","TIC2 2 ppm","TIC3 2 ppm",
                               "TIC1 10 ppm","TIC2 10 ppm","TIC3 10 ppm")]),
               id.vars=c('Component Name','RT'),
               measure.vars=c("TIC1 0.4 ppm","TIC2 0.4 ppm","TIC3 0.4 ppm",
                              "TIC1 2 ppm","TIC2 2 ppm","TIC3 2 ppm",
                              "TIC1 10 ppm","TIC2 10 ppm","TIC3 10 ppm"),
               variable.name = 'TIC_level',
               value.name = 'TIC')
#Melt Conc data from wide column format to long column format.
conc_melt<-melt(setDT(gcdata[,c("Component Name",
                                "uM low","uM mid","uM high")]),
                id.vars='Component Name',
                measure.vars=c("uM low","uM low","uM low",
                               "uM mid","uM mid","uM mid",
                               "uM high","uM high","uM high"),
                variable.name='conc_level',
                value.name='uM_conc')
#Concatenate TIC and Conc on columns
gc_melt<-cbind(conc_melt[,c('Component Name','conc_level','uM_conc')],TIC_melt[,c('RT','TIC')])
gc_melt$RF <- gc_melt$TIC/gc_melt$uM_conc #Calculate Response Factors as TIC/Conc
# GC_boxmethod<-boxcox(gc_melt$RF)
# gc_melt$RF_boxcox<-GC_boxmethod$x.t
# GC_tRF_1st<- mean(gc_melt$RF_boxcox,na.rm=T)+qt(0.005,df=length(gc_melt$RF_boxcox)-1)*sd(gc_melt$RF_boxcox,na.rm=T)*sqrt(1+1/length(gc_melt$RF_boxcox))
RFplot<-ggplot()+
  geom_histogram(aes(x=gc_melt$RF,
                     y=..density..,
                     fill='I/C RF'),alpha=0.5,color='black')+
  # geom_histogram(aes(x=boxcox(gcdata$RFintercept)$x.t,
  #                    y=..density..,
  #                    fill='Linear\nIntercept RF'),
  #                color='black',alpha=0.5)+
  # geom_vline(aes(xintercept=GC_tRF_1st),lwd=1.5,lty=2)+
  scale_x_log10()+
  labs(x=' RF')+
  theme_classic()
#how to save a hi-res image (jpeg in this case, but can do png, tiff, eps, etc.)
ggsave(plot=RFplot,filename='RFIntercept_RFIC_BritaGC.jpg',dpi=1200,device='jpeg')

#list of unique chemical names:
chemunique<-unique(gc_melt$`Component Name`)

#used for randomly selecting y0 values for the inverse predictions
#setting the seed value on the random number generator will allow you to
#reproduce the same pseudorandom dataset, the number choice doesn't matter, 
#it's just picking the 16384th random state. You can set it to any number.
set.seed(16384) 

#generate a random sample of numbers valued between 1-100,
#in a list the size of the number of rows in gc_melt
#divide those integer numbers by 100 to produce a random number between 0 and 1.
rand_idx <- sample(1:100,nrow(gc_melt),size=nrow(gc_melt))/100

#Preallocate empty columns of dataframe to insert results of inverse prediction
gcdata$y0<-NA
gcdata$Conc_CC <- NA
gcdata$Conc_CC_Upper <- NA
gcdata$slope <- NA #linear model only
gcdata$intercept <- NA #linear model only
# gcdata$x2_coeff <- NA #x^2 model only
# gcdata$x_coeff <- NA #x^2 model only
# gcdata$x2_intercept <- NA #x^2 model only

#preallocate an empty list equal to the length of the list of unique chemicals
#to store the successful calibration curve ggplot visualizations in
gglist <- vector(mode='list',length=length(chemunique))

#iterate through each unique chemical's linear cal curve and generate 99% Prediction Interval
for (i in 1:length(chemunique)){
  tmpdf <- gc_melt[gc_melt$`Component Name`==chemunique[[i]],] #subset the data for an individual chemical
  #pick a random observation for y0:
  #choose intensity row index on the subset dataframe based on ith previously determined
  #random value between 0-1, multiplied by the number of rows in the dataframe, rounded up
  #e.g., if rand_idx[[i]] = 0.48 and nrow(tmpdf) = 5, then 0.48*5 = 2.4
  #ceiling() always rounds up if the number is greater than the current integer
  #so ceiling(2.4) = 3, so we pick the 3rd TIC value from the subset dataframe.
  gcdata$y0[[i]]<-tmpdf[ceiling(rand_idx[[i]]*nrow(tmpdf))]$TIC 
# #  linear model:
  #store linear model as an object whose parameters can be called upon later by other functions
  tmplm <- lm(tmpdf,formula=log10(TIC)~log10(uM_conc))
  #store slope and intercept
  gcdata$intercept[[i]] <- tmplm$coefficients[[1]]
  gcdata$slope[[i]] <- tmplm$coefficients[[2]]
  
  #Generate 99% Prediction Interval using predict() function
  tmp_pi<-predict(tmplm, #linear model object
                  level=0.99, #confidence level
                  interval='prediction') #type of interval to calculate
  
#  linear calibrations:
  #generate inverse prediction based on the prediction interval above as temporary variable
  #encapsulated in a try statement since sometimes calibrate()/invest() will fail
  #but you don't want it to kill the remaining for loop execution. 
  #Use invest() for nonlinear inverse calibrations. silent = T so as to not 
  #clutter the console with warnings/error messages if try statement fails.
  tmpcal<-try(calibrate(tmplm, #linear model object
                        y0=log10(gcdata$y0[[i]]), #randomly selected y0
                        interval='inversion', #calculating an inverse calibration interval
                        level=0.99), # confidence level
              silent=T)
  
  #uncomment for x^2 unweighted polynomial model, and uncomment or update 
  #associated calls that use tmplm with tmplm2:
  # tmplm2 <- lm(data=tmpdf,formula=log10(TIC)~log10(uM_conc)+I(log10(uM_conc)^2))
  # #Non-linear calculations for x^2 fits
  
  # tmpcal <- try(calibrate(tmplm2, #linear model object
  #                         y0=log10(gcdata$y0[[i]]), #y-value to predict x-range from
  #                         interval='inversion', #calculating an inverse prediction interval
  #                         extendInt = 'yes', #if y0 is outside of the original PI range, extend
  #                         level=0.99, #confidence level
  #                         maxiter = 1e4, #max number of iterations to attempt for convergence
  #                         tol=0.001), #minimum difference between two fit attempts to say "good enough" on convergence
  #               silent=T)
  
  #store cal curve in temporary variable tmpgg
  tmpgg <- ggplot(tmpdf,
                  aes(x=log10(uM_conc),
                      y=log10(TIC)))+
    geom_point()+
    #use geom_ribbon to visualize our calculated prediction interval
    geom_ribbon(aes(ymin=tmp_pi[,'lwr'],
                    ymax=tmp_pi[,'upr']),
                alpha=0.25)+
    #standard error is different, so we set se=F to avoid confusion with stacking
    #on top of the prediction interval from the lines above.
    geom_smooth(method='lm',
                formula=y~x,
                se=F)+
    #axis labels and title for plot
    labs(x='Log10 uM Concentration',
         y='Log10 TIC',
         title=unique(tmpdf$`Component Name`)[[1]])+
    theme_classic()
  
  #If the inverse calibration fails, its data type/class will be 'try-error'
  #Check the data type/class of tmpcal, if it isn't a try-error, that means the
  #inverse calibration succeeded. We should store those results
  if(class(tmpcal)!='try-error'){
    #extra sanity check because with x^2 fits sometimes the upper and lower bounds flip-flop
    if((tmpcal$upper/tmpcal$estimate) > 1){
      # gcdata$x2_coeff[[i]] <- tmplm2$coefficients[[3]]
      # gcdata$x_coeff[[i]] <- tmplm2$coefficients[[2]]
      # gcdata$x2_intercept[[i]] <- tmplm2$coefficients[[1]]
      gcdata$slope[[i]] <- tmplm$coefficients[[2]]
      gcdata$intercept[[i]] <- tmplm$coefficients[[1]]
      gcdata$Conc_CC[[i]] <- tmpcal$estimate
      gcdata$Conc_CC_Upper[[i]] <- tmpcal$upper
      
      #build the log-log inverse calibration visualization for the ith chemical
      tmpgg <- tmpgg+
        geom_hline(yintercept=log10(gcdata$y0[[i]]),
                  lwd=1,
                  lty=2,
                  color='red')+
        geom_vline(xintercept=log10(tmpcal$lower),
                  lwd=1,
                  lty=2,
                  color='forestgreen')+
        geom_vline(xintercept=log10(tmpcal$upper),
                  lwd=1,
                  lty=2,
                  color='forestgreen')
    }
    else{
      gcdata$Conc_CC_Upper[[i]]<-Inf
    }
  }
  gglist[[i]]<-tmpgg #store each successful run in gglist.
}
write.csv(gcdata, file='gcdata_x2fits.csv')

gc_slope_hist<-ggplot(gcdata,aes(x=Conc_CC_Upper/Conc_CC))+
  geom_histogram(fill='slateblue',color='black',alpha=0.5)+
  theme_classic()+
  labs(x='Regression Slope',
       y='Frequency',
       title = paste0('mean slope = ',
                      round(mean(gcdata$slope,na.rm=T),4),
                      '\nsd = ',
                      round(sd(gcdata$slope,na.rm=T),4)))

#Generate Response Factors from exponentiation of log-log 
#regression intercepts:
gcdata$RF<-10^gcdata$intercept
#Use bestNormalize to determine best transformation 
#small dataset, use Leave-One-Out CV more computationally efficient and robust:
set.seed(16384)
bestRT <- bestNormalize(gcdata$RT,new_transforms = custom_transform, standardize=F, loo = T, allow_orderNorm = F)
#store transformed RT in gcdata:
gcdata$bestN_RT<-bestRT$x.t

#check normality of transform result:
dev.new(width=4,height=4)
qqnorm(bestRT$x.t,main='Normal Q-Q Plot\nBest RT Transform (Yeo-Johnson)')
qqline(bestRT$x.t,lwd=3,col='red')

#Repeat to find best transform for RF:
set.seed(16384)
bestRF <- bestNormalize(gcdata$RF,new_transforms = custom_transform, standardize=F, loo = T, allow_orderNorm = F)
#store transformed RF in gcdata:
gcdata$bestN_RF<-bestRF$x.t

#Check normality of transform result:
dev.new(width=4,height=4)
qqnorm(bestRF$x.t,main='Normal Q-Q Plot\nBest RF Transform (RF^(1/4))')
qqline(bestRF$x.t,lwd=3,col='red')

alpha = 0.01
level = 1-alpha/2 #two-sided T-distribution sample
#Calculate 1st percentile RF and generate upper bound concentrations:
tRF_lower <- mean(gcdata$bestN_RF,na.rm=T)-
  qt(level,df=length(gcdata$bestN_RF)-1)*sd(gcdata$bestN_RF,na.rm=T)*sqrt(1+1/length(gcdata$bestN_RF))
gcdata$Conc_RF_Upper <- gcdata$y0/tRF_lower^4

#Generate linear model and PI from transformed RF and transformed RT:
brita_lm <- lm(gcdata,formula=bestN_RT~bestN_RF)
brita_PI <- predict(brita_lm,
                    interval='prediction',
                    level=0.99)
#Plot transformed RT vs. transformed RF:
gcdata_plot<-ggplot(gcdata,
       aes(y=bestN_RT, x=bestN_RF))+
  geom_point()+
  labs(x=expression(Response~Factor^(1/4)),
       y=expression(((RT+1)^lambda-1)/lambda),
       title='Brita GC-EI, 99% PI')+
  geom_smooth(method='lm',
              formula=y~x,
              se=F)+
  geom_ribbon(aes(ymin=brita_PI[,'lwr'],ymax=brita_PI[,'upr']),
              alpha=0.25)+
  theme_classic()

#Check linear model residual normality (Normal Q-Q plot):
dev.new(width=4,height=4)
qqnorm(brita_lm$residuals,main='YeoJohson(RT) vs. RF^(1/4) Residuals\nNormal Q-Q Plot')
qqline(brita_lm$residuals,lwd=3,col='red')

#Check transformed linear model residual heteroskedasticity:
gcdata$residuals <- brita_lm$residuals
bestN_hetero_RT<-ggplot(gcdata,aes(x=bestN_RT,y=residuals))+
  geom_point()+
  labs(x=expression(((RT+1)^lambda-1)/lambda),
       y='Residuals')+
  theme_classic()
#repeat for RFs
bestN_hetero_RF<-ggplot(gcdata,aes(x=bestN_RF,y=residuals))+
  geom_point()+
  labs(x=expression(Response~Factor^(1/4)),
       y='Residuals')+
  theme_classic()

#Perform inverse calibration on Yeo-Johnson RT vs. RF^(1/4) fit
#to generate lower bound RF^(1/4) from regression 99% PI:
gcdata$RF_4throot_lower<-NA
for(i in 1:length(chemunique)){
  tmpcal<-calibrate(brita_lm,
                    y0=gcdata$bestN_RT[[i]],
                    interval='inversion',
                    level=0.99)
  gcdata$RF_4throot_lower[[i]]<-tmpcal$lower
  #impute 1st percentile tRF from 1D Distribution when inverse 
  #calibration RF is < 1st percentile RF from 1D distribution:
  if (gcdata$RF_4throot_lower[[i]] < tRF_lower){
    gcdata$RF_4throot_lower[[i]] <- tRF_lower
  }
}
#generate upper bound concentration predictions from the tRT vs tRF method:
gcdata$Conc_RT_Upper<- gcdata$y0/gcdata$RF_4throot_lower^4

#Generate error quotients for each method:
gcdata$EQ_CCUpper_CCEst <- gcdata$Conc_CC_Upper/gcdata$Conc_CC
gcdata$EQ_RFUpper_CCEst <- gcdata$Conc_RF_Upper/gcdata$Conc_CC
gcdata$EQ_RTUpper_CCEst <- gcdata$Conc_RT_Upper/gcdata$Conc_CC
gcdata$EQ_RFUpper_RTUpper <- gcdata$Conc_RF_Upper/gcdata$Conc_RT_Upper

err_df<- melt(setDT(gcdata[,c("Component Name",
                              "EQ_CCUpper_CCEst",
                              "EQ_RFUpper_CCEst",
                              "EQ_RTUpper_CCEst")]),
              id.vars = 'Component Name',
              measure.vars = c("EQ_CCUpper_CCEst",
                               "EQ_RFUpper_CCEst",
                               "EQ_RTUpper_CCEst"),
              variable.name = 'Method',
              value.name = 'Error')

improved <- 0
for(i in 1:length(gcdata$EQ_RFUpper_RTUpper)){
  if(gcdata$EQ_RTUpper_CCEst[[i]] < gcdata$EQ_RFUpper_CCEst[[i]]){
    improved <- improved+1
  }
}

gc_err_violins <- ggplot(err_df,
       aes(x=Error,y=Method))+
  geom_boxplot(color='blue',
               width=0.1,
               lwd=1)+
  geom_violin(fill=NA,
              scale='width',
              lwd=1)+
  scale_x_log10()+
  theme_classic()+
  labs(x='Error Quotient',
       y='Method Comparison')+
  scale_y_discrete(labels=c('CC Upper / CC Est','RF Upper / CC Est','RT Upper / CC Est'))
# ggsave(gcdata,filename='BritaGC_ErrorViolinPlots.jpeg',device='jpeg',dpi=1200)
# write.xlsx(gcdata,file='BritaGCSummaryTable.xlsx')
# ggsave(gcdata_plot,filename='BritaGC_Y-J_RTvsRF4th.jpeg',device='jpeg',dpi=1200)

brita_percentiles <- function(data, descriptors, runs, samples)
{
  #store rng seeds for n runs for reproducibility
  #comment out for true random MC
  set.seed(16384)
  
  rf_data <- data[,..descriptors]
  #remove NAs from RF list and preallocate data frames/columns:
  rf_data <- rf_data[complete.cases(rf_data$RF),]
  percentile_dist <- data.frame(percentile_2p5th = double(runs),
                                percentile_97p5th = double(runs))
  percentile_dist$RF_sample <- c()
  percentile_dist$Component_Name <- c()
  unique_chems <- unique(rf_data$Component_Name)
  print(paste0('Begin boostrap Simulation for N = ',runs, ' runs. '))
  for (i in 1:runs)
  {
    #get indexes of which chemicals to randomly sample:
    chem_idx_boot <- sample(length(unique_chems),
                            size=samples,
                            replace=T)
    #resampled list of chemicals and their indices:
    chem_resample <- unique_chems[chem_idx_boot]
    rf_resample <- data.frame(RF = double(samples),
                              Component_Name = chem_resample)
    #choosing one random RF method:
    for(j in 1:length(chem_resample)){
      sample_idx <- sample(nrow(rf_data[rf_data$Component_Name == chem_resample[[j]],'RF']),size=1)
      rf_resample[j,'RF'] <- rf_data[rf_data$Component_Name == chem_resample[[j]],]$RF[[sample_idx]]
    }
    ##calculating geometric mean of RFs method:
    # for(j in 1:length(chem_resample))
    # {
    #   rf_resample[j,'RF'] <- 10^mean(log10(rf_data[rf_data$Component_Name == chem_resample[[j]],]$RF))
    # }
    #calculate 1st percentile RF and 5th percentile RF (one-sided):
    percentile_dist$percentile_2p5th[[i]] <- quantile(rf_resample$RF,0.025)[[1]]
    percentile_dist$percentile_97p5th[[i]] <- quantile(rf_resample$RF,0.975)[[1]]
    percentile_dist$RF_sample[[i]] <- rf_resample$RF
    percentile_dist$Component_Name[[i]] <- rf_resample$Component_Name
    print(paste0(round(100*i/runs,2),'% complete.'))
  }
  return(percentile_dist)
}

percentile_bootstrap_kfoldCV <- function(data, descriptors, folds, runs, samples)
{
  set.seed(16384)
  rf_data <- data[,..descriptors]
  
  #remove NAs from RF list and preallocate data frames/columns:
  rf_data <- rf_data[complete.cases(rf_data$RF),]
  chemlist_CV <- unique(rf_data$Component_Name)
  chem_idx_CV <- sample(length(unique(rf_data$Component_Name)),
                        size=samples,
                        replace=F)
  fold_size <- length(chemlist_CV)/folds
  mc_out <- vector(mode='list',length=folds)
  train_folds <- mc_out
  test_folds <- mc_out
  for(i in 1:folds){
    print(paste0('CV Run #',i))
    rf_test <- rf_data[rf_data$Component_Name %in% chemlist_CV[chem_idx_CV[(1+(i-1)*fold_size):(i*fold_size)]],]
    rf_train <- rf_data[!(rf_data$Component_Name %in% chemlist_CV[chem_idx_CV[(1+(i-1)*fold_size):(i*fold_size)]]),]
    train_percentiles<-brita_percentiles(data=rf_train,
                                         descriptors,
                                         runs,
                                         samples)
    train_folds[[i]] <- rf_train
    test_folds[[i]] <- rf_test
    mc_out[[i]]<-train_percentiles
  }
  cv_out <- list(train_folds,test_folds,mc_out)
  return(cv_out)
}

mySumm <- function(.) {
  predict(., newdata=gc_melt, re.form=NULL)
}
####Collapse bootstrap into median, 95% PI
sumBoot <- function(merBoot) {
  return(
    data.frame(fit = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))),
               lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE))),
               upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE)))
    )
  )
}


RF_transform<-bestNormalize(gc_melt$RF,
                            standardize=F,
                            allow_orderNorm = F,
                            norm_stat_fn = function(x) 1-shapiro.test(x)$statistic)
boxlambda<-RF_transform$chosen_transform$lambda
# colnames(gc_melt)[[1]]<-'Component_Name'
Brita_Default_Bootstrap<-percentile_bootstrap_kfoldCV(data = gc_melt,
                                                      descriptors = c('RF','Component_Name'),
                                                      samples = length(unique(gc_melt$Component_Name)),
                                                      runs = 10000,
                                                      folds = 5)
Brita_DefaultRF_gg<-ggplot()+
  geom_histogram(aes(x=RF),
                 fill='blue',
                 color='black')+
  geom_vline(aes(xintercept=quantile(Brita_Default_Bootstrap$percentile_2p5th,0.5)),
             lwd=1.5,
             lty=2)+
  geom_vline(aes(xintercept=quantile(Brita_Default_Bootstrap$percentile_97p5th,0.5)),
             lwd=1.5,
             lty=2)+
  labs(x='Response Factor',
       y='Frequency')+
  theme_classic()

saveRDS(Brita_Default_Bootstrap,'Brita_Default_Bootstrap.rds')
