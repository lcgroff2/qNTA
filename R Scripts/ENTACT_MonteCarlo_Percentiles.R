#ENTACT_MonteCarlo_Percentiles.R - Louis Groff
#May 14 2021

library(ggplot2)
library(readxl)
library(xlsx)
library(gridExtra)
library(bestNormalize)
library(lme4)
library(investr)

bn_obj_shapiro_pos <- bestNormalize(x=data$RF,
                                    allow_orderNorm = F,
                                    loo=T,
                                    norm_stat_fn = function(x) 1-shapiro.test(x)$statistic)
bn_obj_shapiro_neg <- bestNormalize(x=data_neg$RF,
                                    allow_orderNorm = F,
                                    loo=T,
                                    norm_stat_fn = function(x) 1-shapiro.test(x)$statistic)

bn_obj_marg_pos <- bestNormalize(x=data$RF,
                                 allow_orderNorm = F,
                                 loo=T,
                                 norm_stat_fn = function(x) sqrt(mean(residuals(lm(formula=x~data$logIE_Pred))^2)))
bn_obj_marg_neg <- bestNormalize(x=data_neg$RF,
                                 allow_orderNorm = F,
                                 loo=T,
                                 norm_stat_fn = function(x) sqrt(mean(residuals(lm(formula=x~data_neg$logIE_Pred))^2)))

#lm objects to calculate marginal residuals
lm_box_pos = lm(formula = bn_obj_pos$other_transforms$boxcox$x.t~data$logIE_Pred)
lm_box_neg = lm(formula = bn_obj_neg$other_transforms$boxcox$x.t~data_neg$logIE_Pred)

#Marginal Residual QQ Plots ESI+:
lm_box_qq_marginal_pos<-ggplot()+
  stat_qq(aes(sample=residuals(lm(formula=data$RF_boxcox~data$logIE_Pred),type='response')),size=3)+
  stat_qq_line(aes(sample=residuals(lm(formula=data$RF_boxcox~data$logIE_Pred),type='response')),lwd=2,color='red')+
  labs(y=expression(Marginal~Residual~Quantiles),
       x=expression(Theoretical~Quantiles),
       title = expression(lm~Box-Cox(RF)~vs~log[10](IE)~QQ~Plot))+
  theme_classic()

#Marginal Residual QQ Plots ESI-:
lm_box_qq_marginal_neg<-ggplot()+
  stat_qq(aes(sample=residuals(lm(formula=data_neg$RF_boxcox~data_neg$logIE_Pred),type='response')),size=3)+
  stat_qq_line(aes(sample=residuals(lm(formula=data_neg$RF_boxcox~data_neg$logIE_Pred),type='response')),lwd=2,color='red')+
  labs(y=expression(Marginal~Residual~Quantiles),
       x=expression(Theoretical~Quantiles),
       title = expression(lm~Box-Cox(RF)~vs~log[10](IE)~QQ~Plot))+
  theme_classic()

bn_obj_cond_pos <- bestNormalize(x=data$RF,
                                 allow_orderNorm = F,
                                 loo=T,
                                 norm_stat_fn = function(x) sqrt(mean(residuals(lmer(formula=x~data$logIE_Pred+(1|data$MS_Ready_DTXCID)))^2)))
bn_obj_cond_neg <- bestNormalize(x=data_neg$RF,
                                 allow_orderNorm = F,
                                 loo=T,
                                 norm_stat_fn = function(x) sqrt(mean(residuals(lmer(formula=x~data_neg$logIE_Pred+(1|data_neg$MS_Ready_DTXCID)))^2)))

lmer_box_pos = lmer(formula = data$RF_boxcox~data$logIE_Pred+(1|data$MS_Ready_DTXCID))
lmer_box_neg = lmer(formula = data_neg$RF_boxcox~data_neg$logIE_Pred+(1|data_neg$MS_Ready_DTXCID))

#Conditional Residual QQ Plots ESI+:
lmer_box_qq_conditional_pos<-ggplot()+
  stat_qq(aes(sample=residuals(lmer_box_pos,type='response')),size=3)+
  stat_qq_line(aes(sample=residuals(lmer_box_pos,type='response')),lwd=2,color='red')+
  labs(y=expression(Conditional~Residual~Quantiles),
       x=expression(Theoretical~Quantiles),
       title = expression(lmer~Box-Cox(RF)~vs~log[10](IE)~QQ~Plot))+
  theme_classic()

#Conditional Residual QQ Plots ESI-:
lmer_box_qq_conditional_neg<-ggplot()+
  stat_qq(aes(sample=residuals(lmer_box_neg,type='response')),size=3)+
  stat_qq_line(aes(sample=residuals(lmer_box_neg,type='response')),lwd=2,color='red')+
  labs(y=expression(Conditional~Residual~Quantiles),
       x=expression(Theoretical~Quantiles),
       title = expression(lmer~Box-Cox(RF)~vs~log[10](IE)~QQ~Plot))+
  theme_classic()

grid.arrange(lmer_box_qq_pos,lmer_box_qq_neg,ncol=2)
arrangeGrob(lmer_box_qq_pos,lmer_box_qq_neg,ncol=2)
ggsave(arrangeGrob(lmer_box_qq_pos,lmer_box_qq_neg,ncol=2),
       filename='lmer_residual_qqplots.jpeg',
       device = 'jpeg',
       dpi=1200,
       height=4,
       width=8)
boot_sample <- function(ID, data)
{
  tmpdf <- data[data$MS_Ready_DTXCID == ID,]$RF
  tmpdf <- tmpdf[complete.cases(tmpdf)]
  return(sample(tmpdf, 1, replace=T))
}
ggplot(data=glob_percentiles)+
  geom_histogram(aes(x=percentile_2p5th),
                 color='black',fill='blue',alpha=0.25)+
  theme_classic()+
  labs(x='Resampled 2.5th Percentile Estimate',
       y='Frequency')
  # scale_x_log10()

mc_percentiles <- function(data, descriptors, runs, samples)
{
  #store rng seeds for n runs for reproducibility
  #comment out for true random MC
  set.seed(16384)
  rf_data <- data[,descriptors]
  #remove NAs from RF list and preallocate data frames/matrices/vectors/etc:
  rf_data <- rf_data[complete.cases(rf_data$RF),]
  percentile_dist <- data.frame(percentile_2p5th = double(runs),
                                percentile_97p5th = double(runs))
  percentile_dist$RF_sample <- c()
  percentile_dist$MS_Ready_DTXCID <- c()
  percentile_2p5th <- vector('list',runs)
  percentile_97p5th <- percentile_2p5th
  rf_resample <- matrix(nrow = runs, ncol = samples)
  DTXCID <- matrix(nrow = runs, ncol = samples)
  unique_chems <- unique(rf_data$MS_Ready_DTXCID)
  print(paste0('Begin boostrap Simulation for N = ',runs, ' runs. '))
  for (i in 1:runs)
  {
    #get indexes of which chemicals to randomly sample:
    chem_idx_boot <- sample(length(unique_chems),
                            size=samples,
                            replace=T)
    #resampled list of chemicals and their indices:
    chem_resample <- unique_chems[chem_idx_boot]
    #choosing one random RF per MS_Ready_DTXCID:
    sample_out<-sapply(chem_resample, FUN = function(x) boot_sample(x,rf_data))
    DTXCID[i,] <- names(sample_out)
    rf_resample[i,] <- sample_out
    percentile_2p5th[[i]] <- quantile(rf_resample[i,],0.025)[[1]]
    percentile_97p5th[[i]] <- quantile(rf_resample[i,],0.975)[[1]]
    ##calculating geometric mean of RFs method:
    # for(j in 1:length(chem_resample))
    # {
    #   rf_resample[j,'RF'] <- 10^mean(log10(rf_data[rf_data$MS_Ready_DTXCID == chem_resample[[j]],]$RF))
    # }
    #calculate 1st percentile RF and 5th percentile RF (one-sided):
    percentile_dist$RF_resample[[i]] <- rf_resample[i,]
    percentile_dist$MS_Ready_DTXCID[[i]] <- DTXCID[i,]
    percentile_dist$percentile_2p5th[[i]] <- percentile_2p5th[[i]]
    percentile_dist$percentile_97p5th[[i]] <- percentile_97p5th[[i]]
    print(paste0(round(100*i/runs,2),'% complete.'))
  }
  return(percentile_dist)
}

#ESI+ data setup:
datadir <- 'C:/Users/lgroff/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/Data/ENTACT Semi Quant/Standard Mixtures/JRS Tracer Analysis'
setwd(datadir)
data <- read_xlsx('Supp_Table_Pos.xlsx')
data$RF <- data$Normalized_Intensity/data$Concentration
data <- data[complete.cases(data$RF),]
data$logIE_Pred <- as.numeric(NA)
posdata_2<-readRDS(file='posdata_2.rds')
data$logIE_Pred <- as.numeric(NA)
for(i in 1:nrow(posdata_2))
{
  data[data$DTXSID == posdata_2$DTXSID[[i]],]$logIE_Pred <- unique(posdata_2[posdata_2$DTXSID == posdata_2$DTXSID[[i]],]$logIE_Pred)
}
data<-data[complete.cases(data$logIE_Pred),]
data$RF_boxcox <- (data$RF^(0.285)-1)/0.285
#Repeat for ESI-:
data_neg <- read_xlsx('Supp_Table_Neg.xlsx')
data_neg$RF <- data_neg$Normalized_Intensity/data_neg$Concentration
data_neg <- data_neg[complete.cases(data_neg$RF),]
data_neg$logIE_Pred <- as.numeric(NA)
negdata_2<-readRDS(file='negdata_2.rds')
data_neg$logIE_Pred <- as.numeric(NA)
for(i in 1:nrow(negdata_2))
{
  data_neg[data_neg$DTXSID == negdata_2$DTXSID[[i]],]$logIE_Pred <- unique(negdata_2[negdata_2$DTXSID == negdata_2$DTXSID[[i]],]$logIE_Pred)
}
data_neg<-data_neg[complete.cases(data_neg$logIE_Pred),]

percentiles<-mc_percentiles(data,
                            runs=10000,
                            samples = length(unique(data[complete.cases(data$RF),]$MS_Ready_DTXCID)),
                            descriptors = c('MS_Ready_DTXCID','RF'))

#Write outer function to do k-fold CV with bootstrap sampled percentiles:
percentile_bootstrap_kfoldCV <- function(data, descriptors, folds, runs, samples)
{
  set.seed(16384)
  rf_data <- data[,descriptors]
  
  #remove NAs from RF list and preallocate data frames/columns:
  rf_data <- rf_data[complete.cases(rf_data$RF),]
  chemlist_CV <- unique(rf_data$MS_Ready_DTXCID)
  chem_idx_CV <- sample(length(unique(rf_data$MS_Ready_DTXCID)),
                        size=samples,
                        replace=F)
  #Set numbers of chemicals per fold:
  fold_size <- vector('list',5)
  #If mod division yields zero, all 5 folds are the same length
  if(length(chemlist_CV)%%folds == 0)
  {
    fold_size[1:5] <- length(chemlist_CV)/folds
  }
  else
  {
    #If mod division yields a number between 1-4, we add one sequentially to
    #the length of folds 1-4 until we've accounted for the remainder of
    #chemicals that make our N not divisible by 5.
    fold_size[1:5] <- (length(chemlist_CV)-length(chemlist_CV)%%folds)/5
      for(i in 1:(length(chemlist_CV)%%folds))
      {
        fold_size[[i]] <- fold_size[[i]]+1
      }
  }
  mc_out <- vector(mode='list',length=folds)
  train_folds <- mc_out
  test_folds <- mc_out
  idx <- 1
  for(i in 1:folds)
  {
    print(paste0('CV Run #',i))
    rf_test <- rf_data[rf_data$MS_Ready_DTXCID %in% chemlist_CV[chem_idx_CV[idx:(idx+fold_size[[i]])]],]
    rf_train <- rf_data[!(rf_data$MS_Ready_DTXCID %in% chemlist_CV[chem_idx_CV[idx:(idx+fold_size[[i]])]]),]
    idx <- idx+fold_size[[i]]
    train_percentiles<-mc_percentiles(data=rf_train,
                       descriptors,
                       runs,
                       samples)
    train_folds[[i]] <- rf_train
    test_folds[[i]] <- rf_test
    mc_out[[i]]<-train_percentiles
  }
  cv_out <- list(training = train_folds,
                 test = test_folds,
                 bootstrap = mc_out)
  return(cv_out)
}

#Generate calibration curve for ethionamide, for Figure 2:
eth_lm <- lm(data[data$Preferred_Name=='Ethionamide',],
             formula=log10(Normalized_Intensity)~log10(Concentration))
#Generate 95% prediction interval for ethionamide, for Figure 2:
eth_PI <- predict.lm(eth_lm,
                     interval='prediction',
                     level=0.95)
#Generate Calibration Curve with PI plot for Figure 2:
eth_fig2_gg <- ggplot(data=data[data$Preferred_Name=='Ethionamide',])+
  geom_point(aes(x=log10(Concentration),
                 y=log10(Normalized_Intensity)),
             size=3,
             shape='O')+
  geom_line(aes(x=log10(Concentration),
                y=eth_PI[,'fit']),
            lwd=1)+
  geom_line(aes(x=log10(Concentration),
                y=eth_PI[,'lwr']),
            lty=3,
            lwd=1)+
  geom_line(aes(x=log10(Concentration),
                y=eth_PI[,'upr']),
            lty=3,
            lwd=1)+
  geom_segment(aes(x=-0.5,
                   xend=-0.5,
                   y=6.95,
                   yend=0),
               color='red',
               lwd=1,
               lty=2)+
  geom_segment(aes(x=-3,
                   xend=-0.5,
                   y=6.95,
                   yend=6.95),
               color='red',
               lwd=1,
               lty=2)+
  coord_cartesian(xlim=c(-2,0),
                  ylim=c(5.5,7.5))+
  labs(x=expression(log[10]~(Concentration)),
       y=expression(log[10]~(Intensity)),
       title='Example Chemical')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(size=11),
        plot.title = element_text(hjust=0.5,size=10),
        axis.title.x = element_text(vjust=-0.5,size=10),
        axis.title.y = element_text(vjust=1.5,size=10))+
  annotate('label',
           x=-0.55,
           y=6,
           label=expression(italic(widehat(Conc))[0.975[CC]]),
           color='red',
           size=4)+
  annotate('text',
           x=-1.5,
           y=7.1,
           label=expression(italic(Y)[Obs]),
           color='red',
           size=4)

#save to a jpeg:
ggsave(plot = eth_fig2_gg,
       device = 'jpeg',
       dpi = 1200,
       height = 2.3,
       width = 2.3,
       filename = 'Fig2_CCPlot.jpeg')

# tick <- Sys.time()
CV_output_pos<-percentile_bootstrap_kfoldCV(data=data,
                                            descriptors=c('MS_Ready_DTXCID','RF'),
                                            folds=5,
                                            runs=10000,
                                            samples=length(unique(data$MS_Ready_DTXCID)))
# ticktock <- Sys.time()-tick
# print(ticktock)
# saveRDS(CV_output_neg,'CV_bootstrap_ESIneg.rds')
#Process results from CV Bootstrap RF Model:
train_test_ggs_pos <- vector(mode='list',length=5)
for(i in 1:5)
{
  train_test_ggs_pos[[i]]<-ggplot()+
    geom_histogram(data=CV_output[[1]][[i]],
                   aes(x=RF,fill='Training'),
                   color='black',
                   alpha=0.25)+
    geom_histogram(data=CV_output[[2]][[i]],
                   aes(x=RF,fill='Test'),
                   color='black',
                   alpha=0.25)+
    theme_classic()+
    # scale_x_log10()+
    labs(x='Response Factor',
         y='Frequency',
          title=paste0('CV Iteration ',i))+
    theme(axis.text.x = element_text(angle=45,hjust=1))
}

grid.arrange(train_test_ggs_pos[[1]],
             train_test_ggs_pos[[2]],
             train_test_ggs_pos[[3]],
             train_test_ggs_pos[[4]],
             train_test_ggs_pos[[5]],
             nrow=3,ncol=2)
#To save the figure, save an arrangeGrob of the plots to a
#variable that's callable by ggsave:
# train_test_gg_fig_neg<-arrangeGrob(train_test_ggs_neg[[1]],
#                                    train_test_ggs_neg[[2]],
#                                    train_test_ggs_neg[[3]],
#                                    train_test_ggs_neg[[4]],
#                                    train_test_ggs_neg[[5]],
#                                    nrow=3,ncol=2)
# ggsave(train_test_gg_fig_neg,
#        filename = 'CV_TrainTestSets_log10_ESI-.jpeg',
#        dpi=1200,
#        device='jpeg',
#        height=12,
#        width=10)

#Calculate outer 2.5th percentile RF estimates:
percentile_dist_ggs <- vector(mode='list',length=5)
q2p5_cv <- vector(mode='list',length=5)
q97p5_cv <- q2p5_cv
median_cv <- q2p5_cv
median_out_percent <- q2p5_cv
for(i in 1:5)
{
  q2p5_cv[[i]] <- quantile(CV_output_neg[[3]][[i]]$percentile_97p5th,0.025)
  q97p5_cv[[i]] <- quantile(CV_output_neg[[3]][[i]]$percentile_97p5th,0.975)
  median_cv[[i]] <- quantile(CV_output_neg[[3]][[i]]$percentile_97p5th,0.5)
  median_out_percent[[i]] <- length(CV_output_neg[[2]][[i]][CV_output_neg[[2]][[i]]$RF > median_cv[[i]],]$RF)/length(CV_output_neg[[2]][[i]]$RF)
  percentile_dist_ggs[[i]]<-ggplot()+
    geom_histogram(data=CV_output_neg[[3]][[i]],
                   aes(x=percentile_2p5th),
                   color='black',
                   fill='darkred',
                   alpha=0.25)+
    geom_vline(aes(xintercept=q2p5_cv[[i]]),
               lwd=1.5,lty=2)+
    geom_vline(aes(xintercept=q97p5_cv[[i]]),
               lwd=1.5,lty=2)+
    theme_classic()+
    labs(x='Bootstrapped 97.5th Percentile\nTraining Set RF',
         y='Frequency',
         title=paste0('CV Iteration ',i))+
    theme(axis.text.x = element_text(angle=45,hjust=1))
}
mean(unlist(median_out_percent))

grid.arrange(percentile_dist_ggs[[1]],
             percentile_dist_ggs[[2]],
             percentile_dist_ggs[[3]],
             percentile_dist_ggs[[4]],
             percentile_dist_ggs[[5]],
             nrow=2,ncol=3)
# percentile_dist_gg_fig<-arrangeGrob(percentile_dist_ggs[[1]],
#                                     percentile_dist_ggs[[2]],
#                                     percentile_dist_ggs[[3]],
#                                     percentile_dist_ggs[[4]],
#                                     percentile_dist_ggs[[5]],
#                                     nrow=3,ncol=2)
# ggsave(percentile_dist_gg_fig,
#        filename = 'CV_ESI+_2p5th_PercentileDistributions.jpeg',
#        dpi=1200,
#        device='jpeg',
#        height=12,
#        width=10)

#Plot global RF data with 2.5th percentile RF from bootstrap:
Default_global_RF_ESIpos<-ggplot(data=data)+
  geom_histogram(aes(x=RF),
                 na.rm=T,
                 fill='black',
                 color='white')+
  geom_vline(aes(xintercept = 1.37e6),
             lwd=1,
             lty=2,
             color='red')+
  scale_x_log10(limits = c(1e5,1e9))+
  theme_bw()+
  labs(x=expression(Response~Factor),
       y=expression(Frequency),
       title = 'Distribution of RFs')+
  annotate('label',
           x=1.39e6,
           y=210,
           label=expression(italic(widehat(RF)[0.025])),
           size=4,
           color='red')+
  theme(text = element_text(size=10.5),
        panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5,size=10),
        axis.title.x = element_text(vjust=-0.5,size=10),
        axis.title.y = element_text(vjust=1.75,size=10))

ggsave(plot=Default_global_RF_ESIpos,
       device='jpeg',
       dpi=1200,
       filename='Fig2_RFPlot.jpeg',
       height=2.3,
       width=2.3)

#store lambda transformations and shapiro-wilk normality test outputs:
for(i in 1:length(CV_output_neg[[3]]))
{
  CV_output_neg[[3]][[i]]$box_lambda <- as.numeric(NA)
  CV_output_neg[[3]][[i]]$W_shapiro <- as.numeric(NA)
  CV_output_neg[[3]][[i]]$P_shapiro <- as.numeric(NA)
  for(j in 1:nrow(CV_output_neg[[3]][[i]]))
  {
    tmpbox <- boxcox(CV_output_neg[[3]][[i]]$RF_sample[[j]],
                     standardize=F)
    tmpshap <- shapiro.test(tmpbox$x.t)
    CV_output_neg[[3]][[i]]$box_lambda[[j]] <- tmpbox$lambda
    CV_output_neg[[3]][[i]]$W_shapiro[[j]] <- tmpshap$statistic
    CV_output_neg[[3]][[i]]$P_shapiro[[j]] <- tmpshap$p.value
    print(paste0(100*((i-1)*10000+j)/50000,'% complete.'))
  }
}

#Find IE values corresponding to sampled RFs in training set:
for(i in 1:length(CV_output_neg[[1]]))
{
  CV_output_neg[[1]][[i]]$logIE_Pred <- as.numeric(NA)
  CV_output_neg[[1]][[i]]$RF_boxcox <- (CV_output_neg[[1]][[i]]$RF^mean(CV_output_neg[[3]][[i]]$box_lambda)-1)/mean(CV_output_neg[[3]][[i]]$box_lambda)
  for(j in 1:nrow(CV_output_neg[[1]][[i]]))
  {
    CV_output_neg[[1]][[i]]$logIE_Pred[[j]] <- unique(data_neg[data_neg$MS_Ready_DTXCID == CV_output_neg[[1]][[i]]$MS_Ready_DTXCID[[j]],]$logIE_Pred) 
  }
}
#Repeat for test set:
for(i in 1:length(CV_output_neg[[2]]))
{
  CV_output_neg[[2]][[i]]$logIE_Pred <- as.numeric(NA)
  CV_output_neg[[2]][[i]]$RF_boxcox <- (CV_output_neg[[2]][[i]]$RF^mean(CV_output_neg[[3]][[i]]$box_lambda)-1)/mean(CV_output_neg[[3]][[i]]$box_lambda)
  for(j in 1:nrow(CV_output_neg[[2]][[i]]))
  {
    CV_output_neg[[2]][[i]]$logIE_Pred[[j]] <- unique(data_neg[data_neg$MS_Ready_DTXCID == CV_output_neg[[2]][[i]]$MS_Ready_DTXCID[[j]],]$logIE_Pred)
  }
}
# saveRDS(CV_output_neg,'CV_output_neg_ESIpos.rds')


#PI function for bootmer()
mySumm <- function(.) {
  predict(., newdata=data, re.form=NULL)
}
#Fitted Value quantile estimate for PI bounds and fit from bootMer:
sumBoot <- function(merBoot) {
  return(
    data.frame(fit = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))),
               lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE))),
               upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE)))
    )
  )
}

#Calculate training set mixed effect models for input into bootMer:
train_mixmods_neg<-vector('list',5)
mixmods_boot_ESIneg <- vector('list',5)
for(i in 1:length(CV_output_neg[[1]]))
{
  tmpdf<-data.frame(MS_Ready_DTXCID = CV_output_neg[[1]][[i]]$MS_Ready_DTXCID,
                    RF_boxcox = CV_output_neg[[1]][[i]]$RF_boxcox,
                    logIE_Pred = CV_output_neg[[1]][[i]]$logIE_Pred)
  train_mixmods_neg[[i]] <- lmer(RF_boxcox ~ logIE_Pred + (1|MS_Ready_DTXCID),tmpdf)
  mixmods_boot_ESIneg[[i]] <- bootMer(x = train_mixmods_neg[[i]],
                                      FUN = mySumm,
                                      use.u = F,
                                      seed = 16384,
                                      type = 'parametric',
                                      nsim = 10000)
  print(paste0('Mixed Model Bootstrap #', i, ' Complete'))
}
write.csv(mixmods_boot_ESIpos[[1]]$data,file='ESIpos_TrainingSet1.csv')
apply(simulate(object=train_mixmods_neg[[1]],nsim=10,seed=16384,newdata=CV_output_neg[[1]][[1]],re.form=NULL),1,function(x) quantile(x,0.975))
data$RF_boxcox <- (data$RF^(.285)-1)/.285
ggplot(data=data)+
  geom_point(aes(y=RF_boxcox,
                 x=logIE_Pred))+
  geom_line(aes(x=logIE_Pred,
                y=mixmods_boot_ESIpos[[1]][["mle"]]$theta+mixmods_boot_ESIpos[[1]][["mle"]]$sigma+mixmods_boot_ESIpos[[1]][["mle"]]$beta[[1]]+mixmods_boot_ESIpos[[1]][["mle"]]$beta[[2]]*logIE_Pred))

# saveRDS(mixmods_boot_ESIpos,'mixmods_boot_ESIpos.rds')
# saveRDS(CV_output_neg,'CV_bootstrap_ESIneg.rds')
# saveRDS(mixmods_boot_ESIneg, 'mixmods_boot_ESIneg.rds')
#Generate and store mixed model predicition interval data:
mixmods_PI_ESIneg <- vector('list',5)
mixmods_gg_ESIneg <- vector('list',5)
for(i in 1:length(mixmods_boot_ESIneg))
{
  #Use sumBoot to generate PI:
  tmp_PI <- sumBoot(mixmods_boot_ESIneg[[i]])
  #Generate linear model objects 
  #from PI estimates on y-axis vs. logIE from Training sets:
  lwr_lm <- lm(formula=tmp_PI$lwr ~ CV_output_neg[[1]][[i]]$logIE_Pred)
  fit_lm <- lm(formula=tmp_PI$fit ~ CV_output_neg[[1]][[i]]$logIE_Pred)
  upr_lm <- lm(formula=tmp_PI$upr ~ CV_output_neg[[1]][[i]]$logIE_Pred)
  mixmods_PI_ESIneg[[i]] <- list(PI = tmp_PI,
                                 lwr_lm = lwr_lm,
                                 fit_lm = fit_lm,
                                 upr_lm = upr_lm)
  #Apply bootstrap prediction intervals to Test Sets using linear models:
  CV_output_neg[[2]][[i]]$mix_PI_fit <- mixmods_PI_ESIneg[[i]]$fit_lm$coefficients[[2]]*CV_output_neg[[2]][[i]]$logIE_Pred+mixmods_PI_ESIneg[[i]]$fit_lm$coefficients[[1]]
  CV_output_neg[[2]][[i]]$mix_PI_lwr <- mixmods_PI_ESIneg[[i]]$lwr_lm$coefficients[[2]]*CV_output_neg[[2]][[i]]$logIE_Pred+mixmods_PI_ESIneg[[i]]$lwr_lm$coefficients[[1]]
  CV_output_neg[[2]][[i]]$mix_PI_upr <- mixmods_PI_ESIneg[[i]]$upr_lm$coefficients[[2]]*CV_output_neg[[2]][[i]]$logIE_Pred+mixmods_PI_ESIneg[[i]]$upr_lm$coefficients[[1]]
  #Calculate % out of bounds (is one value, but stores in all rows of the list):
  CV_output_neg[[2]][[i]]$mix_percent_out_lwr <- 100*sum(CV_output_neg[[2]][[i]]$RF_boxcox < CV_output_neg[[2]][[i]]$mix_PI_lwr)/length(CV_output_neg[[2]][[i]]$RF_boxcox)
  CV_output_neg[[2]][[i]]$mix_percent_out_upr <- 100*sum(CV_output_neg[[2]][[i]]$RF_boxcox > CV_output_neg[[2]][[i]]$mix_PI_upr)/length(CV_output_neg[[2]][[i]]$RF_boxcox)
  #Store ggplot of Test Set with PI bounds and Percentages Outside PI:
  mixmods_gg_ESIneg[[i]] <- ggplot(data=CV_output_neg[[2]][[i]])+
    geom_point(aes(x=logIE_Pred,
                   y=RF_boxcox),
               size=3,
               shape='O')+
    geom_line(aes(x=logIE_Pred,
                  y=mix_PI_lwr),
              lwd=1.5)+
    geom_line(aes(x=logIE_Pred,
                  y=mix_PI_upr),
              lwd=1.5)+
    geom_line(aes(x=logIE_Pred,
                  y=mix_PI_fit),
              color='blue',
              lwd=1.5)+
    theme_classic()+
    labs(x=expression(log[10]~Ionization~Efficiency),
         y=expression((RF^lambda-1)/lambda),
         title=paste0('Test Set ',i,
                      '\n',round(unique(CV_output_neg[[2]][[i]]$mix_percent_out_lwr),2),
                      '% Outside Lower PI Bound',
                      '\n',round(unique(CV_output_neg[[2]][[i]]$mix_percent_out_upr),2),
                      '% Outside Upper PI Bound'))
}
# saveRDS(CV_output_neg,'CV_bootstrap_ESIneg.rds')

#apply mean CV prediction intervals to global dataset:
#Transform RF using mean lambda from Default RF CV bootstrap:
# data$RF_boxcox <- (data$RF^(0.285)-1)/0.285
#generate matrix of N x 5 for the 5 CV folds for fit, lwr and upr PI:
global_PI_lwr <- matrix(nrow = length(data$RF),ncol=5)
global_PI_fit <- global_PI_lwr
global_PI_upr <- global_PI_lwr

#Use linear regression coefficients to calculate the placement of the PI lines,
#Use global log(IE) as x-vector:
for(i in 1:5)
{
  global_PI_lwr[,i] <- mixmods_PI_ESIpos[[i]]$lwr_lm$coefficients[[2]]*data$logIE_Pred+mixmods_PI_ESIpos[[i]]$lwr_lm$coefficients[[1]]
  global_PI_fit[,i] <- mixmods_PI_ESIpos[[i]]$fit_lm$coefficients[[2]]*data$logIE_Pred+mixmods_PI_ESIpos[[i]]$fit_lm$coefficients[[1]]
  global_PI_upr[,i] <- mixmods_PI_ESIpos[[i]]$upr_lm$coefficients[[2]]*data$logIE_Pred+mixmods_PI_ESIpos[[i]]$upr_lm$coefficients[[1]]
}
#Average the 5 CV vectors across the matrix rows to get the mean fit, lwr, upr:
data$mix_fit <- apply(global_PI_fit,1,function(x) mean(x))
data$mix_lwrPI <- apply(global_PI_lwr,1,function(x) mean(x))
data$mix_uprPI <- apply(global_PI_upr,1,function(x) mean(x))

data$mix_uprPI_smooth <- 195.1*data$logIE_Pred+203.1
data$mix_lwrPI_smooth <- 195.7*data$logIE_Pred-308.7
data$mix_fit_smooth <- 195.4*data$logIE_Pred-52.9
glob_gg_ESIpos <- ggplot(data=data)+
  geom_point(aes(x=logIE_Pred,
                 y=RF_boxcox),
             size=3,
             shape='O')+
  geom_line(aes(x=logIE_Pred,
                y=mix_lwrPI_smooth),
            lwd=1.5)+
  geom_line(aes(x=logIE_Pred,
                y=mix_uprPI_smooth),
            lwd=1.5)+
  geom_line(aes(x=logIE_Pred,
                y=mix_fit_smooth),
            color='blue',
            lwd=1.5)+
  theme_classic()+
  labs(x=expression(log[10]~Ionization~Efficiency),
       y=expression((RF^lambda-1)/lambda),
       title=paste0('Bootstrapped Global Dataset',
                    '\n',round(100*sum(data$RF_boxcox < data$mix_lwrPI_smooth)/length(data$RF_boxcox),2),
                    '% Outside Lower PI Bound',
                    '\n',round(100*sum(data$RF_boxcox > data$mix_uprPI_smooth)/length(data$RF_boxcox),2),
                    '% Outside Upper PI Bound'))

data_neg$mix_uprPI_smooth <- 0.1084*data_neg$logIE_Pred+7.8276
data_neg$mix_lwrPI_smooth <- 0.1065*data_neg$logIE_Pred+6.7936
data_neg$mix_fit_smooth <- 0.1072*data_neg$logIE_Pred+7.3108
glob_gg_ESIneg <- ggplot(data=data_neg)+
  geom_point(aes(x=logIE_Pred,
                 y=RF_boxcox),
             size=3,
             shape='O')+
  geom_line(aes(x=logIE_Pred,
                y=mix_lwrPI_smooth),
            lwd=1.5)+
  geom_line(aes(x=logIE_Pred,
                y=mix_uprPI_smooth),
            lwd=1.5)+
  geom_line(aes(x=logIE_Pred,
                y=mix_fit_smooth),
            color='blue',
            lwd=1.5)+
  theme_classic()+
  labs(x=expression(log[10]~Ionization~Efficiency),
       y=expression((RF^lambda-1)/lambda),
       title=paste0('Bootstrapped Global Dataset',
                    '\n',round(100*sum(data_neg$RF_boxcox < data_neg$mix_lwrPI_smooth)/length(data_neg$RF_boxcox),2),
                    '% Outside Lower PI Bound',
                    '\n',round(100*sum(data_neg$RF_boxcox > data_neg$mix_uprPI_smooth)/length(data_neg$RF_boxcox),2),
                    '% Outside Upper PI Bound'))

grid.arrange(mixmods_gg_ESIpos[[1]],
             mixmods_gg_ESIpos[[2]],
             mixmods_gg_ESIpos[[3]],
             mixmods_gg_ESIpos[[4]],
             mixmods_gg_ESIpos[[5]],
             glob_gg_ESIpos,
             nrow=3,ncol=2)
# gg_mixmod_grid_ESIpos<-arrangeGrob(mixmods_gg_ESIpos[[1]],
#                                    mixmods_gg_ESIpos[[2]],
#                                    mixmods_gg_ESIpos[[3]],
#                                    mixmods_gg_ESIpos[[4]],
#                                    mixmods_gg_ESIpos[[5]],
#                                    glob_gg_ESIpos,
#                                    nrow=3,ncol=2)
# ggsave(gg_mixmod_grid_ESIpos,
#        filename = 'CV_Mixmod_Regression_PI95_ESIpos.jpeg',
#        dpi=1200,
#        device='jpeg',
#        height=12,
#        width=10)

# data$RF_boxcox <- (data$RF^(0.285)-1)/0.285
RF_vs_IE_plot<-ggplot(data=data)+
  geom_point(aes(x=logIE_Pred,
                 y=RF_boxcox),
             shape='O',
             size=2)+
  geom_line(aes(x=logIE_Pred,
                y=mix_uprPI_smooth),
            lty=2,
            lwd=1)+
  geom_line(aes(x=logIE_Pred,
                y=mix_lwrPI_smooth),
            lty=2,
            lwd=1)+
  geom_segment(aes(y= 440,
                   yend=440,
                   x = 0,
                   xend=3.8),
               lwd=1,
               lty=2,
               color='red')+
  geom_segment(aes(y= 440,
                   yend=-100,
                   x = 3.8,
                   xend=3.8),
               lwd=1,
               lty=2,
               color='red')+
  coord_cartesian(ylim=c(0,1100),
                  xlim=c(1.25,4.5))+
  labs(x=expression(log[10](Predicted~IE)),
       y=expression(Transformed~RF~(italic(tRF))),
       title = 'RF vs. IE Calibration',
       size=8)+
  theme(text = element_text(size=11),
        panel.grid = element_blank(),
        panel.background = element_rect(fill='white',color='black'),
        plot.title = element_text(hjust=0.5, size=10),
        axis.title.x = element_text(vjust=-0.5, size=10),
        axis.title.y = element_text(vjust=1,size=10))+
  annotate('label',
           x=2.0,
           y=585,
           label=expression(italic(widehat(tRF)[0.025[IE]])),
           size=4,
           color='red')
 
ggsave(plot=RF_vs_IE_plot,
       device='jpeg',
       dpi=1200,
       height=2.3,
       width=2.3,
       filename='Fig2_tRFplot.jpeg')

#Plot and save RF and Boxcox(RF) Q-Q Plots for full dataset.
gg_RF_qq<-ggplot(data=data)+
  stat_qq(aes(sample=RF),
          na.rm=T,
          size=3,
          shape='O')+
  stat_qq_line(aes(sample=RF),
               na.rm=T,
               lwd=1.5,
               color='red')+
  theme_classic()+
  labs(x='Theoretical Quantiles',
       y='RF Quantiles',
       title=paste0('ESI+ RF QQ Plot (n = ',nrow(data),')'))

gg_boxcox_qq <- ggplot(data=data)+
  stat_qq(aes(sample=RF_boxcox),
          na.rm=T,
          size=3,
          shape='O')+
  stat_qq_line(aes(sample=RF_boxcox),
               na.rm=T,
               lwd=1.5,
               color='red')+
  theme_classic()+
  labs(x='Theoretical Quantiles',
       y=expression((RF^lambda-1)/lambda~Quantiles),
       title=paste0('ESI+ Box-Cox RF QQ Plot (n = ',nrow(data),')'))

#Repeat for Negative Mode:
gg_RF_qq_Neg<-ggplot(data=data_neg)+
  stat_qq(aes(sample=RF),
          size=3,shape='O')+
  stat_qq_line(aes(sample=RF),
               lwd=1.5,color='red')+
  theme_classic()+
  labs(x='Theoretical Quantiles',
       y='RF Quantiles',
       title=paste0('ESI- Normal QQ Plot (n = ',nrow(data_neg),')'))

gg_boxcox_qq_Neg <-ggplot(data=data_neg)+
  stat_qq(aes(sample=RF_boxcox),
          size=3,shape='O')+=
  stat_qq_line(aes(sample=RF_boxcox),
               lwd=1.5,color='red')+
  theme_classic()+
  labs(x='Theoretical Quantiles',
       y=expression((RF^lambda-1)/lambda~Quantiles),
       title=paste0('ESI- Box-Cox(RF) QQ Plot (n = ',nrow(data_neg),')'))

# grid.arrange(gg_RF_qq,gg_boxcox_qq,nrow=1)
# gg_RF_boxcox_qqs<-arrangeGrob(gg_RF_qq,gg_boxcox_qq,nrow=1)
ggsave(arrangeGrob(gg_RF_qq_Neg,
                   gg_boxcox_qq_Neg,
                   lm_box_qq_marginal_neg,
                   lmer_box_qq_conditional_neg,
                   nrow=2,ncol=2),
       filename = 'ENTACT_ESI-_RF_Transform_QQPlots.jpeg',
       dpi=1200, device='jpeg',
       height=8,width=8)

# grid.arrange(gg_RF_qq_Neg,gg_boxcox_qq_Neg,nrow=1)
# gg_RF_boxcox_qqs_Neg<-arrangeGrob(gg_RF_qq_Neg,gg_boxcox_qq_Neg,nrow=1)
# ggsave(gg_RF_boxcox_qqs_Neg,
#        filename = 'ENTACT_ESI-_RF_Transform_QQPlots.jpeg',
#        dpi=1200,device='jpeg',
#        height=4,width=8)

for(i in 1:length(CV_Test_Sets_Neg))
{
  tmpIE_neg <- vector('list',length(CV_Test_Sets_Neg[[i]]$RF))
  for(j in 1:length(CV_Test_Sets_Neg[[i]]$RF))
  {
    tmpIE_neg[[j]] <- negdata_2[negdata_2$New_ID == CV_Test_Sets_Neg[[i]]$New_ID[[j]],]$logIE_Pred
  }
  CV_Test_Sets_Neg[[i]]$logIE_Pred<-unlist(tmpIE_neg)
}

# saveRDS(CV_Training_Sets_Neg,'CV_Training_Sets_Neg.rds')
# saveRDS(CV_Test_Sets_Neg,'CV_Test_Sets_Neg.rds')
# saveRDS(CV_bootstrap_IE_Neg,'CV_bootstrap_IE_Neg.rds')

posdata_1$Conc_CC <- -16384
posdata_1$Conc_CC_upr <- -16384
posdata_1$Conc_CC_lwr <- -16384
posdata_1$logIE_Pred <- -16384
for (i in 1:length(unique(posdata_1$New_ID)))
{ 
  tmpdf <- posdata_1[posdata_1$New_ID == unique(posdata_1$New_ID)[[i]],]
  tmpdf <- tmpdf[complete.cases(tmpdf$Normalized_Intensity),]
  tmplm <- lm(tmpdf,
              formula = log10(Normalized_Intensity) ~ log10(Concentration))
  posdata_1[posdata_1$New_ID == unique(posdata_1$New_ID)[[i]],]$Conc_CC <- 10^(log10(posdata_1[posdata_1$New_ID == unique(posdata_1$New_ID)[[i]],]$Normalized_Intensity)-tmplm$coefficients[[1]])/tmplm$coefficients[[2]]
  if(nrow(tmpdf)>3)
  {
    tmp_PI <- matrix(ncol=3,nrow=nrow(tmpdf))
    for (j in 1:nrow(tmp_PI))
    {
      tmpcal <- try(invest(tmplm,
                           interval='inversion',
                           level=0.95,
                           maxiter=1e4,
                           y0 = log10(tmpdf$Normalized_Intensity[[j]])),
                    silent=T)
      if(class(tmpcal) != 'try-error')
      {
        tmp_PI[j,1] <- tmpcal$lower
        tmp_PI[j,2] <- tmpcal$estimate
        tmp_PI[j,3] <- tmpcal$upper
        posdata_1[posdata_1$New_ID == unique(posdata_1$New_ID)[[i]],]$Conc_CC_lwr[[j]] <- tmp_PI[j,1]
        posdata_1[posdata_1$New_ID == unique(posdata_1$New_ID)[[i]],]$Conc_CC_upr[[j]] <- tmp_PI[j,3]
      }
      else
      {
        posdata_1[posdata_1$New_ID == unique(posdata_1$New_ID)[[i]],]$Conc_CC_lwr[[j]] <- NA
        posdata_1[posdata_1$New_ID == unique(posdata_1$New_ID)[[i]],]$Conc_CC_lwr[[j]] <- NA  
      }
    }
  }
  else
  {
    posdata_1[posdata_1$New_ID == unique(posdata_1$New_ID)[[i]],]$Conc_CC_lwr <- NA
    posdata_1[posdata_1$New_ID == unique(posdata_1$New_ID)[[i]],]$Conc_CC_upr <- NA    
  }
  posdata_1[posdata_1$New_ID == unique(posdata_1$New_ID)[[i]],'logIE_Pred'] <- unique(posdata_2[posdata_2$New_ID == unique(posdata_1$New_ID)[[i]],]$logIE_Pred)
}
posdata_1$EQ_CCupr_CCest <- posdata_1$Conc_CC_upr/posdata_1$Conc_CC

#Concentration estimates for ESI+
#Preallocate Calibration curve, upper/lower PI columns
data$Conc_CC <- as.numeric(NA)
data$Conc_CC_upr <- as.numeric(NA)
data$Conc_CC_lwr <- as.numeric(NA)
for (i in 1:length(unique(data$MS_Ready_DTXCID)))
{ 
  #subset data frame for one chemical's data:
  tmpdf <- data[data$MS_Ready_DTXCID == unique(data$MS_Ready_DTXCID)[[i]],]
  tmpdf <- tmpdf[complete.cases(tmpdf$Normalized_Intensity),]
  #generate calibration curve linear model:
  tmplm <- lm(tmpdf,
              formula = log10(Normalized_Intensity) ~ log10(Concentration))
  #generate cal. curve estimates:
  data[data$MS_Ready_DTXCID == unique(data$MS_Ready_DTXCID)[[i]],]$Conc_CC <- 10^((log10(data[data$MS_Ready_DTXCID == unique(data$MS_Ready_DTXCID)[[i]],]$Normalized_Intensity)-tmplm$coefficients[[1]])/tmplm$coefficients[[2]])
  #Given inverse calibration degrees of freedom requirements, only attempt
  #inverse calibration on N > 3 data points:
  if(nrow(tmpdf)>3)
  {
    #Store temporary matrix for PI
    tmp_PI <- matrix(ncol=3,nrow=nrow(tmpdf))
    for (j in 1:nrow(tmp_PI))
    {
      #encapsulate invest in a silent try statement so that failures don't break loop:
      # 1e4 iterations to try for convergence, 95% PI, y0 is for each intensity,
      # inversion interval is the x-axis equivalent PI from abscissas of y-axis PI:
      tmpcal <- try(calibrate(tmplm,
                              interval='inversion',
                              level=0.95,
                              maxiter=1e4,
                              y0 = log10(tmpdf$Normalized_Intensity[[j]])),
                    silent=T)
      #store results if the invest function ran successfully (i.e., not a try-error):
      if(class(tmpcal) != 'try-error')
      {
        # print(j)
        tmp_PI[j,1] <- 10^tmpcal$lower
        tmp_PI[j,2] <- 10^tmpcal$estimate
        tmp_PI[j,3] <- 10^tmpcal$upper
        data[data$MS_Ready_DTXCID == unique(data$MS_Ready_DTXCID)[[i]],]$Conc_CC_lwr[[j]] <- tmp_PI[j,1]
        data[data$MS_Ready_DTXCID == unique(data$MS_Ready_DTXCID)[[i]],]$Conc_CC_upr[[j]] <- tmp_PI[j,3]
      }
    }
  }
}
#Store error quotient of CC_upper_PI/CC_estimate:
data$EQ_CCupr_CCest <- data$Conc_CC_upr/data$Conc_CC
write.csv(data,
          'posdata_calibrate_CC.csv')
#ESI+ Default RF bounded concentrations:
#CV-averaged 95% PI bounds for RF bootstrap:
boot_RF_2p5_pos <- 1.37e6
boot_RF_97p5_pos <- 2.43e8
#Generate Conc_RF estimates from the bootstrap RF bounds:
data$Conc_RF_upr <- data$Normalized_Intensity/boot_RF_2p5_pos
data$Conc_RF_lwr <- data$Normalized_Intensity/boot_RF_97p5_pos
data$EQ_RF_CC <- data$Conc_RF_upr/data$Conc_CC

tmplm_pos <- lm(data[data$Preferred_Name == "1,3-Diphenylguanidine",],
                formula=log10(Normalized_Intensity)~log10(Concentration))
tmp_PI_pos <- predict.lm(tmplm_pos,
                         interval='prediction',
                         level=0.95)
# tmpcal <- calibrate(tmplm,
#                     interval='inversion',
#                     level=0.95,
#                     y0=log10(data[data$Preferred_Name == "1,3-Diphenylguanidine",]$Normalized_Intensity[[4]]))

worst_CC_ESIpos<-ggplot(data=data[data$Preferred_Name == "1,3-Diphenylguanidine",])+
  geom_point(aes(x=log10(Concentration),
                 y=log10(Normalized_Intensity)),
             shape='O',size=3)+
  labs(x=expression(log[10]~Concentration),
       y=expression(log[10]~Intensity),
       title='1,3-Diphenylguanidine')+
  geom_abline(aes(slope=tmplm_pos$coefficients[[2]],
                  intercept=tmplm_pos$coefficients[[1]]),
              lwd=1.5, color='blue')+
  theme_classic()+
  geom_ribbon(aes(x=log10(Concentration),
                  ymin=tmp_PI_pos[,'lwr'],
                  ymax=tmp_PI_pos[,'upr']),
              alpha=0.25)+
  annotate('label',
           x=-1.46,
           y=8.15,
           label=paste0('ESI+\n R² =',round(summary(tmplm_pos)$r.squared,4),
                        '\nslope = ',round(tmplm_pos$coefficients[[2]],3),
                        '\nintercept =',round(tmplm_pos$coefficients[[1]],3)))+
  theme(text = element_text(size=15),
        panel.grid = element_blank(),
        panel.background = element_rect(fill='white',color='black'),
        plot.title = element_text(hjust=0.5),
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=1))

tmplm_neg <- lm(data_neg[data_neg$Preferred_Name == "Tris(1,3-dichloro-2-propyl) phosphate",],
            formula=log10(Normalized_Intensity)~log10(Concentration))
tmp_PI_neg <- predict.lm(tmplm_neg,
                     interval='prediction',
                     level=0.95)

worst_CC_ESIneg<-ggplot(data=data_neg[data_neg$Preferred_Name == "Tris(1,3-dichloro-2-propyl) phosphate",])+
  geom_point(aes(x=log10(Concentration),
                 y=log10(Normalized_Intensity)),
             shape='O',size=3)+
  labs(x=expression(log[10]~Concentration),
       y=expression(log[10]~Intensity),
       title="Tris(1,3-dichloro-2-propyl) phosphate")+
  geom_abline(aes(slope=tmplm_neg$coefficients[[2]],
                  intercept=tmplm_neg$coefficients[[1]]),
              lwd=1.5, color='blue')+
  theme_classic()+
  geom_ribbon(aes(x=log10(Concentration),
                  ymin=tmp_PI_neg[,'lwr'],
                  ymax=tmp_PI_neg[,'upr']),
              alpha=0.25)+
  annotate('label',
           x=-1.44,
           y=6.025,
           label=paste0('ESI-\nR² =',round(summary(tmplm_neg)$r.squared,4),
                        '\nslope = ',round(tmplm_neg$coefficients[[2]],3),
                        '\nintercept =',round(tmplm_neg$coefficients[[1]],3)))+
  theme(text = element_text(size=15),
        panel.grid = element_blank(),
        panel.background = element_rect(fill='white',color='black'),
        plot.title = element_text(hjust=0.5),
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=1))
ggplot()+
  geom_point(aes(y=seq(from=0,to=100,length.out = length(data$EQ_RF_CC)),
             x=sort(data$EQ_RF_CC)),
             size=4,shape='O',color='darkgreen')+
  geom_point(aes(y=seq(from=0,to=100,length.out = length(data$EQ_IE_CC)),
                 x=sort(data$EQ_IE_CC)),
             size=4,shape='O',color='purple')+
  geom_point(aes(y=seq(from=0,to=100,length.out = length(data[complete.cases(data$EQ_CCupr_CCest),]$EQ_CCupr_CCest)),
                 x=sort(data[complete.cases(data$EQ_CCupr_CCest),]$EQ_CCupr_CCest)),
             size=4,shape='O',color='blue')+
  scale_x_log10(limits=c(1e-1,1e3))+
  labs(x='Error Quotient',
       y='Percentile Rank (%)')

# worst_CCs_ESI<-arrangeGrob(worst_CC_ESIpos,
#                            worst_CC_ESIneg,
#                            nrow=1,
#                            ncol=2)
# ggsave(worst_CCs_ESI,
#        height=4,
#        width=10,
#        dpi=1200,
#        device='jpeg',
#        filename='ENTACT_1B_WorstCalCurves.jpeg')

#Calculate ESI+ Conc_IE_95
data$Conc_IE_upr <- as.numeric(NA)
data$imputed <- as.logical(NA)
cnt<-0
for (i in 1:nrow(data))
{
  #Taking min(tRF) as close to the LOD, 
  #If IE-calculated RF is <= min(tRF)/sqrt(2), impute min(tRF)/sqrt(2):
  min_impute_pos <- ((min(data$RF)/sqrt(2))^(0.285)-1)/0.285
  if(data$mix_lwrPI_smooth[[i]] <= min_impute_pos)
  {
    cnt<-cnt+1
    data$imputed[[i]] <- T
    #convert RF back to linear space and calculate concentration as I/RF:
    min_impute_pos_lin <- (min_impute_pos*0.285+1)^(1/0.285)
    data$Conc_IE_upr[[i]] <- data$Normalized_Intensity[[i]]/min_impute_pos_lin  
  }
  else
  {
    data$imputed[[i]] <- F
    #convert RF back to linear space and calculate concentration as I/RF:
    tmp_RF_IE <- (data$mix_lwrPI_smooth[[i]]*0.285+1)^(1/0.285)
    data$Conc_IE_upr[[i]] <- data$Normalized_Intensity[[i]]/tmp_RF_IE
  }
}
data$Conc_IE_lwr = data$Normalized_Intensity/((data$mix_uprPI_smooth*0.285+1)^(1/0.285))
#Print Count and percentage of imputed RFs:
print(paste0(cnt,'/',length(data$RF),' RFs < min(tRF)/sqrt(2)'))
print(paste0(round(100*cnt/length(data$RF_IE_lwr),2),'% of RFs < min(tRF)/sqrt(2)'))
data$EQ_IE_CC <- data$Conc_IE_upr/data$Conc_CC
data$EQ_RF_IE <- data$Conc_RF_upr/data$Conc_IE_upr

############
#Imputed Response Factors are not stored
#just the resulting concentrations in conc_IE_upr
############

min(CV_output_neg$bootstrap[[1]]$percentile_2p5th)
median(CV_output_neg$bootstrap[[1]]$percentile_2p5th)
max(CV_output_neg$bootstrap[[1]]$percentile_2p5th)
sum((unique(data$MS_Ready_DTXCID) %in% unique(data_neg$MS_Ready_DTXCID)))
min(glob_percentiles_neg$percentile_97p5th)
max(glob_percentiles_neg$percentile_97p5th)

impute_gg_pos <- ggplot(data)+
  geom_point(aes(y=RF_boxcox,x=logIE_Pred),
             shape='O',size=3)+
  geom_hline(aes(yintercept=min_impute_pos),
             lty=2,lwd=1.5)+
  geom_line(aes(x=logIE_Pred,y=mix_lwrPI_smooth),
            lwd=1.5)+
  #Color code imputed data points in red:
  geom_point(data=data[data$mix_lwrPI_smooth <= min_impute_pos,],
             aes(x=logIE_Pred,y=RF_boxcox),
             color='red')+
  annotate('label',
           x=4,
           y=min_impute_pos,
           label=expression(((RF[min]/sqrt(2))^lambda-1)/lambda))+
  labs(x=expression(log[10](Predicted~IE)),
       y=expression((RF^lambda-1)/lambda),
       title='ESI+')+
  theme_classic()+
  geom_hline(aes(yintercept=(1.39e6^(.285)-1)/.285),
             lty=3,lwd=1.5)


impute_gg_neg <- ggplot(data=data_neg)+
  geom_point(aes(y=RF_boxcox,x=logIE_Pred),
             shape='O',size=3)+
  geom_hline(aes(yintercept=min_impute_neg),
             lty=2,lwd=1.5)+
  geom_line(aes(x=logIE_Pred,y=mix_lwrPI_smooth),
            lwd=1.5)+
  #Color code imputed data_neg points in red:
  geom_point(data=data_neg[data_neg$mix_lwrPI_smooth <= min_impute_neg,],
             aes(x=logIE_Pred,y=RF_boxcox),
             color='red')+
  annotate('label',
           x=2,
           y=1.005*min_impute_neg,
           label=expression(((RF[min]/sqrt(2))^lambda-1)/lambda))+
  labs(x=expression(log[10](Predicted~IE)),
       y=expression((RF^lambda-1)/lambda),
       title='ESI-')+
  theme_classic()

RF_IE_gg<-ggplot()+
  geom_point(aes(y=data$EQ_RF_IE,
                 x=(data$logIE_Pred-mean(data$logIE_Pred))/sd(data$logIE_Pred),
                 color='ESI+'),
             alpha=0.25,
             size=3)+
  geom_point(aes(y=data_neg$EQ_RF_IE,
                 x=(data_neg$logIE_Pred-mean(data_neg$logIE_Pred))/sd(data_neg$logIE_Pred),
                 color='ESI-'),
             alpha=0.25,
             size=3)+
  geom_vline(aes(xintercept=-3))+
  geom_vline(aes(xintercept=-2))+
  geom_vline(aes(xintercept=-1))+
  geom_vline(aes(xintercept=0))+
  geom_vline(aes(xintercept=1))+
  geom_vline(aes(xintercept=2))+
  geom_vline(aes(xintercept=3))+
  labs(x=expression(Z~Score~log[10]~(Predicted~IE)),
       y=expression(EQ[RF]/EQ[IE]))+
  theme_classic()+
  scale_y_log10(limits = c(0.1,100))+
  scale_x_continuous(limits=c(-3.25,3.25),
                     breaks=c(-3,-2,-1,0,1,2,3),
                     labels = c(-3,-2,-1,0,1,2,3))+
  theme(legend.title=element_blank(),
        legend.position = c(0.15,0.9),
        legend.background = element_blank())
RF_IE_gg

EQ_IE_vs_RF_gg <- ggplot()+
  geom_point(aes(y=data$EQ_RF_CC,
                 x=data$EQ_IE_CC,
                 color='ESI+'),
             shape='O',size=3)+
  geom_point(aes(y=data_neg$EQ_RF_CC,
                 x=data_neg$EQ_IE_CC,
                 color='ESI-'),
             shape='O',size=3)+
  geom_abline(aes(slope=1,
                  intercept=0),
              lwd=1.5)+
  scale_x_log10(limits=c(0.01,1000),
                breaks = c(1e-2,1e-1,1e0,1e1,1e2,1e3),
                labels = c(1e-2,1e-1,1e0,1e1,1e2,1e3))+
  scale_y_log10(limits=c(0.01,1000),
                breaks = c(1e-2,1e-1,1e0,1e1,1e2,1e3),
                labels = c(1e-2,1e-1,1e0,1e1,1e2,1e3))+
  labs(y=expression(widehat(Conc)[0.975[RF]]/widehat(Conc)[CC]),
       x=expression(widehat(Conc)[0.975[IE]]/widehat(Conc)[CC]))+
  theme_classic()+
  theme(legend.title = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.9,0.15))
EQ_IE_vs_RF_gg

# impute_plots <- arrangeGrob(impute_gg_pos,
#                             impute_gg_neg,
#                             ncol=2,nrow=1)
# ggsave(impute_plots,
#        dpi = 1200,
#        device = 'jpeg',
#        filename = 'ENTACT_Imputation_Plot.jpeg',
#        height = 4, width = 10)


#Concentration estimates for ESI-
#Preallocate Calibration curve, upper/lower PI columns
data_neg$Conc_CC <- as.numeric(NA)
data_neg$Conc_CC_upr <- as.numeric(NA)
data_neg$Conc_CC_lwr <- as.numeric(NA)
for (i in 1:length(unique(data_neg$MS_Ready_DTXCID)))
{ 
  #subset data_neg frame for one chemical's data_neg:
  tmpdf <- data_neg[data_neg$MS_Ready_DTXCID == unique(data_neg$MS_Ready_DTXCID)[[i]],]
  tmpdf <- tmpdf[complete.cases(tmpdf$Normalized_Intensity),]
  #generate calibration curve linear model:
  tmplm <- lm(tmpdf,
              formula = log10(Normalized_Intensity) ~ log10(Concentration))
  #generate cal. curve estimates:
  data_neg[data_neg$MS_Ready_DTXCID == unique(data_neg$MS_Ready_DTXCID)[[i]],]$Conc_CC <- 10^((log10(data_neg[data_neg$MS_Ready_DTXCID == unique(data_neg$MS_Ready_DTXCID)[[i]],]$Normalized_Intensity)-tmplm$coefficients[[1]])/tmplm$coefficients[[2]])
  #Given inverse calibration degrees of freedom requirements, only attempt
  #inverse calibration on N > 3 data_neg points:
  if(nrow(tmpdf)>3)
  {
    #Store temporary matrix for PI
    tmp_PI <- matrix(ncol=3,nrow=nrow(tmpdf))
    for (j in 1:nrow(tmp_PI))
    {
      #encapsulate invest in a silent try statement so that failures don't break loop:
      # 1e4 iterations to try for convergence, 95% PI, y0 is for each intensity,
      # inversion interval is the x-axis equivalent PI from abscissas of y-axis PI:
      tmpcal <- try(calibrate(tmplm,
                              interval='inversion',
                              level=0.95,
                              maxiter=1e4,
                              y0 = log10(tmpdf$Normalized_Intensity[[j]])),
                    silent=T)
      #store results if the invest function ran successfully (i.e., not a try-error):
      if(class(tmpcal) != 'try-error')
      {
        # print(j)
        tmp_PI[j,1] <- 10^tmpcal$lower
        tmp_PI[j,3] <- 10^tmpcal$upper
        data_neg[data_neg$MS_Ready_DTXCID == unique(data_neg$MS_Ready_DTXCID)[[i]],]$Conc_CC_lwr[[j]] <- tmp_PI[j,1]
        data_neg[data_neg$MS_Ready_DTXCID == unique(data_neg$MS_Ready_DTXCID)[[i]],]$Conc_CC_upr[[j]] <- tmp_PI[j,3]
      }
    }
  }
}
#Store error quotient of CC_upper_PI/CC_estimate:
data_neg$EQ_CCupr_CCest <- data_neg$Conc_CC_upr/data_neg$Conc_CC
write.csv(data,'posdata_withimputations.csv')
#ESI+ Default RF bounded concentrations:
#CV-averaged 95% PI bounds for RF bootstrap:
boot_RF_2p5_neg <- 2.77e5
boot_RF_97p5_neg <- 4.52e7
#Generate Conc_RF estimates from the bootstrap RF bounds:
data_neg$Conc_RF_upr <- data_neg$Normalized_Intensity/boot_RF_2p5_neg
data_neg$Conc_RF_lwr <- data_neg$Normalized_Intensity/boot_RF_97p5_neg
data_neg$EQ_RF_CC <- data_neg$Conc_RF_upr/data_neg$Conc_CC

#Compute the percentage of cases that yield imputations where RF_2.5 > RF_IE_2.5

#Calculate ESI- Conc_IE_95
data_neg$Conc_IE_upr <- as.numeric(NA)
data_neg$imputed <- as.logical(NA)
cnt<-0
for (i in 1:nrow(data_neg))
{
  #If IE-calculated RF is < 2.5th percentile default RF, use default RF to calc Conc_IE:
  min_impute_neg <- ((min(data_neg$RF)/sqrt(2))^(-0.106)-1)/-0.106
  if(data_neg$mix_lwrPI_smooth[[i]] <= min_impute_neg)
  {
    data_neg$imputed[[i]] <- T
    cnt<-cnt+1
    min_impute_neg_lin <- (min_impute_neg*(-0.106)+1)^(1/-0.106)
    data_neg$Conc_IE_upr[[i]] <- data_neg$Normalized_Intensity[[i]]/min_impute_neg_lin  
  }
  else
  {
    data_neg$imputed <- F
    tmp_RF_IE <- (data_neg$mix_lwrPI_smooth[[i]]*-0.106+1)^(1/-0.106)
    data_neg$Conc_IE_upr[[i]] <- data_neg$Normalized_Intensity[[i]]/tmp_RF_IE
  }
  #Print Count and percentage of imputed RFs:
}
data_neg$Conc_IE_lwr = data_neg$Normalized_Intensity/((data_neg$mix_uprPI_smooth*-0.106+1)^(1/-0.106))
print(paste0(cnt,'/',length(data_neg$RF),' RFs < 2.5th percentile default RF'))
print(paste0(round(100*cnt/length(data_neg$RF_IE_lwr),2),'% of RFs < 2.5th percentile default RF'))
data_neg$EQ_IE_CC <- data_neg$Conc_IE_upr/data_neg$Conc_CC
data_neg$EQ_RF_IE <- data_neg$Conc_RF_upr/data_neg$Conc_IE_upr

#Calculate ESI+ Default RF bootstrap on global data:
glob_percentiles<-mc_percentiles(data,
                                 runs=10000,
                                 samples = length(unique(data[complete.cases(data$RF),]$MS_Ready_DTXCID)),
                                 descriptors = c('MS_Ready_DTXCID','RF'))
quantile(glob_percentiles$percentile_2p5th,probs=c(0.025,0.5,0.975))
print(100*sum(data$RF < quantile(glob_percentiles$percentile_2p5th,0.5))/nrow(data))
quantile(glob_percentiles$percentile_97p5th,probs=c(0.025,0.5,0.975))
print(100*sum(data$RF > quantile(glob_percentiles$percentile_97p5th,0.5))/nrow(data))

#Calculate ESI- Default RF bootstrap on global data:
glob_percentiles_neg<-mc_percentiles(data_neg,
                                     runs=10000,
                                     samples = length(unique(data_neg[complete.cases(data_neg$RF),]$MS_Ready_DTXCID)),
                                     descriptors = c('MS_Ready_DTXCID','RF'))
quantile(glob_percentiles_neg$percentile_2p5th,probs=c(0.025,0.5,0.975))
print(100*sum(data_neg$RF < quantile(glob_percentiles_neg$percentile_2p5th,0.5))/nrow(data_neg))
quantile(glob_percentiles_neg$percentile_97p5th,probs=c(0.025,0.5,0.975))
print(100*sum(data_neg$RF > quantile(glob_percentiles_neg$percentile_97p5th,0.5))/nrow(data_neg))

#Calculate mixed models on global RF data:
mixmods_glob_ESIpos <- lmer(data=data,formula=RF_boxcox~logIE_Pred+(1|MS_Ready_DTXCID))
mixmods_glob_boot_pos <- bootMer(x=mixmods_glob_ESIpos,
                                 FUN = mySumm,
                                 type = 'parametric',
                                 use.u = FALSE,
                                 nsim = 10000,
                                 seed = 16384)
mixmods_PI_glob_pos <- sumBoot(mixmods_glob_boot_pos)
data$mix_fit <- mixmods_PI_glob_pos$fit
data$mix_lwrPI <- mixmods_PI_glob_pos$lwr
data$mix_uprPI <- mixmods_PI_glob_pos$upr
lm(data,formula=mix_lwrPI~logIE_Pred)
lm(data,formula=mix_uprPI~logIE_Pred)
lm(data,formula=mix_fit~logIE_Pred)

mixmods_glob_ESIneg <- lmer(data=data_neg,formula=RF_boxcox~logIE_Pred+(1|MS_Ready_DTXCID))
mixmods_glob_boot_neg <- bootMer(x=mixmods_glob_ESIneg,
                                 FUN = mySumm,
                                 type = 'parametric',
                                 use.u = FALSE,
                                 nsim = 10000,
                                 seed = 16384)
mixmods_PI_glob_neg <- sumBoot(mixmods_glob_boot_neg)
data_neg$mix_fit <- mixmods_PI_glob_neg$fit
data_neg$mix_lwrPI <- mixmods_PI_glob_neg$lwr
data_neg$mix_uprPI <- mixmods_PI_glob_neg$upr
lm(data_neg,formula=mix_lwrPI~logIE_Pred)
lm(data_neg,formula=mix_uprPI~logIE_Pred)
lm(data_neg,formula=mix_fit~logIE_Pred)
  #Save results:
# saveRDS(data,'RF_IE_Conc_Estimates_ESIpos.rds')
# saveRDS(data_neg,'RF_IE_Conc_Estimates_ESIneg.rds')
# saveRDS(CV_output_neg,'CV_bootstrap_ESIneg.rds')
# saveRDS(mixmods_boot_ESIneg,'mixmods_boot_ESIneg.rds')
# write.xlsx(data,'Supp_Table_ESIpos.xlsx')
# write.xlsx(data_neg,'Supp_Table_ESIneg.xlsx')

# Residual analysis for mixed-effect models:
# coeffs <- c(-58.3,201.5)
# tmpfit <- coeffs[[2]]*data$logIE_Pred+coeffs[[1]]
# tmpresid <- tmpfit-data$RF_boxcox
# 
# ggplot(data)+
#   geom_point(aes(x=logIE_Pred,
#                  y=RF_boxcox))+
#   geom_abline(aes(slope=coeffs[[2]],
#                   intercept = coeffs[[1]]))
# tmpshap <- shapiro.test(tmpresid)
# tmpresid_qq <- ggplot()+
#   stat_qq(aes(sample=tmpresid),
#           size=3,shape='O')+
#   stat_qq_line(aes(sample=tmpresid),
#                lwd=1.5,color='red')+
#   labs(y='Residual Quantiles',
#        x='Theoretical Quantiles',
#        title = paste0('ESI+ Mean CV Residual QQ-Plot\n','W-statistic = ',
#                       round(tmpshap$statistic,3),
#                       ' p-value = ',
#                       signif(tmpshap$p.value,digits = 3)))+
#   theme_classic()
# 
# ggsave(tmpresid_qq,
#        height=4,
#        width=4.5,
#        device='jpeg',
#        dpi=1200,
#        filename='Residual_QQ_Plot_MeanCV_ESI+.jpeg')

errvec <- c()
for (i in 1:length(unique(data$MS_Ready_DTXCID)))
{ 
  tmpdf <- data[data$MS_Ready_DTXCID == unique(data$MS_Ready_DTXCID)[[i]],]
  tmpdf <- tmpdf[complete.cases(tmpdf$Normalized_Intensity),]
  concmat <- matrix(nrow=nrow(tmpdf),ncol=nrow(tmpdf)-1)
  for (j in 1:(nrow(tmpdf)))
  {
    concmat[j,] <- tmpdf[seq(from=1,to=nrow(tmpdf)) != j,]$Normalized_Intensity/tmpdf$RF[j]/tmpdf[seq(from=1,to=nrow(tmpdf)) != j,]$Concentration
  }
  errvec <- append(errvec,values=as.vector(concmat))
}

RF_surr_gg<-ggplot()+
  geom_histogram(aes(x=errvec,
                     fill='Conc_RF_Surr/Conc_Known'),
                 color='black',
                 alpha=0.25)+
  labs(x='Error Quotient',
       y='Frequency')+
  scale_x_log10()+
  theme_classic()+
  theme(panel.grid = element_blank(),
        text = element_text(size=15),
        plot.title = element_text(hjust=0.5),
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=1.5))+
  geom_vline(aes(xintercept=quantile(errvec,0.975)),
             color='black',lty=2,lwd=1.5)+
  geom_vline(aes(xintercept=quantile(errvec,0.025)),
             color='black',lty=2,lwd=1.5)+
  annotate('text',label='Pred. Error ± 2.5x',x=1,y=1950)

Conc_CC_gg<-ggplot()+
  geom_histogram(aes(x=data$EQ_CCupr_CCest,
                     fill='Conc_0.975_CC/Conc_CC'),
                 color='black',
                 alpha=0.25)+
  scale_x_log10(limits=c(1,1000))+
  theme_classic()+
  labs(x='Error Quotient',
       y='Frequency')+
  theme(panel.grid = element_blank(),
        text = element_text(size=15),
        plot.title = element_text(hjust=0.5),
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=1.5))

RF_surr_comparison<-arrangeGrob(RF_surr_gg,
                                Conc_CC_gg,
                                nrow=2)
ggsave(plot=RF_surr_gg,
       dpi=1200,
       height=5,
       width=8,
       device='jpeg',
       filename='RF_Surr_ErrorQuotients.jpeg')

mixmods_logRF_pos = lmer(data=data, formula=log10(RF)~logIE_Pred+(1|MS_Ready_DTXCID))
mixmods_logRF_boot_pos <- bootMer(x=mixmods_logRF_pos,
                                  FUN = mySumm,
                                  type = 'parametric',
                                  use.u = FALSE,
                                  nsim = 10000,
                                  seed = 16384)
mixmods_PI_logRF_pos <- sumBoot(mixmods_logRF_boot_pos)
data$logmix_fit <- mixmods_PI_logRF_pos$fit
data$logmix_lwrPI <- mixmods_PI_logRF_pos$lwr
data$logmix_uprPI <- mixmods_PI_logRF_pos$upr
lm_logRF_lwrPI_smooth <- lm(data,formula=logmix_lwrPI~logIE_Pred)
lm_logRF_uprPI_smooth <- lm(data,formula=logmix_uprPI~logIE_Pred)
lm_logRF_fit_smooth <- lm(data,formula=logmix_fit~logIE_Pred)
data$logmix_lwrPI_smooth <- lm_logRF_lwrPI_smooth$coefficients[[2]]*data$logIE_Pred+lm_logRF_lwrPI_smooth$coefficients[[1]]
data$logmix_uprPI_smooth <- lm_logRF_uprPI_smooth$coefficients[[2]]*data$logIE_Pred+lm_logRF_uprPI_smooth$coefficients[[1]]
data$logmix_fit_smooth <- lm_logRF_fit_smooth$coefficients[[2]]*data$logIE_Pred+lm_logRF_fit_smooth$coefficients[[1]]
logmix_gg_pos <- ggplot(data=data)+
  geom_point(aes(x=logIE_Pred,
                 y=log10(RF)),
             size=3,
             shape='O')+
  geom_line(aes(x=logIE_Pred,
                y=logmix_lwrPI_smooth),
            lwd=1.5,lty=2)+
  geom_line(aes(x=logIE_Pred,
                y=logmix_uprPI_smooth),
            lwd=1.5,lty=2)+
  geom_line(aes(x=logIE_Pred,
                y=logmix_fit_smooth),
            lwd=1.5,color='blue')+
  theme_classic()+
  labs(y=expression(log[10]~Response~Factor),
       x=expression(log[10]~Predicted~IE),
       title = paste0(round(100*sum(log10(data$RF)<data$logmix_lwrPI_smooth)/length(data$RF),2),
                      '% Outside Lower 95% PI Bound\n',
                      round(100*sum(log10(data$RF)>data$logmix_uprPI_smooth)/length(data$RF),2),
                      '% Outside Upper 95% PI Bound'))
100-5.43
ggsave(logmix_gg_pos,
       device='jpeg',
       dpi=1200,
       height=4,
       width=5,
       filename='lg2_bootMer_95PI_logRF_ESI+.jpeg')

data$Conc_logmix_IE_upr <- as.numeric(NA)
cnt<-0
for (i in 1:nrow(data))
{
  #Taking min(tRF) as close to the LOD, 
  #If IE-calculated RF is <= min(tRF)/sqrt(2), impute min(tRF)/sqrt(2):
  min_impute_pos <- log10(min(data$RF)/sqrt(2))
  if(data$logmix_lwrPI_smooth[[i]] <= min_impute_pos)
  {
    cnt<-cnt+1
    #convert RF back to linear space and calculate concentration as I/RF:
    min_impute_pos_lin <- 10^(min_impute_pos)
    data$Conc_logmix_IE_upr[[i]] <- data$Normalized_Intensity[[i]]/min_impute_pos_lin  
  }
  else
  {
    #convert RF back to linear space and calculate concentration as I/RF:
    tmp_RF_IE <- 10^(data$logmix_lwrPI_smooth[[i]])
    data$Conc_logmix_IE_upr[[i]] <- data$Normalized_Intensity[[i]]/tmp_RF_IE
  }
}
data$Conc_logmix_IE_lwr = data$Normalized_Intensity/(10^(data$logmix_uprPI_smooth))
#Print Count and percentage of imputed RFs:
print(paste0(cnt,'/',length(data$RF),' RFs < min(tRF)/sqrt(2)'))
print(paste0(round(100*cnt/length(data$RF_IE_lwr),2),'% of RFs < min(tRF)/sqrt(2)'))
100-100*(sum(data$Conc_logmix_IE_upr<data$Concentration)+sum(data$Conc_logmix_IE_lwr>data$Concentration))/length(data$RF)
