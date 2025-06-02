#ENTACT_MonteCarlo_Percentiles.R - Louis Groff
#May 14 2021

library(ggplot2)
library(readxl)
library(xlsx)
library(gridExtra)
library(caret)

ID_list <- unique(data$New_ID)
sample_apply <- function(ID_list, data)
{
  tmpdf <- data[data$New_ID == ID_list,]$RF
  tmpdf <- tmpdf[complete.cases(tmpdf)]
  return(sample(tmpdf,1))
}
sample_out<-sapply(ID_list,FUN = function(x) sample_apply(x,data))

mc_percentiles <- function(data, descriptors, runs, samples)
{
  #store rng seeds for n runs for reproducibility
  #comment out for true random MC
  set.seed(16384)

  rf_data <- data[,descriptors]
  #remove NAs from RF list and preallocate data frames/columns:
  rf_data <- rf_data[complete.cases(rf_data$RF),]
  percentile_dist <- data.frame(percentile_5th = double(runs))
  percentile_dist$RF_sample <- c()
  percentile_dist$New_ID <- c()
  percentile_5th <- vector('list',runs)
  rf_resample <- matrix(nrow = runs, ncol = samples)
  New_ID <- matrix(nrow = runs, ncol = samples)
  unique_chems <- unique(rf_data$New_ID)
  print(paste0('Begin boostrap Simulation for N = ',runs, ' runs. '))
  for (i in 1:runs)
  {
    #get indexes of which chemicals to randomly sample:
    chem_idx <- sample(length(unique_chems),
                       size=samples,
                       replace=T)
    #resampled list of chemicals and their indices:
    chem_resample <- unique_chems[chem_idx]
    #choosing one random RF per New_ID:
    sample_out<-sapply(chem_resample, FUN = function(x) sample_apply(x,rf_data))
    New_ID[i,] <- names(sample_out)
    rf_resample[i,] <- sample_out
    percentile_5th[[i]] <- quantile(rf_resample[i,],0.05)[[1]] 
    print(paste0(round(100*i/runs,2),'% complete.'))
    percentile_dist$RF_resample[[i]] <- rf_resample[i,]
    percentile_dist$New_ID[[i]] <- New_ID[i,]
    percentile_dist$percentile_5th[[i]] <- percentile_5th[[i]]
  }
  return(percentile_dist)
}
datadir <- 'E:/EPA/Data/ENTACT Semi Quant/Standard Mixtures'
# datadir <- 'C:/Users/lgroff/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/Data/ENTACT Semi Quant/Standard Mixtures/JRS Tracer Analysis'
setwd(datadir)
data <- read_xlsx('Supp_Table1_Pos.xlsx')
data$RF <- data$Normalized_Intensity/data$Concentration
percentiles<-mc_percentiles(data,
                            runs=10000,
                            samples = length(unique(data[complete.cases(data$RF),]$New_ID)),
                            descriptors= c('New_ID','RF'))

# calculate global distribution percentiles for comparison:
rf_resample
#plot the results:
example_RF_resample<-ggplot()+
  geom_histogram(data=data,
               aes(x=RF),
               color='black',
               fill='darkgreen',
               alpha=0.25)+
  annotate('text',x=2e8,y=200,label='Global',color='darkgreen')+
  geom_histogram(data=data.frame(percentiles$RF_resample[1]),
               aes(x=RF),
               color='black',
               fill='darkblue',
               alpha=0.25)+
  annotate('text',x=3e8,y=75,label='Resample',color='darkblue')+
  theme_classic()+
  theme(axis.text.x = element_text(angle=45,hjust=1))+
  labs(x='Response Factor',
       y='Frequency')

example_RF_resample_geom<-ggplot()+
  geom_histogram(data=data,
                 aes(x=RF),
                 color='black',
                 fill='darkgreen',
                 alpha=0.25)+
  annotate('text',x=2e8,y=200,label='Global',color='darkgreen')+
  geom_histogram(data=data.frame(percentiles_geom$RF_resample[1]),
                 aes(x=RF),
                 color='black',
                 fill='darkblue',
                 alpha=0.25)+
  annotate('text',x=3e8,y=75,label='Resample',color='darkblue')+
  theme_classic()+
  theme(axis.text.x = element_text(angle=45,hjust=1))+
  labs(x='Response Factor',
       y='Frequency')
q_2p5<-quantile(percentiles$percentile_5th,0.025)
q_97p5<-quantile(percentiles$percentile_5th,0.975)
q_2p5_geom <- quantile(percentiles_geom$percentile_5th,0.025)
q_97p5_geom <- quantile(percentiles_geom$percentile_5th,0.975)
bootstrap_results_5th<-ggplot()+
  geom_histogram(data=percentiles,
               aes(x=percentile_5th),
               color='black',
               fill='darkred',
               alpha=0.25)+
  geom_vline(aes(xintercept=q_2p5),
             lwd=1.5,lty=2,color='darkred')+
  geom_vline(aes(xintercept=q_97p5),
             lwd=1.5,lty=2,color='darkred')+
  labs(x=expression(5^th~Percentile~Bootstrapped~RF),
       y='Frequency')+
  theme_classic()+
  theme(axis.text.x = element_text(angle=45, hjust=1))
bootstrap_results_5th_geom<-ggplot()+
  geom_histogram(data=percentiles_geom,
                 aes(x=percentile_5th),
                 color='black',
                 fill='purple',
                 alpha=0.25)+
  geom_vline(aes(xintercept=q_2p5_geom),
             lwd=1.5,lty=2,color='purple')+
  geom_vline(aes(xintercept=q_97p5_geom),
             lwd=1.5,lty=2,color='purple')+
  labs(x=expression(5^th~Percentile~Bootstrapped~RF),
       y='Frequency')+
  theme_classic()+
  theme(axis.text.x = element_text(angle=45, hjust=1))
grid.arrange(bootstrap_results_5th,example_RF_resample,
             bootstrap_results_5th_geom,example_RF_resample_geom,
             ncol=2,nrow=2)
# grid.arrange(plot_1st,plot_5th,ncol=2)
# percentile_plot<-arrangeGrob(plot_1st,plot_5th,ncol=2)
# ggsave(percentile_plot,file='MC_Percentile_Dist_Bootstrap_N1996_10kRuns.jpeg',
#        device='jpeg',dpi=1200)

#Calculating percentage of RFs outside median value of 5th percentile:
#Global RF Distribution:
glob_5th <- quantile(data$RF,na.rm=T,0.05)
glob_mean <- mean(na.rm=T,data$RF)
glob_median <- median(na.rm=T,data$RF)
glob_percent_out <- nrow(data[complete.cases(data$RF) & (data$RF < glob_5th),'RF'])/nrow(data[complete.cases(data$RF),'RF'])
#Bootstrap 1 RF Per Chemical:
boot_mean <- mean(percentiles$percentile_5th)
boot_median <- median(percentiles$percentile_5th)
boot_2p5 <- quantile(percentiles$percentile_5th,0.025)
boot_97p5 <- quantile(percentiles$percentile_5th,0.975)
boot_percent_out <- nrow(data[complete.cases(data$RF) & (data$RF < boot_median),'RF'])/nrow(data[complete.cases(data$RF),'RF'])
#Bootstrap Geometric Mean:
boot_mean_geom <- mean(percentiles_geom$percentile_5th)
boot_median_geom <- median(percentiles_geom$percentile_5th)
boot_2p5_geom <- quantile(percentiles_geom$percentile_5th,0.025)
boot_97p5_geom <- quantile(percentiles_geom$percentile_5th,0.975)
boot_percent_out_geom <- nrow(data[complete.cases(data$RF) & (data$RF < boot_median_geom),'RF'])/nrow(data[complete.cases(data$RF),'RF'])

data=data
runs=10
samples=551
folds=5
percentile_bootstrap_kfoldCV <- function(data, descriptors, folds, runs, samples)
{
  set.seed(16384)
  rf_data <- data[,descriptors]
  #remove NAs from RF list and preallocate data frames/columns:
  rf_data <- rf_data[complete.cases(rf_data$RF),]
  rf_folds<- createFolds(rf_data$RF,
                         k=folds,
                         list=T)
  fold_idx<-seq(1:folds)
  mc_out <- vector(mode='list',length=folds)
  train_folds <- mc_out
  test_folds <- mc_out
  for(i in 1:folds){
    rf_test <- rf_data[rf_folds[[i]],]
    rf_train <- data.frame(RF = double(nrow(rf_data)-nrow(rf_test)),
                           New_ID = character(nrow(rf_data)-nrow(rf_test)))
    rf_train$RF <- rf_data[unlist(rf_folds[fold_idx != i]),]$RF
    rf_train$New_ID <- rf_data[unlist(rf_folds[fold_idx != i]),]$New_ID
    train_percentiles<-mc_percentiles(data=rf_train,
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
#Structure of output:
#CV_output[[1]] is the training set RF data for all CV iterations
#(e.g. CV_output[[1]][[3]] is the training set for the 3rd CV iteration)
#CV_output[[2]] is the test set RF data for all CV iterations 
#(e.g. CV_output[[2]][[3]] is the test set for the 3rd CV iteration)
#CV_output[[3]] is the bootstrap output for all CV iterations
#(e.g, CV_output[[3]][[3]] is the bootstrap output for the 3rd CV iteration)

CV_output<-percentile_bootstrap_kfoldCV(data=data,
                                        descriptors=c('New_ID','RF'),
                                        folds=5,
                                        runs=10000,
                                        samples=551)

#Process results from CV Bootstrap RF Model:
train_test_ggs <- vector(mode='list',length=5)
for(i in 1:5)
{
  train_test_ggs[[i]]<-ggplot()+
    geom_histogram(data=CV_output[[1]][[i]],
                   aes(x=RF,fill='Training'),
                   color='black',
                   alpha=0.25)+
    geom_histogram(data=CV_output[[2]][[i]],
                   aes(x=RF,fill='Test'),
                   color='black',
                   alpha=0.25)+
    
    theme_classic()+
    labs(x='Response Factor',
         y='Frequency',
          title=paste0('CV Iteration ',i))+
    theme(axis.text.x = element_text(angle=45,hjust=1))
}

grid.arrange(train_test_ggs[[1]],
             train_test_ggs[[2]],
             train_test_ggs[[3]],
             train_test_ggs[[4]],
             train_test_ggs[[5]],
             nrow=2,ncol=3)
# train_test_gg_fig<-arrangeGrob(train_test_ggs[[1]],
#                                train_test_ggs[[2]],
#                                train_test_ggs[[3]],
#                                train_test_ggs[[4]],
#                                train_test_ggs[[5]],
#                                nrow=2,ncol=3)
# ggsave(train_test_gg_fig,
#        filename = 'CV_TrainTestSets.jpeg',
#        dpi=1200,
#        device='jpeg',
#        height=8,
#        width=18)

percentile_dist_ggs <- vector(mode='list',length=5)
q2p5_cv <- vector(mode='list',length=5)
q97p5_cv <- q2p5_cv
median_cv <- q2p5_cv
median_out_percent <- q2p5_cv
for(i in 1:5)
{
  q2p5_cv[[i]] <- quantile(CV_output[[3]][[i]]$percentile_5th,0.025)
  q97p5_cv[[i]] <- quantile(CV_output[[3]][[i]]$percentile_5th,0.975)
  median_cv[[i]] <- median(CV_output[[3]][[i]]$percentile_5th)
  median_out_percent[[i]] <- length(CV_output[[2]][[i]][CV_output[[2]][[i]]$RF < median_cv[[i]],]$RF)/length(CV_output[[2]][[i]]$RF)
  percentile_dist_ggs[[i]]<-ggplot()+
    geom_histogram(data=CV_output[[3]][[i]],
                   aes(x=percentile_5th),
                   color='black',
                   fill='darkred',
                   alpha=0.25)+
    geom_vline(aes(xintercept=q2p5_cv[[i]]),
               lwd=1.5,lty=2)+
    geom_vline(aes(xintercept=q97p5_cv[[i]]),
               lwd=1.5,lty=2)+
    theme_classic()+
    labs(x='Bootstrapped 5th Percentile\nTraining Set RF',
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
#                                     nrow=2,ncol=3)
# ggsave(percentile_dist_gg_fig,
#        filename = 'CV_PercentileDistributions.jpeg',
#        dpi=1200,
#        device='jpeg',
#        height=8,
#        width=18)
# ggplot()

