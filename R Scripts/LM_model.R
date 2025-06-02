library(readr)
library(readxl)
library(caret)
library(rsample)
library(MASS)
library(dplyr)
library(gridExtra)

#SEED!
set.seed(1313)

#Functions
get_model_lm <- function(data, endpoint, descriptors) 
{
  # fitControl <- trainControl(method = "cv", 
  #                            number = 5, 
  #                            savePredictions = 'final',
  #                            allowParallel = TRUE)
  
  data_model=dplyr::select(data, tidyselect::any_of(descriptors))
  
  fit_model <- train(x=data_model,
                       y= endpoint,
                       #preProcess = c("center", "scale"),
                       method = "lm")
  
  #metric=min(fit_model$results["RMSE"])
  #metric=fit_model$results
  
  return(fit_model)
}
datadir2<-'C:/Users/lgroff/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/Data/ENTACT Semi Quant/Standard Mixtures/Anneli Updates'
setwd(datadir2)
# setwd("//aa.ad.epa.gov/ORD/RTP/USERS/A-D/CLowe/Net MyDocuments/R Data/Bayes_SQ")

df <- read_excel(paste0(datadir2,"/Supplemental_Tables_v3.xlsx"), skip=1, 
                 sheet = "Table S3", col_types = c("text", 
                                                   "text", "text", "text", "text", "text", 
                                                   "text", "text", "text", "text", "numeric", 
                                                   "numeric", "numeric", "numeric", 
                                                   "numeric", "numeric", "numeric", 
                                                   "numeric", "numeric", "numeric", 
                                                   "numeric", "numeric", "numeric", 
                                                   "numeric", "numeric", "numeric", 
                                                   "numeric", "numeric"))

#Let's get rid of most of the columns...
df <- df[,c("DTXSID with Isomer Order","Transformed Response Factor (RF(1/6))","log10 Predicted IE")]
colnames(df) <- c('DTXSID_wISO','tRF_6throot','logIE_Pred')
#NAs are pesky and must be removed.
df <- df[which(is.na(df$logIE_Pred) == F),]

#This is the magic data splitting function that maintains the distribution 
#of your endpoint in both training and test sets
data_split <- initial_split(df,
                            strata = "tRF_6throot",
                            prop = 0.75)

#Make training and test sets
data_train <- training(data_split) %>% as.data.frame()
data_test  <- testing(data_split) %>% as.data.frame()

#Graphics showing endpoint disbtributions
g_train <- ggplot(data = data_train,aes(x=(tRF_6throot))) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
  ggtitle("Training Set") +
  theme_classic(base_size = 16) +
  labs(x =expression(Response~Factor^(1/6)))
# g_train
# ggsave(g_train,file='TrainingSet_DensityHistogram.jpeg',
#        device='jpeg',dpi=1200)
g_test <- ggplot(data = data_test,aes(x=(tRF_6throot))) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
  ggtitle("Test Set") +
  theme_classic(base_size = 16) + 
  labs(x = expression(Response~Factor^(1/6)))
# g_test
# ggsave(g_test,file='ExternalTestSet_DensityHistogram.jpeg',
#        device='jpeg',dpi=1200)
grid.arrange(g_train,g_test,ncol=2)
#group histograms and ggsave:
train_test_density<-arrangeGrob(g_train,g_test,ncol=2)
# ggsave(plot=train_test_density,file='Train_Test_DensityPlots.jpeg',
#        device='jpeg',dpi=1200)
#Every good (non decision tree) model scales the data first..
preProc <- preProcess(data_train)

training <- predict(preProc, data_train)
testing <- predict(preProc, data_test)
#training <- data_train
#testing <- data_test

#Build model
model <- get_model_lm(training,
                       training$tRF_6throot,
                       c("logIE_Pred"))

#Predictions from model
training$RF_Predicted <- predict(model)
#Q-Q style plot
pred_vs_known_training <- ggplot()+
  geom_point(aes(y=training$RF_Predicted,
                 x=training$tRF_6throot),
             color='blue',
             size=3,
             shape='O')+
  geom_abline(slope=1,
              intercept=0,
              lwd=1.5)+
  labs(x=expression(Actual~Response~Factor^(1/6)),
       y=expression(Predicted~Response~Factor^(1/6)),
       title='ESI+ Training Set')+
  theme_classic()
# ggsave(pred_vs_known_training,file='tRF_pred_vs_tRF_known.jpeg',
#        device='jpeg',dpi=1200, height=4,width = 5)
pred_vs_known_training
#Generate performance metrics
rlm_predictions <- predict(model$finalModel,newdata = testing[c("tRF_6throot","logIE_Pred")])
model
model_testing<-postResample(pred = rlm_predictions, obs = testing$tRF_6throot)
model_testing
#Plot regression on training data
plot_regress_training<-ggplot()+
  geom_point(aes(x=training$logIE_Pred,
                 y=training$tRF_6throot),
             color='blue',
             size=3,
             shape='O')+
  geom_abline(intercept = model$finalModel$coefficients[[1]][1],
              slope = model$finalModel$coefficients[[2]][1],
              lwd=1.5)+
  theme_classic()+
  labs(x=expression(log[10]~Predicted~IE),
       y=expression(Response~Factor^(1/6)),
       title='ESI+ Training Set')+
  annotate('text', 
           x=1.75,y=-2.25,
           label=paste0('R^2 = ',round(model[[4]][[3]],4)))
plot_regress_training
#Plot regression on testing data
plot_regress_testing<-ggplot()+
  geom_point(aes(x=testing$logIE_Pred,
                 y=testing$tRF_6throot),
             color='blue',
             size=3,
             shape='O')+
  geom_abline(intercept = model$finalModel$coefficients[[1]][1],
              slope = model$finalModel$coefficients[[2]][1],
              lwd=1.5)+
  theme_classic()+
  labs(x=expression(log[10]~Predicted~IE),
       y=expression(Response~Factor^(1/6)),
       title='ESI+ External Test Set')+
  annotate('text', 
           x=1.5,y=-2.25,
           label=paste0('R^2 = ',round(model_testing[[2]],4)))
plot_regress_testing
grid.arrange(plot_regress_training,plot_regress_testing,ncol=2)
# for saving the grid image:
# pos_5cv_plot<-arrangeGrob(plot_regress_training,plot_regress_testing,ncol=2)
# ggsave(pos_5cv_plot,file='pos_5cv_plot.jpeg',device='jpeg',dpi=1200)
