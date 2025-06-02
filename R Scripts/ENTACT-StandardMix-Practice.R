library(ggplot2)
library(readxl)
library(investr)
setwd('C:/Users/lgroff/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/Data/ENTACT Semi Quant/Standard Mixtures')
posdata <- read_xlsx(path='All_PosData_for_Evaluation.xlsx')
negdata <- read_xlsx(path='All_NegData_for_Evaluation.xlsx')
posdata <- posdata[complete.cases(posdata$mean_Intensity),] #remove rows with NaN intensities
negdata <- negdata[complete.cases(negdata$mean_Intensity),]
#Plot posmode data:
ggplot(posdata,aes(x=Spike_Conc,y=mean_Intensity))+
  geom_point(aes(color=Level,size=1.5),alpha=0.5)+
  scale_x_log10()+
  scale_y_log10()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(color='black'))+
  guides(size=FALSE)+
  xlab('Standard Concentration (mM)')+
  ylab('Mean Intensity')+
  ggtitle('ENTACT Standard Mixtures, ESI+')+
  geom_smooth(method='lm', lwd=2, se = FALSE, formula = y~x, linetype='dotdash', aes(color='Regression'))+
  scale_color_manual(name='Legend',values=c('chartreuse3','blue2','darkorange','black'))

DTXUnique <- unique(posdata$DTXSID)
plts <- replicate(length(DTXUnique),data.frame())
mods <- replicate(length(DTXUnique),data.frame())
pltdata <- posdata[c('mean_Intensity','Spike_Conc','DTXSID')]
pltdata$slope <- NaN
pltdata$intercept <- NaN
for (i in 1:length(DTXUnique)){
  tmpx = c()
  tmpy = c()
  tmpdata <- pltdata[pltdata$DTXSID==DTXUnique[i],]
  mods[[i]] <- lm(log10(tmpdata$mean_Intensity) ~ log10(tmpdata$Spike_Conc))
  plts[[i]] <- ggplot(tmpdata,aes(x=Spike_Conc,y=mean_Intensity))+
    geom_point(color='blue2', aes(size = 1.5))+
    guides(size=FALSE)+
    scale_x_log10()+
    scale_y_log10()+
    xlab('Log10(Standard Concentration (mM)')+
    ylab('Log10(Mean Intensity)')+
    ggtitle(paste0(DTXUnique[i],', ESI+\n log(y) = ',round(mods[[i]]$coefficients[2],2),'*log(x) + ',round(mods[[i]]$coefficients[1],2)))+
    geom_smooth(method='lm',lwd=2, formula = y ~ x, linetype='dashed',color='black')
  pltdata$slope[pltdata$DTXSID==DTXUnique[i]] <- mods[[i]]$coefficients[2]
  pltdata$intercept[pltdata$DTXSID==DTXUnique[i]] <- mods[[i]]$coefficients[1]
  ggsave(filename=paste0(DTXUnique[i],'.png'),plot=plts[[i]])
}
mod <- lm(measured ~ actual, data = arsenic)
(res <- calibrate(mod, y0 = 3, interval = "inversion", level = 0.95))
plotFit(mod, interval = "confidence", level = 0.9, shade = TRUE, col.conf = "skyblue")
abline(h = 3, v = c(res$lower, res$estimate, res$upper), lty = 2)
