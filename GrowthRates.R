#Calculate growth rates and carrying capacity



#Use growthrates package to calculate r and K

#Step 1: fit each growth curve individually
#Clear environment
rm(list=ls())
library(growthrates)
mydata<-read.csv("Data/CellCountsmax.csv")

mydata$count<-mydata$X16.1631.30.3 #change genotype each time
mydata<-mydata[-14,] #subtract off rows if no data
mydata<-mydata[-13,]
mydata<-mydata[-12,]
mydata<-mydata[-11,]

View(mydata)

#Set initial parameters
p1 <- c(y0=0.01, mumax = 0.2, K=50)
lower1<-c(y0=1e-6, mumax =0, K=10)
upper1<-c(y0=0.05, mumax=5, K=300)


logfit<-fit_growthmodel(FUN=grow_logistic, p=p1, time=mydata$time, y=mydata$count, lower=lower1, upper=upper1)
coef(logfit)
rsquared(logfit)

plot(logfit)



#Now do analyses with data from growthrates package
#Import the new data file with r and K for each genotype
#Clear environment
rm(list=ls())


mydata2<-read.csv("Data/GrowthRates.csv")
mydata2$Temp<-as.factor(mydata2$Temp)
mydata2$Genotype<-as.factor(mydata2$Genotype)


model1<-lm(newr~Temp*Genotype, data=mydata2)
anova(model1)
plot(model1) #Check assumptions
library(car)
qqp(residuals(model1), "norm")

model2<-lm(newK~Temp*Genotype, data=mydata2)
anova(model2)
plot(model2)
qqp(residuals(model2), "norm")


#Normality is meh. Variances are pretty unequal. Try log-transform
mydata2$logr<-log(mydata2$newr)
mydata2$logK<-log(mydata2$newK)

model3<-lm(logr~Temp*Genotype, data=mydata2)
model4<-lm(logK~Temp*Genotype, data=mydata2)
plot(model3)
plot(model4)
library(car)
qqp(residuals(model3), "norm")
qqp(residuals(model4), "norm")

#Homogeneity of variances better with log data. Still pretty normal.
#Use log-trans data for results. Going to plot raw data though.

#Get Results
anova(model3)
anova(model4)

#Post-hoc tests
model3
library(emmeans)
rmeans<-emmeans(model3, pairwise~Temp*Genotype, adjust="tukey")
rmeans

Kmeans<-emmeans(model4, pairwise~Temp*Genotype, adjust="tukey")
Kmeans

#Tukey tests
tx<-with(mydata2, interaction(Temp, Genotype))
rmod<-lm(logr~tx, data=mydata2)
library(agricolae)
HSD.test(rmod, "tx", console=TRUE)

Kmod<-lm(logK~tx, data=mydata2)
HSD.test(Kmod, "tx", console=TRUE)

#Are r and K correlated with each other?
plot(logr~logK, data=mydata2)
correlation<-cor.test(mydata2$logr, mydata2$logK, method="pearson")
correlation

#Is there a difference in time to reach K?
model5<-lm(TimetoK~Temp*Genotype, data=mydata2)
anova(model5)


#graphs
library(tidyr)
library(dplyr)
graphdata2 <- mydata2 %>%
  group_by(Temp, Gs) %>%
  summarize(meanr=mean(newr), meanK=mean(newK), logmeanr=mean(logr), logmeanK=mean(logK), 
            ser=sd(newr)/sqrt(length((newr))), seK=sd(newK)/sqrt(length((newK))), selogr=sd(logr)/sqrt(length((logr))), selogK=sd(logK)/sqrt(length((logK))))
graphdata2

###BarPlot of growth rates
library(ggplot2)
ggplot(graphdata2, aes(x=Gs, y=meanr, fill=factor(Temp), group=factor(Temp))) + #basic plot with TidalHeight as a grouping factor
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=14), axis.title = element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=14))+
  geom_bar(stat="identity", position="dodge", size=0.6) + #determines the bar width
  geom_errorbar(aes(ymax=meanr+ser, ymin=meanr-ser), stat="identity", position=position_dodge(width=0.9), width=0.1) + #adds error bars
  labs(x="Genotype", y="Growth Rate (r)", fill="Temperature") + #labels the x and y axes
  scale_fill_manual(values=c("26"="dodgerblue2","30"="red3")) #fill colors for the bars
ggsave("rGenotypesTemp barplot.png", dpi=800, height=5, width=7)


###BarPlot of carrying capacity
library(ggplot2)
ggplot(graphdata2, aes(x=Gs, y=meanK, fill=factor(Temp), group=factor(Temp))) + #basic plot with TidalHeight as a grouping factor
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=14), axis.title = element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=14))+
  geom_bar(stat="identity", position="dodge", size=0.6) + #determines the bar width
  geom_errorbar(aes(ymax=meanK+seK, ymin=meanK-seK), stat="identity", position=position_dodge(width=0.9), width=0.1) + #adds error bars
  labs(x="Genotype", y="Carrying Capacity (K)", fill="Temperature") + #labels the x and y axes
  scale_fill_manual(values=c("26"="dodgerblue2","30"="red3")) #fill colors for the bars
ggsave("KGenotypesTemp barplot.png", dpi=800, height=5, width=7)


###BoxPlot of log growth rates
library(ggplot2)
ggplot(mydata2, aes(x=Gs, y=logr, fill=Temp)) +
  geom_boxplot(position=position_dodge(1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Genotype", y="log growth rate (r)", fill="Temperature") +
  scale_fill_manual(values=c("26"="dodgerblue2","30"="red3"))

ggsave("rGenotypesTemp boxplot.png", dpi=300, height=5, width=7)

###BoxPlot of log K
library(ggplot2)
ggplot(mydata2, aes(x=Gs, y=logK, fill=Temp)) +
  geom_boxplot(position=position_dodge(1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Genotype", y="log carrying capacity (r)", fill="Temperature") +
  scale_fill_manual(values=c("26"="dodgerblue2","30"="red3"))

ggsave("KGenotypesTemp boxplot.png", dpi=300, height=5, width=7)


###BoxPlot of r 
library(ggplot2)
ggplot(mydata2, aes(x=Gs, y=newr, fill=Temp)) +
  geom_boxplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Genotype", y="Growth Rate (r)", fill="Temperature") +
  scale_fill_manual(values=c("26"="dodgerblue2","30"="red3"))

ggsave("rGenotypesTemp boxplot.png", dpi=600, height=5, width=7)

###BoxPlot of log K
library(ggplot2)
ggplot(mydata2, aes(x=Gs, y=newK, fill=Temp)) +
  geom_boxplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Genotype", y="Carrying Capacity (K)", fill="Temperature") +
  scale_fill_manual(values=c("26"="dodgerblue2","30"="red3"))

ggsave("KGenotypesTemp boxplot.png", dpi=600, height=5, width=7)


#Time to Reach K
library(ggplot2)
ggplot(mydata2, aes(x=Gs, y=TimetoK, fill=Temp)) +
  geom_boxplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Genotype", y="Time to Reach K", fill="Temperature") +
  scale_fill_manual(values=c("26"="dodgerblue2","30"="red3"))


#Try to make lattice of all fits
rm(list=ls())

#This gets fits for every replicate at every temp/genotype combo
mydata<-read.csv("Data/CellCountsforlattice.csv")
mydata$Temp<-as.factor(mydata$Temp)
mydata$Replicate<-as.factor(mydata$Replicate)

library(growthrates)
library(lattice)

xyplot(count~time|Genotype+Temp, data=mydata, groups = Replicate, pch=16, cex=0.5)

p1 <- c(y0=0.01, mumax = 0.2, K=50)
lower1<-c(y0=1e-6, mumax =0, K=10)
upper1<-c(y0=0.05, mumax=5, K=300)

many_fits<-all_growthmodels(count~grow_logistic(time, parms) |Temp + Replicate + Genotype, data=mydata,
                            p=p1, lower=lower1, upper=upper1, which = c("y0", "mumax", "K"))



png("/Users/cpt48249/Box Sync/Casey Documents/Projects-Current/Symbio Growth Rates/RGrowth/AllFits.png", width = 3600, height = 4800)
par(mfrow=c(6,5))
par(mar= c(1, 1, 1, 1))
par(pty="s")
par(mai = c(0.1,0.1,0.1,0.1))

many_fits<-all_growthmodels(count~grow_logistic(time, parms) |Temp + Replicate + Genotype, data=mydata,
                            p=p1, lower=lower1, upper=upper1, which = c("y0", "mumax", "K"))
plot(many_fits)
dev.off()


#Make line graphs with means of replicates
rm(list=ls())

#Import file that already has predicted values
mydata<-read.csv("Data/MeanCellCounts.csv")
mydata$Temp<-as.factor(mydata$Temp)
mydata$Genotype<-as.factor(mydata$Genotype)

library(growthrates)
library(lattice)

#Subset by genotype
library(tidyr)
library(dplyr)
mydata2<- mydata %>%
  filter(Genotype == "G5")
View(mydata2)

library(ggplot2)
ggplot(data=mydata2, aes(time, count, color=Temp)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=14), axis.title = element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=14))+
  xlim(0,40)+
  ylim(0,210)+
  geom_point() + 
  geom_line(aes(y=predictions)) +
  scale_color_manual(values=c("26"="dodgerblue2","30"="red3"))+
  labs(x="Time", y="Cell Density (cells/mL)", fill="Temperature")

ggsave("G5.png", dpi=800, height=5, width=7)





