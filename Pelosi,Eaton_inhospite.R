### Analyses for Pelosi, J.A, K.M.Eaton, S.Mychajliw, C.P. terHorst, and M.A. Coffroth 
### Thermally tolerant symbionts may explain Caribbean octocoral resilience to heat stress
### Last updated May 5, 2021

library(ggplot2)
library(nlme)
library(dplyr)
library(agricolae)

### Symbiont Growth Rates -- see Casey's script (GrowthRates.R) ###########################################

## 1. Import data

symbiont_growth <- read.csv("NEW_GROWTH_RATE_JUL7.csv")

## 2. Set temperatures as factors and remove rows with N/A

symbiont_growth$Treat_temp <- as.factor(symbiont_growth$Treat_temp)
symbiont_growth$Hist_temp <- as.factor(symbiont_growth$Hist_temp)
symbiont_growth <- na.omit(symbiont_growth)

## 3. Plot the data and save graph 

ggplot(data = symbiont_growth, mapping = aes(x = Replicate.., y = mu, fill = Treat_temp)) + geom_boxplot() +
  scale_fill_manual(values = c("dodgerblue2", "red3")) + theme_classic() + ylab(expression("Growth Rate ("*mu*d^-1*")")) +
  xlab("Culture") + theme(legend.position = "none") + theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14))

ggsave("in_vitro_growth_oct15.png", dpi = 300, height = 5, width = 7)

## 4. Run ANOVA

growth_fit1 <- aov(mu ~ Treat_temp*Replicate.., data = symbiont_growth)
summary(growth_fit1)

## 5. Check assumptions

# Check homogeneity of variances
resid_1<-symbiont_growth$mu - predict.lm(growth_fit1)
plot(predict.lm(growth_fit1), resid_1)
abline(a=0, b= 0, col = "red")
# Homogeneity of variances look fine 

# Normality of residuals
stdRes <-rstandard(growth_fit1)
#qqplot
qqnorm(stdRes,ylab="Standardized Residuals", xlab="Theoretical Quantiles")
qqline(stdRes, col=2,lwd=2)
#histogram of residuals
hist(stdRes)
# Normality looks fine

## 6. Run Tukey HSD for pairwise comparisons and multiple corrections

TukeyHSD(growth_fit1)

## 7. Run LME model to determine if historical temperature has significant effect 
symbiont_growth_lme <-lme(mu~Hist_temp*Treat_temp, random= ~1|Replicate.., data=symbiont_growth, method = "REML") 
summary(symbiont_growth_lme)
anova.lme(symbiont_growth_lme, type="sequential", adjustSigma = FALSE)

# Check assumptions 
hist(residuals(symbiont_growth_lme))
plot(fitted(symbiont_growth_lme), residuals(symbiont_growth_lme))

## 8. Run LM model for only fixed effects 
symbiont_growth_no_culture <- gls(mu~Hist_temp*Treat_temp, data = symbiont_growth)
summary(symbiont_growth_no_culture)

# Check assumptions 
hist(residuals(symbiont_growth_no_culture))
plot(fitted(symbiont_growth_no_culture), residuals(symbiont_growth_no_culture))

## 9. Test significance of random effect  
anova(symbiont_growth_lme, symbiont_growth_no_culture)

### Polyp-Symbiont Cell Density #######################################################################

## Data from scored during experiment (by eye)

## 1. Import data
infection_scored <- read.csv("infection_scores_11feb2021.csv")
infection_scored <- infection_scored %>% 
  filter(Percent_infected != 'NA') %>% 
  mutate(log_Percent_infected=log(1+Percent_infected))


## 2. Set variables as factors 

infection_scored$Container <- as.factor(infection_scored$Container)
infection_scored$Tempt <- as.factor(infection_scored$Tempt)
infection_scored$Genotype <- as.factor(infection_scored$Genotype)
infection_scored$Period <- as.factor(infection_scored$Period)

infection_1 <- infection_scored %>% 
  filter(Period == 1)

infection_2 <- infection_scored %>% 
  filter(Period == 2)

infection_3 <- infection_scored %>% 
  filter(Period == 3)

## 3. Run ANOVA

ggplot(data = infection_1, mapping = aes(x = Genotype, y = Percent_infected, fill=Tempt)) + geom_boxplot() + ylab("% Visually Infected Polyps") +
  theme_classic() + xlab("Genotype") + theme(legend.position = "None", axis.text = element_text(size=14), axis.title = element_text(size=16)) +
  scale_fill_manual(values = c("dodgerblue2", "red3"))

ggsave("Fig3a_feb14_21.pdf", height = 8, width = 8)

ggplot(data = infection_2, mapping = aes(x = Genotype, y = Percent_infected, fill=Tempt)) + geom_boxplot() + ylab("% Visually Infected Polyps") +
  theme_classic() + xlab("Genotype") + theme(legend.position = "None", axis.text = element_text(size=14), axis.title = element_text(size=16)) +
  scale_fill_manual(values = c("dodgerblue2", "red3"))

ggsave("Fig3b_feb14_21.pdf", height=8, width =8)

infection_fit <- aov(Percent_infected~Genotype*Tempt*Period, data = infection_scored)
summary(infection_fit)

period_1 <- aov(Percent_infected~Genotype*Tempt, data = infection_1)
summary(period_1)

#Check assumptions

#homogeneity of variances
infection_resids <- infection_1$Percent_infected-predict.lm(period_1)
plot(predict.lm(period_1), infection_resids)
abline(a=0, b= 0, col = "red") #sort of cone-like 

# Normality of residuals
stdRes_infect <-rstandard(period_1)
#qqplot
qqnorm(stdRes_infect,ylab="Standardized Residuals", xlab="Theoretical Quantiles")
qqline(stdRes_infect, col=2,lwd=2)
#histogram of residuals
hist(stdRes_infect)
# close to normal 
summary(infection_fit)
TukeyHSD(infection_fit)

period_2 <- aov(Percent_infected~Genotype*Tempt, data = infection_2)
summary(period_2)
TukeyHSD(period_2)

#Check assumptions

#homogeneity of variances
infection_resids <- infection_2$Percent_infected-predict.lm(period_2)
plot(predict.lm(period_2), infection_resids)
abline(a=0, b= 0, col = "red") #sort of cone-like 

# Normality of residuals
stdRes_infect <-rstandard(period_2)
#qqplot
qqnorm(stdRes_infect,ylab="Standardized Residuals", xlab="Theoretical Quantiles")
qqline(stdRes_infect, col=2,lwd=2)
#histogram of residuals
hist(stdRes_infect)
# close to normal 

period_3 <- aov(Percent_infected~Genotype*Tempt, data = infection_3)
summary(period_3)

#Check assumptions

#homogeneity of variances
infection_resids <- infection_3$Percent_infected-predict.lm(period_3)
plot(predict.lm(period_3), infection_resids)
abline(a=0, b= 0, col = "red") #sort of cone-like 

# Normality of residuals
stdRes_infect <-rstandard(period_3)
#qqplot
qqnorm(stdRes_infect,ylab="Standardized Residuals", xlab="Theoretical Quantiles")
qqline(stdRes_infect, col=2,lwd=2)
#histogram of residuals
hist(stdRes_infect)
# close to normal

## 4. Check assumptions

#homogeneity of variances
infection_resids <- infection_scored$Percent_infected-predict.lm(infection_fit)
plot(predict.lm(infection_fit), infection_resids)
abline(a=0, b= 0, col = "red") #sort of cone-like 

# Normality of residuals
stdRes_infect <-rstandard(infection_fit)
#qqplot
qqnorm(stdRes_infect,ylab="Standardized Residuals", xlab="Theoretical Quantiles")
qqline(stdRes_infect, col=2,lwd=2)
#histogram of residuals
hist(stdRes_infect)
# close to normal 
summary(infection_fit)
TukeyHSD(infection_fit)

## Cell density in polyps ##############################################################

## 1. Import data 
infection <- read.csv("BCO_DMO_Cell cnts in polyps_updated.csv")

## REMOVE G6!

infection <- filter(infection, Culture != "G6")

## 2. Set variables as factors 
infection$Container <- as.factor(infection$Container)
infection$Culture <- factor(infection$Culture)
infection$Tempt <- factor(infection$Tempt)

## 3. Plot data and save graph 
ggplot(data = infection, mapping = aes(x = Culture, y = Cells.polyp, fill = Tempt)) + geom_boxplot() +
  scale_fill_manual(values = c("dodgerblue2", "red3")) + ylab("Mean # cells/polyp") + theme_classic() + xlab("Genotype") +
  theme(legend.position = "None", axis.text = element_text(size = 14), axis.title = element_text(size =16))

ggsave("polyp_cell_density_feb14_21.pdf", height = 8, width =8)

## 4. Run ANOVA 
infection_fit <- aov(Cells.polyp~Culture*Tempt, data = infection)
summary(infection_fit)

## 5. Check assumptions

#homogeneity of variances
infection_resids <- infection$Cells.polyp-predict.lm(infection_fit)
plot(predict.lm(infection_fit), infection_resids)
abline(a=0, b= 0, col = "red")
# homogeneity of var is fine 

# Normality of residuals
stdRes_infect <-rstandard(infection_fit)
#qqplot
qqnorm(stdRes_infect,ylab="Standardized Residuals", xlab="Theoretical Quantiles")
qqline(stdRes_infect, col=2,lwd=2)
#histogram of residuals
hist(stdRes_2)
# Not normal! Use log to transform with log()

# transform data with log
infect_sqrt_fit <- aov(sqrt(Cells.polyp)~Tempt*Culture, data = infection)
summary(infect_sqrt_fit)

# Normality of residuals
stdRes_infect_sqrt <-rstandard(infect_sqrt_fit)
#qqplot
qqnorm(stdRes_infect_sqrt,ylab="Standardized Residuals", xlab="Theoretical Quantiles")
qqline(stdRes_infect_sqrt, col=2,lwd=2)
#histogram of residuals
hist(stdRes_infect_sqrt)
# Better normality with sqrt transform 

## 6. Run Tukey HSD for pairwise comparisons and multiple corrections

TukeyHSD(infect_sqrt_fit)

infect_tukey <- TukeyHSD(infect_sqrt_fit)
infect_tukey

infect_HSD <- HSD.test(infect_sqrt_fit, trt = c("Culture","Tempt"), console = TRUE)

infect_groups <- infect_HSD %>% 
  .$groups %>% 
  as_tibble(rownames="class")

infect_groups <- infect_groups[order(infect_groups$class),]

ggplot()+ geom_boxplot(data = infection, mapping = aes(x = Culture, y = Cells.polyp, fill = Tempt)) +
  xlab("Culture") +ylab("Mean No. of Cells/Polyp") + scale_fill_manual(values = c("dodgerblue2", "red3")) + theme_classic() + theme(legend.position = "none") +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
  geom_text(data = infect_groups, aes(label = groups, y= 37500, x=c("G1", "G1", "G2", "G2", "G3", "G3", "G4", "G4", "G5", "G5")), position = position_dodge2(width = 0.7), size = 6)

ggsave("Fig3c_feb14-21.pdf", width = 8, height = 6)

## 7. Run LME model to determine if historical temperature has significant effect
# Add Hist_temp as column 
infection$Hist_temp[infection$Culture=="G1"] <- 26
infection$Hist_temp[infection$Culture=="G2"] <- 26
infection$Hist_temp[infection$Culture=="G3"] <- 26
infection$Hist_temp[infection$Culture=="G4"] <- 30
infection$Hist_temp[infection$Culture=="G5"] <- 30

infection_lme <-lme(sqrt(Cells.polyp)~Hist_temp*Tempt, random= ~1|Culture, data=infection, method = "REML") 
summary(infection_lme)
anova.lme(infection_lme, type="sequential", adjustSigma = FALSE)

# Check assumptions 
hist(residuals(infection_lme))
plot(fitted(infection_lme), residuals(infection_lme))

## 8. Run LM model for only fixed effects 
infection_no_culture <- gls(sqrt(Cells.polyp)~Hist_temp*Tempt, data = infection)
summary(infection_no_culture)

# Check assumptions 
hist(residuals(infection_no_culture))
plot(fitted(infection_no_culture), residuals(infection_no_culture))

## 9. Test significance of random effect  
anova(infection_lme, infection_no_culture)

### Polyp survivorship ##################################################################################

## 1. Import data 
unweighted <- read.csv("../../polyp_analyses/unweighted average polyp survivorship_2-10-21.csv")

## 2. Add column for historical temp and remove extra columns

unweighted$Hist_temp[unweighted$Treatment=="G1"] <- 26
unweighted$Hist_temp[unweighted$Treatment=="G2"] <- 26
unweighted$Hist_temp[unweighted$Treatment=="G3"] <- 26
unweighted$Hist_temp[unweighted$Treatment=="G4"] <- 30
unweighted$Hist_temp[unweighted$Treatment=="G5"] <- 30

unweighted <- unweighted %>% 
  select(Treatment, Temperature, Day, Container, Survivorship, Hist_temp, Initial_number) %>% 
  na.omit()

## 3. Plot survivorship over time
unweighted$Temperature <- as.factor(unweighted$Temperature)

unweighted_summary <- unweighted %>% 
  group_by(Treatment, Temperature, Day) %>% 
  summarise(Mean = mean(Survivorship), SE = sd(Survivorship)/sqrt(length(Survivorship)))

ggplot(data = unweighted_summary, mapping = aes(x = Day, y = Mean, color = Temperature)) + geom_point() +
  facet_wrap(~Treatment, scales = "free") + scale_color_manual(values = c("dodgerblue2", "red3")) +
  geom_errorbar(aes(ymin = Mean-SE, ymax = Mean+SE)) + theme_classic() + ylab("Mean Survivorship") +geom_line() +
  theme(legend.position = "None", axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
  ylim(0,1) + xlim(0,52)

#ggsave("Polyp_survivorship_over_time_SMOOTH_oct16.png", dpi = 800, height = 7, width = 10)
ggsave("Polyp_survivorship_over_time_CONNECT_feb10.pdf", height = 7 , width = 10)
#ggsave("Polyp_survivorship_over_time_CONNECT_jan18.pdf", height = 7 , width = 10)

## 4. Create subset for month1 survivorship 

month_1 <- unweighted %>% 
  filter(Day == 29|Day == 31)

## 5. Set variables as factors 

month_1$Treatment <- as.factor(month_1$Treatment)
month_1$Temperature <- as.factor(month_1$Temperature)
month_1$Hist_temp <- as.factor(month_1$Hist_temp)

ggplot(data = month_1, mapping = aes(x = Treatment, y = Survivorship, fill = Temperature))+ geom_boxplot() +
  xlab("Culture") + scale_fill_manual(values = c("dodgerblue2", "red3")) + theme_classic() + theme(legend.position = "none") +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) + ylim(0,1)

ggsave("month1_survivorship_ANOVA_2-10-21.pdf", dpi = 800, height =5, width =7)

## 6. Run ANOVA
# Where "Treatment" is symiont genotype
month_1_fit <- aov(Survivorship~Treatment*Temperature, data = month_1)

summary(month_1_fit)

## 7. Check assumptions 

# Check homogeneity of variances
resid<-month_1$Survivorship - predict.lm(month_1_fit)
plot(predict.lm(month_1_fit), resid)
abline(a=0, b= 0, col = "red")
#residuals are greater at lower values, cone shape toward higher values 

# Normality of residuals
stdRes <-rstandard(month_1_fit)
#qqplot
qqnorm(stdRes,ylab="Standardized Residuals", xlab="Theoretical Quantiles")
qqline(stdRes, col=2,lwd=2)
hist(stdRes)
# Normality is a bit off, try using a transformation (arcsine, normally used for proportions)

month_1_arcsine_fit <- aov(asin(Survivorship)~Treatment*Temperature, data = month_1)
summary(month_1_arcsine_fit)

# Check homogeneity of variances
resid<-asin(month_1$Survivorship) - predict.lm(month_1_arcsine_fit)
plot(predict.lm(month_1_arcsine_fit), resid)
abline(a=0, b= 0, col = "red")
#better

# Normality of residuals
stdRes <-rstandard(month_1_arcsine_fit)
#qqplot       
qqnorm(stdRes,ylab="Standardized Residuals", xlab="Theoretical Quantiles")
qqline(stdRes, col=2,lwd=2)
hist(stdRes)

# transformation slightly improves normality, will use this transformation 

## 8. Run Tukey Post-Hoc test 

month_1_tukey <- TukeyHSD(month_1_arcsine_fit)
month_1_tukey

month_1_HSD <- HSD.test(month_1_arcsine_fit, trt = c("Treatment","Temperature"), console = TRUE)

month_1_groups <- month_1_HSD %>% 
  .$groups %>% 
  as_tibble(rownames="class")

month_1_groups <- month_1_groups[order(month_1_groups$class),]

ggplot()+ geom_boxplot(data = month_1, mapping = aes(x = Treatment, y = Survivorship, fill = Temperature)) +
  xlab("Culture") + scale_fill_manual(values = c("dodgerblue2", "red3")) + theme_classic() + theme(legend.position = "none") +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) + ylim(0,1.05) +
  geom_text(data = month_1_groups, aes(label = groups, y= 1.05, x=c("G1", "G1", "G2", "G2", "G3", "G3", "G4", "G4", "G5", "G5")), position = position_dodge2(width = 0.7), size = 6)

ggsave("Fig4b_feb10-21.pdf", width = 8, height = 8)

## 7. Run LME model to determine if historical temperature has significant effect

month1_lme <-lme(asin(Survivorship)~Hist_temp*Temperature, random= ~1|Treatment, data=month_1, method = "REML") 
summary(month1_lme)
anova.lme(month1_lme, type="sequential", adjustSigma = FALSE)

# Check assumptions 
hist(residuals(month1_lme))
plot(fitted(month1_lme), residuals(month1_lme))

## 8. Run LM model for only fixed effects 
month1_no_culture <- gls(asin(Survivorship)~Hist_temp*Temperature, data = month_1)
summary(month1_no_culture)
anova(month1_no_culture)

# Check assumptions 
hist(residuals(month1_no_culture))
plot(fitted(month1_no_culture), residuals(month1_no_culture))

## 9. Test significance of random effect  
anova(month1_lme, month1_no_culture)

## Test for density-dep. 
## Try this: model1<- lm(survival~Genotype*Temperature + Initial_Number)

month_1_w_initial <- lm(Survivorship~Treatment*Temperature + Initial_number, data = month_1)
anova(month_1_w_initial)

# Check homogeneity of variances
resid<-month_1$Survivorship - predict.lm(month_1_w_initial)
plot(predict.lm(month_1_w_initial), resid)
abline(a=0, b= 0, col = "red")
#residuals are greater at lower values, cone shape toward higher values 

# Normality of residuals
stdRes <-rstandard(month_1_w_initial)
#qqplot
qqnorm(stdRes,ylab="Standardized Residuals", xlab="Theoretical Quantiles")
qqline(stdRes, col=2,lwd=2)
hist(stdRes)
# Normality is a bit off, try using a transformation (arcsine, normally used for proportions)
month_1_arcsine_fit <- lm(asin(Survivorship)~Treatment*Temperature+Initial_number, data = month_1)
anova(month_1_arcsine_fit)

# Check homogeneity of variances
resid<-asin(month_1$Survivorship) - predict.lm(month_1_arcsine_fit)
plot(predict.lm(month_1_arcsine_fit), resid)
abline(a=0, b= 0, col = "red")
#better

# Normality of residuals
stdRes <-rstandard(month_1_arcsine_fit)
#qqplot
qqnorm(stdRes,ylab="Standardized Residuals", xlab="Theoretical Quantiles")
qqline(stdRes, col=2,lwd=2)
hist(stdRes)

### Month 2 Survival #####################################################################################

## 1. Create subset for month2 survivorship 

month_2 <- unweighted %>% 
  filter(Day == 50|Day == 52)

## Remove contaminated containers: 

month_2 <- month_2 %>% 
  filter(Container != 122 & Container != 104 & Container != 107 & Container != 143 & Container != 214)

## 2. Set variables as factors 

month_2$Treatment <- as.factor(month_2$Treatment)
month_2$Temperature <- as.factor(month_2$Temperature)
month_2$Hist_temp <- as.factor(month_2$Hist_temp)

ggplot(data = month_2, mapping = aes(x=Treatment, y = Survivorship, fill = Temperature)) + geom_boxplot() +
  xlab("Culture") + scale_fill_manual(values = c("dodgerblue2", "red3")) + 
    theme_classic() + theme(legend.position = "none") + theme(axis.text = element_text(size = 14), axis.title = element_text(size=16))

ggsave("month_2_survivorship_oct15.png", dpi = 800, height =5, width = 7)

## 3. Run ANOVA

month_2_fit <- aov(Survivorship~Treatment*Temperature, data = month_2)
summary(month_2_fit)

## 4. Check assumptions 

# Check homogeneity of variances
resid<-month_2$Survivorship - predict.lm(month_2_fit)
plot(predict.lm(month_2_fit), resid)
abline(a=0, b= 0, col = "red")
#fine 

# Normality of residuals
stdRes <-rstandard(month_2_fit)
#qqplot
qqnorm(stdRes,ylab="Standardized Residuals", xlab="Theoretical Quantiles")
qqline(stdRes, col=2,lwd=2)
hist(stdRes)
# Normality is fine

## 5. Run Tukey Post-Hoc test 

TukeyHSD(month_2_fit)

month_2_tukey <- TukeyHSD(month_2_fit)
month_2_tukey

month_2_HSD <- HSD.test(month_2_fit, trt = c("Treatment","Temperature"), console = TRUE)

month_2_groups <- month_2_HSD %>% 
  .$groups %>% 
  as_tibble(rownames="class")

month_2_groups <- month_2_groups[order(month_2_groups$class),]

ggplot()+ geom_boxplot(data = month_2, mapping = aes(x = Treatment, y = Survivorship, fill = Temperature)) +
  xlab("Culture") + scale_fill_manual(values = c("dodgerblue2", "red3")) + theme_classic() + theme(legend.position = "none") +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) + ylim(0,1.05) +
  geom_text(data = month_2_groups, aes(label = groups, y= 1.05, x=c("G1", "G1", "G2", "G2", "G3", "G3", "G4", "G4", "G5", "G5")), position = position_dodge2(width = 0.7), size = 6)

ggsave("Fig4c_jan18-21.pdf", width = 8, height = 8)


## 6. Run LME model to determine if historical temperature has significant effect

month2_lme <-lme(Survivorship~Hist_temp*Temperature, random= ~1|Treatment, data=month_2, method = "REML") 
summary(month2_lme)
anova.lme(month2_lme, type="sequential", adjustSigma = FALSE)

# Check assumptions 
hist(residuals(month2_lme))
plot(fitted(month2_lme), residuals(month2_lme))

## 8. Run LM model for only fixed effects 
month2_no_culture <- gls(Survivorship~Hist_temp*Temperature, data = month_2)
summary(month2_no_culture)
anova(month2_no_culture)

# Check assumptions 
hist(residuals(month2_no_culture))
plot(fitted(month2_no_culture), residuals(month2_no_culture))

## 9. Test significance of random effect  
anova(month2_lme, month2_no_culture)

## Test for density-dep. 
## Try this: model1<- lm(survival~Genotype*Temperature + Initial_Number)

month_2_w_initial <- lm(Survivorship~Treatment*Temperature + Initial_number, data = month_2)
anova(month_2_w_initial)


## NATURAL WATER TEMPS ###############################################################

EL_A <- read.csv("Elbow-EL-A Temps_means.csv")
EL_A$Year <- as.factor(EL_A$Year)
EL_A_avg <- EL_A %>% 
  filter(Year != 2016 & Year != 2017 & Year != 2018 & Year != 2019 & Year != 2020) %>% 
  group_by(Month, Day, Year) %>% 
  summarize(average = mean(Temp))

ggplot(data = EL_A_avg, mapping = aes (x= Month, y = average, color = Year)) + geom_smooth() + xlim(0,12) + theme_classic() + 
  scale_color_brewer(palette= "Spectral") + xlim(1,12) + geom_line(aes(y=30), color = "red") + geom_line(aes(y=26), color = "blue") +
  ylab("Average Temperature (C)")

ggsave("EL_A_plot.pdf", width = 8, height = 8)

EL_B <- read.csv("Elbow-EL-B Temps_means.csv")
EL_B$Year <- as.factor(EL_B$Year)
EL_B_avg <- EL_B %>% 
  #filter(Year != 2016 & Year != 2017 & Year != 2018 & Year != 2019 & Year != 2020) %>% 
  group_by(Month, Day, Year) %>% 
  summarize(average = mean(Temp.C.))

ggplot(data = EL_B_avg, mapping = aes (x= Month, y = average, color = Year)) + geom_smooth() + xlim(1,12) + theme_classic() + 
  scale_color_brewer(palette= "Spectral")

ggsave("EL_B_plot.pdf", width = 8, height = 8)

EL_C <- read.csv("Elbow-EL-C Temps_means.csv")
EL_C$Year <- as.factor(EL_C$Year)
EL_C_avg <- EL_C %>% 
  #filter(Year != 2016 & Year != 2017 & Year != 2018 & Year != 2019 & Year != 2020) %>% 
  group_by(Month, Day, Year) %>% 
  summarize(average = mean(Temp.C.))

ggplot(data = EL_C_avg, mapping = aes (x= Month, y = average, color = Year)) + geom_smooth() + xlim(1,12) + theme_classic() + 
  scale_color_brewer(palette= "Spectral") + ylim(20,32)

ggsave("EL_C_plot.pdf", width = 8, height = 8)

## OSMO4	TN Reef-forereef	24.7610 -80.7417		Depth:	18.3m
tennreef <- read.csv("TennReef_5.3.21")
tennreef$Year <- as.factor(tennreef$Year)
tennreef_avg <- tennreef %>% 
  group_by(Month, Day, Year) %>% 
  summarize(average = mean(tempt))

ggplot(data = tennreef_avg, mapping = aes (x= Month, y = average, color = Year)) + geom_smooth() + xlim(0,12) + theme_classic() + 
  scale_color_brewer(palette= "Spectral") + xlim(1,12) + geom_line(aes(y=30), color = "red") + geom_line(aes(y=26), color = "blue") +
  ylab("Average Temperature (C)")

ggsave("tennreef_18m.pdf", width = 8, height = 8)
