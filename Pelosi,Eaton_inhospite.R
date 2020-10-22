### Analyses for Pelosi, Eaton, Mychajliw, terHorst, and Coffroth 
### Thermal tolerance in an octocoral symbiont and the implications for
###    Caribbean corals through in hospite studies  
###

library(ggplot2)
library(nlme)
library (dplyr)

### Symbiont Growth Rates

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

### Polyp-Symbiont Cell Density 

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

ggsave("polyp_cell_density_oct16.png", dpi = 800, height = 5, width =7)

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

### Polyp survivorship 

## 1. Import data 
unweighted <- read.csv("../../polyp_analyses/unweighted average polyp survivorship.csv")

## 2. Add column for historical temp and remove extra columns

unweighted$Hist_temp[unweighted$Treatment=="G1"] <- 26
unweighted$Hist_temp[unweighted$Treatment=="G2"] <- 26
unweighted$Hist_temp[unweighted$Treatment=="G3"] <- 26
unweighted$Hist_temp[unweighted$Treatment=="G4"] <- 30
unweighted$Hist_temp[unweighted$Treatment=="G5"] <- 30

unweighted <- unweighted %>% 
  select(Treatment, Temperature, Day, Container, Survivorship, Hist_temp) %>% 
  na.omit()

## 3. Plot survivorship over time
unweighted$Temperature <- as.factor(unweighted$Temperature)

unweighted_summary <- unweighted %>% 
  group_by(Treatment, Temperature, Day) %>% 
  summarise(Mean = mean(Survivorship), SE = sd(Survivorship)/sqrt(length(Survivorship)))

ggplot(data = unweighted_summary, mapping = aes(x = Day, y = Mean, color = Temperature)) + geom_point() +
  facet_wrap(~Treatment, scales = "free") + scale_color_manual(values = c("dodgerblue2", "red3")) +
  geom_errorbar(aes(ymin = Mean-SE, ymax = Mean+SE)) + theme_classic() + ylab("Mean Survivorship") +geom_line() +
  theme(legend.position = "None", axis.text = element_text(size = 14), axis.title = element_text(size = 16))

ggsave("Polyp_survivorship_over_time_SMOOTH_oct16.png", dpi = 800, height = 7, width = 10)
ggsave("Polyp_survivorship_over_time_CONNECT_oct15.png", dpi = 300, height = 7 , width = 10)

## 4. Create subset for month1 survivorship 

month_1 <- unweighted %>% 
  filter(Day == 28|Day == 29)

## 5. Set variables as factors 

month_1$Treatment <- as.factor(month_1$Treatment)
month_1$Temperature <- as.factor(month_1$Temperature)
month_1$Hist_temp <- as.factor(month_1$Hist_temp)

ggplot(data = month_1, mapping = aes(x = Treatment, y = Survivorship, fill = Temperature))+ geom_boxplot() +
  xlab("Culture") + scale_fill_manual(values = c("dodgerblue2", "red3")) + theme_classic() + theme(legend.position = "none") +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16))

ggsave("month1_survivorship_oct15.png", dpi = 800, height =5, width =7)

## 6. Run ANOVA

month_1_fit <- aov(Survivorship~reatment*Temperature, data = month_1)

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

TukeyHSD(month_1_arcsine_fit)

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

### Month 2 Survival 


## 1. Create subset for month1 survivorship 

month_2 <- unweighted %>% 
  filter(Day == 59|Day == 60)

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




