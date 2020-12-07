#############################################################################
### Analysis of IS magnitude data for Appendix 2
#############################################################################

# load libraries

library(dplyr); 

# set working directory

setwd('D:/Documents/Data/SwitchingAndIndVar/DataAnalysis/')

###########################################################
### relationship between variance and IS
###########################################################

### get magnitude data

Mag.data <- read.csv('MagnitudeData.csv')

### subset studies for which IS could be estimated

IS.data <- Mag.data %>% filter(EstimatedIS > 0)

### minimum maximum and mean of IS

min(IS.data$ObservedIS)

max(IS.data$ObservedIS)

mean(IS.data$ObservedIS)

### correlation between variance and IS

summary(lm(ObservedIS ~ ObservedVariance, data = IS.data))  
cor.test(x = IS.data$ObservedVariance, y = IS.data$ObservedIS)

### observed versus predicted IS values

cor(x = IS.data$EstimatedIS, y = IS.data$ObservedIS)
summary(lm(EstimatedIS ~ ObservedIS, data = IS.data))

### Coverage of observed IS's by 95% CI's 

Covered <- ifelse(IS.data$ObservedIS >= IS.data$LowerIS & IS.data$ObservedIS <= IS.data$UpperIS, 1, 0)

sum(Covered)/66

sum(Covered)

