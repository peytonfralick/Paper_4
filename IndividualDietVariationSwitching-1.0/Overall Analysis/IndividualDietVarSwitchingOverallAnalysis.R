###############################################################
### Analysis of data extracted from switching experiments
###############################################################

### clear workspace

rm(list = ls())

### load libraries

library(dplyr)

### set working directory

setwd('/media/kyle/KEC_External/Documents/Data/SwitchingAndIndVar/DataAnalysis/')

#################################################################
### Analysis 1: Predicting the relative abundance at which the 
### maximum diet variation occurs
#################################################################

### load data

max.data <- read.csv('MaximumData.csv')

### split data into studies that only considered 3 relative abundances
### and those that considered more than three relative abundances

max.3.data <- max.data %>% filter(RelativeDensityNumber == '3')

### calculate R^2 of 1:1 line between predicted relative abundance and observed

### write function to calculate R^2 of 1:1 line

calcr2 <- function(expected, observed){
  SS.res <- sum((expected - observed)^2)
  SS.tot <- sum((observed - mean(observed))^2)
  1 - (SS.res/SS.tot)
}

calcr2(max.3.data$ExpectedMaximum, max.3.data$ObservedMaximum)

### r2 = 0.33

### do the same for studies with more than 3 relative abundances

max.mt3.data <- max.data %>% filter(RelativeDensityNumber == ">3")

### get an r2

calcr2(expected = max.mt3.data$ExpectedMaximum, observed = max.mt3.data$ObservedMaximum)

### r2 = 0.79

##############################################################################
### Analysis 2: Predicting the magnitude of variation
##############################################################################

### load magnitude data

mag.data <- read.csv('MagnitudeData.csv')

### 161 total observations

### correlation coefficient between predicted variance and observed variance

cor(x = mag.data$EstimatedVariance, y = mag.data$ObservedVariance)

### rho = 0.59

### calculate proportion of observed variances falling within the 95% CI's

Covered <- ifelse(mag.data$ObservedVariance >= mag.data$LowerBound & mag.data$ObservedVariance <= mag.data$UpperBound, 1, 0)

sum(Covered)/161

### 79.5% of observations within 95% CI's

### number underestimated

sum(ifelse(mag.data$ObservedVariance > mag.data$UpperBound, 1, 0))

### 25

### number overestimated

sum(ifelse(mag.data$ObservedVariance < mag.data$LowerBound, 1, 0))

### 8









