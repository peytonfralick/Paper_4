### get data from Butler1991MitesChydorus

### set working directory

setwd('/media/kyle/KEC_External/Documents/Data/SwitchingAndIndVar/DataFromGraphs/')

### load libraries

library(dplyr)

### load data

mite.data <- read.csv('Butler1991MitesChydorus.csv')

colnames(mite.data)[1] <- 'RelativeDensity'

### function to calculate standard deviation

St.Dev <- function(x) {
  SS <- sum((x - mean(x))^2)
  SD <- sqrt(SS/length(x))  
  SD
}

mite.data <- mite.data %>% group_by(RelativeDensity) %>% summarise(MeanProportion = mean(Proportion), StandardDeviation = St.Dev(Proportion))

### relative preference is 0.51

### model where the maximum variance occurs

# calculate variance from sd

mite.data$Variance <- (mite.data$StandardDeviation)^2

# fit spline

mite.smooth <- smooth.spline(x = mite.data$RelativeDensity, y = mite.data$Variance)

plot(x = mite.data$RelativeDensity, y = mite.data$Variance)
lines(x = mite.data$RelativeDensity, y = mite.smooth$y)

mite.predict <- predict(mite.smooth, x = seq(0, 1, by = 0.01))

plot(x = seq(0, 1, by = 0.01), y = mite.predict$y)

mite.predict$x[which(mite.predict$y == max(mite.predict$y))]

### try to predict the magnitude of variance

### function to calculate p's for each relative abundance

calcp <- function(rel.pref, rel.density) {
  rel.pref*rel.density/(rel.pref*rel.density + ((1 - rel.pref) * (1 - rel.density)))
}

mite.p <- calcp(rel.pref = 0.51, rel.density = mite.data$RelativeDensity)

### calculate variance -- p(1-p)/n. n here is 15 from reported number of prey eaten.

mite.var.est <- mite.p*(1-mite.p)/15

plot(x = mite.var.est, y = mite.data$Variance)
abline(a = 0, b = 1)

### calculate CI's for estimates of the magnitude of the variance

GenCI <- function(p.vec, n.samples, n.events, n.sim) {
  
  quant.data <- matrix(nrow = length(p.vec), ncol = 3)
  
  for(j in 1:length(p.vec)){
    mat <- matrix(nrow = n.sim, ncol = n.samples)
    for(i in 1:n.sim){
      mat[i,] <- rbinom(n = n.samples, size = n.events, prob = p.vec[j])
    }
    prop.mat <- mat/n.events
    var.vec <- apply(prop.mat, 1, var)
    quant.data[j, ] <- c(p.vec[j], quantile(var.vec, probs = c(0.025, 0.975)))
  }
  quant.data
}

GenCI(p.vec = mite.p, n.samples = 4, n.events = 15, n.sim = 10000)

### Estimate IS and create predicted 95% CI's for IS

mite.data.ind <- read.csv('Butler1991MitesChydorus.csv')

colnames(mite.data.ind)[1] <- 'RelativeDensity'

mite.data.PS <- mite.data.ind %>% mutate(PropAlt = 1 - Proportion)

mite.data.PS <- mite.data.PS %>% group_by(RelativeDensity) %>% mutate(Diff.Prop = abs(Proportion - mean(Proportion)), 
                                                                      Diff.Alt = abs(PropAlt - mean(PropAlt)))

mite.data.PS <- mite.data.PS %>% mutate(PS = 1 - Diff.Prop)

mite.data.IS <- mite.data.PS %>% group_by(RelativeDensity) %>% summarise(IS = mean(PS))

### calculate confidence intervals for IS estimates using modified version of CI function above

GenCI_IS <- function(p.vec, n.samples, n.events, n.sim) {
  
  quant.data <- matrix(nrow = length(p.vec), ncol = 4)
  
  for(j in 1:length(p.vec)){
    mat <- matrix(nrow = n.sim, ncol = n.samples)
    PS <- matrix(nrow = n.sim, ncol = n.samples)
    mean.prop <- vector(length = n.sim)
    IS <- vector(length = n.sim)
    
    for(i in 1:n.sim){
      mat[i,] <- rbinom(n = n.samples, size = n.events, prob = p.vec[j])
      prop.mat <- mat/n.events
      mean.prop[i] <- mean(prop.mat[i])
      PS[i,] <- 1 - abs(prop.mat[i,] - mean.prop[i])
      IS[i] <- mean(PS[i,])
    }
    
    quant.data[j, ] <- c(p.vec[j], quantile(IS, probs = c(0.025, 0.5, 0.975)))
  }
  quant.data
}

GenCI_IS(p.vec = mite.p, n.samples = 4, n.events = 15, n.sim = 10000)
