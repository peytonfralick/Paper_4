### get data from Murdoch1969Acanthina

### set working directory

setwd('/media/kyle/KEC_External/Documents/Data/SwitchingAndIndVar/DataFromGraphs/')

### load libraries

library(dplyr)

### load data

whelk.data <- read.csv('Murdoch1969Acanthina.csv')

### calculate standard deviation

St.Dev <- function(x) {
  SS <- sum((x - mean(x))^2)
  SD <- sqrt(SS/length(x))  
  SD
}

whelk.data <- whelk.data %>% group_by(RelativeDensity) %>% summarise(MeanProportion = mean(Proportion), StandardDeviation = St.Dev(Proportion))

### calculate variance from standard deviation

whelk.data$Variance <- (whelk.data$StandardDeviation)^2

### relative preference -- 0.53 from proportion at equal abundance

### model variance with spline to find where maximum occurs

whelk.smooth <- smooth.spline(x = whelk.data$RelativeDensity, y = whelk.data$Variance, df = 4)

plot(x = whelk.data$RelativeDensity, y = whelk.data$Variance)
lines(x = whelk.data$RelativeDensity, y = whelk.smooth$y)

whelk.predict <- predict(whelk.smooth, x = seq(0,1, by = 0.01))

plot(x = seq(0, 1, by = 0.01), y = whelk.predict$y)

whelk.predict$x[which(whelk.predict$y == max(whelk.predict$y))]

### now want to try to predict the magnitude of variance

### first calculate the p's 

calcp <- function(rel.pref, rel.density) {
  rel.pref*rel.density/(rel.pref*rel.density + ((1 - rel.pref) * (1 - rel.density)))
}

whelk.p <- calcp(rel.pref = 0.53, rel.density = whelk.data$RelativeDensity)

### variance is p(1-p)/n. n here is 24 from total number of prey eaten.

whelk.var.est <- whelk.p*(1 - whelk.p)/24

plot(x = whelk.var.est, y = whelk.data$Variance)
abline(a = 0, b = 1)

### find confidence intervals for variance estimates

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

GenCI(p.vec = whelk.p, n.samples = 5, n.events = 5, n.sim = 10000)

### Estimate IS and create predicted 95% CI's for IS

### calculate IS for each relative abundance

whelk.data.ind <- read.csv('Murdoch1969Acanthina.csv')

whelk.data.PS <- whelk.data.ind %>% mutate(PropAlt = 1 - Proportion)

whelk.data.PS <- whelk.data.PS %>% group_by(RelativeDensity) %>% mutate(Diff.Prop = abs(Proportion - mean(Proportion)), 
                                                                        Diff.Alt = abs(PropAlt - mean(PropAlt)))

whelk.data.PS <- whelk.data.PS %>% mutate(PS = 1 - Diff.Prop)

whelk.data.IS <- whelk.data.PS %>% group_by(RelativeDensity) %>% summarise(IS = mean(PS))

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

GenCI_IS(p.vec = whelk.p, n.samples = 5, n.events = 5, n.sim = 10000)
