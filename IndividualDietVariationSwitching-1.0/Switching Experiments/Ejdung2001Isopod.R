### get data from Ejdung2001Isopod

### set working directory

setwd('D:/Documents/Data/SwitchingAndIndVar/DataFromGraphs/')

### load libraries

library(dplyr)

### load data

isopod.data <- read.csv('Ejdung2001Isopod.csv')

colnames(isopod.data)[1] <- 'RelativeDensity'

### alter data to have standard deviation and proportion

St.Dev <- function(x) {
  SS <- sum((x - mean(x))^2)
  SD <- sqrt(SS/length(x))  
  SD
}

isopod.data <- isopod.data %>% group_by(RelativeDensity) %>% summarise(MeanProportion = mean(Proportion), StandardDeviation = St.Dev(Proportion))

### relative preference is 0.9 from equal abundance treatment

### convert SD to variance

isopod.data$Variance <- (isopod.data$StandardDeviation)^2

### calculate p

calcp <- function(rel.pref, rel.density) {
  rel.pref*rel.density/(rel.pref*rel.density + ((1 - rel.pref) * (1 - rel.density)))
}

isopod.p <- calcp(rel.pref = 0.9, rel.density = isopod.data$RelativeDensity)

# variance is p(1-p)/n. n here is 10 from prey depletion experiment

isopod.var.est <- isopod.p*(1 - isopod.p)/10

plot(x = isopod.var.est, y = isopod.data$Variance)
abline(a = 0, b = 1)

### calculate CI's for each variance estimate

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

GenCI(p.vec = isopod.p, n.samples = 6, n.events = 10, n.sim = 10000)

### Estimate IS and create predicted 95% CI's for IS

### Calculate IS for each relative abundance

isopod.data.ind <- read.csv('Ejdung2001Isopod.csv')

colnames(isopod.data.ind)[1] <- 'RelativeDensity'

isopod.data.PS <- isopod.data.ind %>% mutate(PropAlt = 1 - Proportion)

isopod.data.PS <- isopod.data.PS %>% group_by(RelativeDensity) %>% mutate(Diff.Prop = abs(Proportion - mean(Proportion)), 
                                                                          Diff.Alt = abs(PropAlt - mean(PropAlt)))

isopod.data.PS <- isopod.data.PS %>% mutate(PS = 1 - Diff.Prop)

isopod.data.IS <- isopod.data.PS %>% group_by(RelativeDensity) %>% summarise(IS = mean(PS))

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

GenCI_IS(p.vec = isopod.p, n.samples = 6, n.events = 10, n.sim = 10000)
