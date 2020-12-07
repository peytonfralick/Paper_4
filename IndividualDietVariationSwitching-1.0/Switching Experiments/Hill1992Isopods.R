### get data from Hill1992Isopods

### set working directory

setwd('/media/kyle/KEC_External/Documents/Data/SwitchingAndIndVar/DataFromGraphs/')

### load libraries

library(dplyr)

### load data

isopod.data <- read.csv('Hill1992Isopods.csv')

colnames(isopod.data)[1] <- 'RelativeDensity'

### alter the data

St.Dev <- function(x) {
  SS <- sum((x - mean(x))^2)
  SD <- sqrt(SS/length(x))  
  SD
}

isopod.data <- isopod.data %>% group_by(RelativeDensity) %>% summarise(MeanProportion = mean(Proportion), StandardDeviation = St.Dev(Proportion))

### transform SD to variance

isopod.data$Variance <- (isopod.data$StandardDeviation)^2 

### try to predict levels of variance

# calculate p

calcp <- function(rel.pref, rel.density) {
  rel.pref*rel.density/(rel.pref*rel.density + ((1 - rel.pref) * (1 - rel.density)))
}

isopod.p <- calcp(rel.pref = 0.54, rel.density = isopod.data$RelativeDensity)

### variance is p(1-p)/n. Here n is 10 from reported mean number of prey eaten.

isopod.est.var <- isopod.p*(1-isopod.p)/10

plot(x = isopod.est.var, y = isopod.data$Variance)
abline(a = 0, b = 1)

### confidence intervals for variance estimates

GenCI.DiffSamples <- function(p.vec, n.samples, n.events, n.sim) {
  
  quant.data <- matrix(nrow = length(p.vec), ncol = 3)
  
  for(j in 1:length(p.vec)){
    mat <- matrix(nrow = n.sim, ncol = n.samples[j])
    for(i in 1:n.sim){
      mat[i,] <- rbinom(n = n.samples[j], size = n.events, prob = p.vec[j])
    }
    prop.mat <- mat/n.events
    var.vec <- apply(prop.mat, 1, var)
    quant.data[j, ] <- c(p.vec[j], quantile(var.vec, probs = c(0.025, 0.975)))
  }
  quant.data
}

GenCI.DiffSamples(p.vec = isopod.p, n.samples = c(4,4,5), n.events = 10, n.sim = 10000)

### Estimate IS and create predicted 95% CI's for IS

### calculate IS for each relative abundance

isopod.data.ind <- read.csv('Hill1992Isopods.csv')

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
    mat <- matrix(nrow = n.sim, ncol = n.samples[j])
    PS <- matrix(nrow = n.sim, ncol = n.samples[j])
    mean.prop <- vector(length = n.sim)
    IS <- vector(length = n.sim)
    
    for(i in 1:n.sim){
      mat[i,] <- rbinom(n = n.samples[j], size = n.events, prob = p.vec[j])
      prop.mat <- mat/n.events
      mean.prop[i] <- mean(prop.mat[i])
      PS[i,] <- 1 - abs(prop.mat[i,] - mean.prop[i])
      IS[i] <- mean(PS[i,])
    }
    
    quant.data[j, ] <- c(p.vec[j], quantile(IS, probs = c(0.025, 0.95, 0.975)))
  }
  quant.data
}

GenCI_IS(p.vec = isopod.p, n.samples = c(4,4,5), n.events = 10, n.sim = 10000)
