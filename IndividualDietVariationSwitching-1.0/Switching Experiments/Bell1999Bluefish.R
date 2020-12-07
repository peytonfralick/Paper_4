### extracting data from Bell1999Bluefish
### set working directory

setwd('/media/kyle/KEC_External/Documents/Data/SwitchingAndIndVar/DataFromGraphs/')

### load libraries

library(dplyr)

### load data

bluefish.data <- read.csv('Bell1999Bluefish.csv')

colnames(bluefish.data)[1] <- 'RelativeDensity'

### relative preference is 0.18 from the 1:1 trial

### estimate the magnitude of variance at each of the densities

# get the variances from the standard errors
# var = (SE * sqrt(n))^2. n in this study is 3.

bluefish.data$Variance <- (bluefish.data$StandardError * sqrt(3))^2

# calculate what the expected variance is. Variance of a proportion is
# p(1-p)/n where n is the number of trials/number prey eaten.

# calculate p for each of the relative densities

calcp <- function(rel.pref, rel.density) {
  rel.pref*rel.density/(rel.pref*rel.density + ((1 - rel.pref) * (1 - rel.density)))
}

estimated.p <- calcp(rel.pref = 0.18, rel.density = bluefish.data$RelativeDensity)

# number of total feeding events for each trial is 8, from data reported on number of attacks and proportion of successes

pred.var <- estimated.p*(1-estimated.p)/8

# examine relationship between observed and predicted variances

plot(x = pred.var, y = bluefish.data$Variance)
abline(a = 0, b = 1)

### get confidence intervals for variance estimates

# function to calculate confidence intervals

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

GenCI(p.vec = estimated.p, n.samples = 3, n.events = 8, n.sim = 10000)
























