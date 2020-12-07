### get data from Johansonn1992Dragonfly

### set working directory

setwd('D:/Documents/Data/SwitchingAndIndVar/DataFromGraphs/')

### load libraries

library(dplyr)

### load data

dragon.data <- read.csv('Johansson1992Dragonfly.csv')

colnames(dragon.data)[1] <- 'RelativeDensity'

### relative preference at equal abundance = 0.62

### convert standard deviations into variances

dragon.data$Variance <- (dragon.data$StandardDeviation)^2

### model variance to find where maximum occurs

dragon.smooth <- smooth.spline(x = dragon.data$RelativeDensity, y = dragon.data$Variance)

plot(x = dragon.data$RelativeDensity, y = dragon.data$Variance)
lines(x = dragon.data$RelativeDensity, y = dragon.smooth$y)

dragon.predict <- predict(dragon.smooth, x = seq(0,1, by = 0.01))

plot(x = seq(0,1, by = 0.01), y = dragon.predict$y)

dragon.predict$x[which(dragon.predict$y == max(dragon.predict$y))]

### try to predict the variances 

# need p's 

calcp <- function(rel.pref, rel.density) {
  rel.pref*rel.density/(rel.pref*rel.density + ((1 - rel.pref) * (1 - rel.density)))
}

dragon.p <- calcp(rel.pref = 0.62, rel.density = dragon.data$RelativeDensity)

# calculate variance p(1-p)/n. n here is 4 from functional response experiment.

dragon.variance <- dragon.p*(1-dragon.p)/4

plot(x = dragon.variance, y = dragon.data$StandardDeviation)
abline(a = 0, b = 1)


### calculate confidence interval for the variance estimates

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

GenCI(p.vec = dragon.p, n.samples = 4, n.events = 4, n.sim = 10000)


















