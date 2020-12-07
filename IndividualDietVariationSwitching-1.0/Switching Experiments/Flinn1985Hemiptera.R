### pulling data from Flinn1985Hemiptera

### set working directory

setwd('D:/Documents/Data/SwitchingAndIndVar/DataFromGraphs/')

### load libraries

library(dplyr)

### load data

bug.data <- read.csv('Flinn1985Hemiptera.csv')

# first column name weird again

colnames(bug.data)[1] <- 'RelativeDensities'

### relative preference is 0.12

### convert standard errors to variances
# variance = (SE*sqrt(n))^2
# n in this study is 6

bug.data$variance <- (bug.data$StandardError * sqrt(6))^2

# now want to fit a spline to the data

bug.spline <- smooth.spline(x = bug.data$RelativeDensities, y = bug.data$variance)

plot(x = bug.data$RelativeDensities, y = bug.data$variance)
lines(x = bug.data$RelativeDensities, y = bug.spline$y)

bug.predict <- predict(bug.spline, x = seq(0, 1, by = 0.01))

plot(x = seq(0,1, by = 0.01), y = bug.predict$y)

bug.predict$x[which(bug.predict$y == max(bug.predict$y))]

### now want to try to predict the magnitude of the variance at each ratio of prey

### want to calculate p for each of the relative densities

calcp <- function(rel.pref, rel.density) {
  rel.pref*rel.density/(rel.pref*rel.density + ((1 - rel.pref) * (1 - rel.density)))
}

bug.p.est <- calcp(rel.pre = 0.12, rel.density = bug.data$RelativeDensities)

### calculate predicted variance. n here is 8 from functional response experiment

bug.pred.var <- bug.p.est*(1-bug.p.est)/8

plot(x = bug.pred.var, y = bug.data$variance)
abline(a = 0, b = 1)

### calculate confidence intervals for each variance estimate

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

GenCI(p.vec = bug.p.est, n.samples = 6, n.events = 8, n.sim = 10000)














