### get data from Katz1985Whelks

### set working directory

setwd('D:/Documents/Data/SwitchingAndIndVar/DataFromGraphs/')

### load libraries

library(dplyr)

### load data

whelk.data <- read.csv('Katz1985Whelks.csv')

colnames(whelk.data)[1] <- 'RelativeDensity'

### relative preference is 0.97

### calculte variance from standard error
### VAR = (SE * sqrt(n))^2

whelk.data$Variance <- (whelk.data$StandardError * sqrt(5))^2

### no need to model maximum variance

### try to estimate variance at each relative density

calcp <- function(rel.pref, rel.density) {
  rel.pref*rel.density/(rel.pref*rel.density + ((1 - rel.pref) * (1 - rel.density)))
}

whelk.p <- calcp(rel.pref = 0.97, rel.density = whelk.data$RelativeDensity)

### variance is p(1-p)/n. n here is 6 from functional response experiment.

whelk.var.est <- whelk.p*(1-whelk.p)/6

plot(x = whelk.var.est, y = whelk.data$Variance)
abline(a = 0, b = 1)

### get confidence intervals for each of the variances

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

GenCI(p.vec = whelk.p, n.samples = 5, n.events = 6, n.sim = 10000)




















