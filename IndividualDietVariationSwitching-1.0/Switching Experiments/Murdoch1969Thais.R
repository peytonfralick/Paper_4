### get data from Murdoch1969Thais

### set working directory

setwd('/media/kyle/KEC_External/Documents/Data/SwitchingAndIndVar/DataFromGraphs/')

### load libraries

library(dplyr)

### load data

whelk.data <- read.csv('Murdoch1969Thais.csv')

### Relative pref. of whelk -- 0.92

### calculate variance

### variance is (SE * sqrt(n))^2. n here is 5

whelk.data$Variance <- (whelk.data$StandardError * sqrt(5))^2

### maximum variance occurs lowest observed relative abundance. Spline would say maximum variance
### should occur at zero. Use the lowest observed relative abundance as the observed 
### maximum variation.

### now try to estimate the magnitude of variance

### first calculate p

calcp <- function(rel.pref, rel.density) {
  rel.pref*rel.density/(rel.pref*rel.density + ((1 - rel.pref) * (1 - rel.density)))
}

whelk.p <- calcp(rel.pref = 0.9, rel.density = whelk.data$RelativeDensity)

### calculate variance 

### variance = p(1-p)/n. n here is 5 from functional response.

whelk.est.var <- whelk.p*(1-whelk.p)/5

plot(x = whelk.est.var, y = whelk.data$Variance)
abline(a = 0, b = 1)

### calculate confidence intervals for the variance estimates

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




















































