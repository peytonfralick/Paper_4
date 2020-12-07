### Get data from Cuthbert2018Amphipod

### set working directory

setwd('/media/kyle/KEC_External/Documents/Data/SwitchingAndIndVar/DataFromGraphs/')

### load libraries

library(dplyr)

### load data

amphi.data <- read.csv('Cuthbert2018Amphipod.csv')

### weird ... column name needs changed, need to drop a bunch of rows

colnames(amphi.data)[1] <- 'RelativeDensity'

amphi.data <- amphi.data[-8:-24,]

# relative preference is 0.8 where relative abundance is 0.5

### convert standard error to variance

# variance = (SE * sqrt(n))^2, n here is 6

amphi.data <- amphi.data %>% mutate(Variance = (StandardError * sqrt(6))^2)

# model variance to find the maximum

amphi.smooth <- smooth.spline(x = amphi.data$RelativeDensity, y = amphi.data$Variance)

plot(x = amphi.data$RelativeDensity, y = amphi.data$Variance)
lines(x = amphi.data$RelativeDensity, y = amphi.smooth$y)

amphi.predict <- predict(amphi.smooth, x = seq(0, 1, by = 0.01))

plot(x = seq(0,1, by = 0.01), y = amphi.predict$y)

amphi.predict$x[which(amphi.predict$y == max(amphi.predict$y))]

### next step: try to predict the variances

# function to calculate p's at each relative abundance

calcp <- function(rel.pref, rel.density) {
  rel.pref*rel.density/(rel.pref*rel.density + ((1 - rel.pref) * (1 - rel.density)))
}

amphi.p.est <- calcp(rel.pre = 0.8, rel.density = amphi.data$RelativeDensity)

# variance estimate is p(1-p)/n
# n for this is approximately 12/8

var.est <- amphi.p.est*(1-amphi.p.est)/8

plot(x = var.est, y = amphi.data$Variance)
abline(a = 0, b = 1)

### calculate CI for variance

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

GenCI(p.vec = amphi.p.est, n.samples = 6, n.events = 8, n.sim = 10000)



