### extract data from Murdoch1975 Guppies

### set working directory

setwd('/media/kyle/KEC_External/Documents/Data/SwitchingAndIndVar/DataFromGraphs/')

### load packages

library(dplyr)

### load data

guppy.data <- read.csv('Murdoch1975Guppies.csv')

### relative preference is 0.38 from paper

### something weird happened with the first column name

colnames(guppy.data)[1] <- 'RelativeDensity'

### transform Standard Errors into variances -- (SE*sqrt(n))^2

guppy.data$Variance <- (guppy.data$StandardError/2 * sqrt(11))^2 

### fit spline to model maximum variance

guppy.spline <- smooth.spline(x = guppy.data$RelativeDensity, y = guppy.data$Variance, df = 4)

plot(x = guppy.data$RelativeDensity, y = guppy.data$Variance)
lines(x = guppy.spline$x, y = guppy.spline$y)

guppy.predict <- predict(guppy.spline, x = seq(0, 1, 0.01))

plot(guppy.predict$x, guppy.predict$y, type = 'l')

guppy.predict$x[which(guppy.predict$y == max(guppy.predict$y))]

### now we want to calculate the expected/predicted variance at each point

calcp <- function(rel.pref, rel.density) {
  rel.pref*rel.density/(rel.pref*rel.density + ((1 - rel.pref) * (1 - rel.density)))
}

est.p <- calcp(rel.pref = 0.38, rel.density = guppy.data$RelativeDensity)

### n is 20 because fish were observed for 20 feeding observations each

n <- 20

pred.var <- est.p*(1-est.p)/n

plot(x = pred.var, y = guppy.data$Variance)
abline(a = 0, b = 1)

### get confidence intervals for estimates of variance

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

GenCI(p.vec = est.p, n.samples = 11, n.events = 20, n.sim = 10000)




