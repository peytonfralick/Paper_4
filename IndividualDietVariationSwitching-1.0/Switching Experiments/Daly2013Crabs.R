### get data from Daly2013Crabs

### set working directory

setwd('/media/kyle/KEC_External/Documents/Data/SwitchingAndIndVar/DataFromGraphs/')

### load libraries

library(dplyr)

### load data

crab.data <- read.csv('Daly2013Crabs.csv')

colnames(crab.data)[1] <- 'Experiment'

### break up experiments

crab.cobble <- crab.data %>% filter(Experiment == 'Cobble')

crab.shell <- crab.data %>% filter(Experiment == 'Shell')

# relative preference cobble is 0.52

# relative preference shell is 0.49

### calculate variances

# n's for cobble

n.cobble <- c(4,4,9,9,5)

# n's for shell

n.shell <- c(4,5,12,5,6) 

### variance is (se * sqrt(n))^2

# cobble

crab.cobble$Variance <- (crab.cobble$StandardError * sqrt(n.cobble))^2

# shell

crab.shell$Variance <- (crab.shell$StandardError * sqrt(n.shell))^2

### model variance to find where the maximum was 

# cobble

cobble.smooth <- smooth.spline(x = crab.cobble$RelativeDensity, y = crab.cobble$Variance, df = 3)

plot(x = crab.cobble$RelativeDensity, y = crab.cobble$Variance)
lines(x = crab.cobble$RelativeDensity, y = cobble.smooth$y)

cobble.predict <- predict(cobble.smooth, x = seq(0,1, by = 0.01))

plot(x = seq(0,1, by = 0.01), y = cobble.predict$y)

cobble.predict$x[which(cobble.predict$y == max(cobble.predict$y))]

# shell

shell.smooth <- smooth.spline(x = crab.shell$RelativeDensity, y = crab.shell$Variance)

plot(x = crab.shell$RelativeDensity, y = crab.shell$Variance)
lines(x = crab.shell$RelativeDensity, y = shell.smooth$y)

shell.predict <- predict(shell.smooth, x = seq(0,1, by = 0.01))

plot(x = seq(0,1, by = 0.01), y = shell.predict$y)

shell.predict$x[which(shell.predict$y == max(shell.predict$y))]

### try to predict variances

### for cobble

# calculate p

calcp <- function(rel.pref, rel.density) {
  rel.pref*rel.density/(rel.pref*rel.density + ((1 - rel.pref) * (1 - rel.density)))
}

cobble.p <- calcp(rel.pref = 0.52, rel.density = crab.cobble$RelativeDensity)

# calculate variance -- p(1-p)/n. n is 4 from functional response experiments

cobble.variance <- cobble.p*(1 - cobble.p)/4

plot(x = cobble.variance, y = crab.cobble$Variance)
abline(a = 0, b = 1)


### for shell

shell.p <- calcp(rel.pref = 0.49, rel.density = crab.shell$RelativeDensity)

# calculate variance -- p(1-p)/n. n is 4 from functional response experiments.

shell.variance <- shell.p*(1-shell.p)/4

plot(x = shell.variance, y = crab.shell$Variance)
abline(a = 0, b = 1)

# calculate CI's for each of the estimated variances

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

### Cobble

GenCI.DiffSamples(p.vec = cobble.p, n.samples = n.cobble, n.events = 4, n.sim = 10000)

### Shell

GenCI.DiffSamples(p.vec = shell.p, n.samples = n.shell, n.events = 4, n.sim = 10000 )















