### getting data from Buckel 2000 on bluefish

### set working directory

setwd('/media/kyle/KEC_External/Documents/Data/SwitchingAndIndVar/DataFromGraphs/')

### load libraries

library(dplyr)

### load data for bluefish with no seagrass

Bluefish.NoSeagrass <- read.csv('Buckel2000BluefishNoSeagrass.csv')

### load data for bluefish with seagrass

Bluefish.Seagrass <- read.csv('Buckel2000BluefishSeagrass.csv')

colnames(Bluefish.Seagrass)[1] <- 'RelativeDensity'

# relative preferences from 1:1 treatment

# no seagrass = 0.58

# seagrass = 0.54

### Need to calculate variances from standard errors

# standard error is standard deviation divided by sqrt of the sample size SE = sigma/sqrt(n)
# So, variance is (SE*sqrt(n))^2

Bluefish.NoSeagrass$Variance <- (Bluefish.NoSeagrass$StandardError*sqrt(3))^2

Bluefish.Seagrass$Variance <- (Bluefish.Seagrass$StandardError*sqrt(3))^2

### fit splines to data to find where the maximum variance occurs

Bluefish.NoSeagrass.spline <- smooth.spline(x = Bluefish.NoSeagrass$RelativeDensity, y = Bluefish.NoSeagrass$Variance)

plot(x = Bluefish.NoSeagrass$RelativeDensity, y = Bluefish.NoSeagrass$Variance)
lines(x = Bluefish.NoSeagrass.spline$x, y = Bluefish.NoSeagrass.spline$y)

Bluefish.NoSeagrass.predict <- predict(Bluefish.NoSeagrass.spline, x = seq(0, 1, by = 0.01))

plot(x = seq(0, 1, by = 0.01), y = Bluefish.NoSeagrass.predict$y)

### where does maximum occur?

Bluefish.NoSeagrass.predict$x[which(Bluefish.NoSeagrass.predict$y == max(Bluefish.NoSeagrass.predict$y))]

### now for seagrass

Bluefish.Seagrass.spline <- smooth.spline(x = Bluefish.Seagrass$RelativeDensity, y = Bluefish.Seagrass$Variance)

plot(x = Bluefish.Seagrass$RelativeDensity, y = Bluefish.Seagrass$Variance)
lines(x = Bluefish.Seagrass.spline$x, y = Bluefish.Seagrass.spline$y)

Bluefish.Seagrass.predict <- predict(Bluefish.Seagrass.spline, x = seq(0, 1, by = 0.01))

plot(x = seq(0, 1, by = 0.01), y = Bluefish.NoSeagrass.predict$y)

### where does the maximum occur?

Bluefish.Seagrass.predict$x[which(Bluefish.Seagrass.predict$y == max(Bluefish.Seagrass.predict$y))]

### Now want to try to predict the magnitude of the variances 

# function to calculate p's

calcp <- function(rel.pref, rel.density) {
  rel.pref*rel.density/(rel.pref*rel.density + ((1 - rel.pref) * (1 - rel.density)))
}

# no seagrass

Bluefish.NoSeagrass.p <- calcp(rel.pref = 0.58, rel.density = Bluefish.NoSeagrass$RelativeDensity)

# seagrass

Bluefish.Seagrass.p <- calcp(rel.pref = 0.54, rel.density = Bluefish.Seagrass$RelativeDensity)

### what are the n's? Functional response experiments suggest 12 ... 

n <- 12

### predicted variances 

# no seagrass

NoSeagrass.pred.var <- Bluefish.NoSeagrass.p*(1-Bluefish.NoSeagrass.p)/12

plot(x = NoSeagrass.pred.var, y = Bluefish.NoSeagrass$Variance)
abline(a = 0, b = 1)

# seagrass

Seagrass.pred.var <- Bluefish.Seagrass.p*(1 - Bluefish.Seagrass.p)/12

plot(x = Seagrass.pred.var, y = Bluefish.Seagrass$Variance)
abline(a = 0, b = 1)

### calculate CI's for the predicted variances

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

# No Seagrass

GenCI(p.vec = Bluefish.NoSeagrass.p, n.samples = 3, n.events = 12, n.sim = 10000)

# Seagrass

GenCI(p.vec = Bluefish.Seagrass.p, n.samples = 3, n.events = 12, n.sim = 10000)





















