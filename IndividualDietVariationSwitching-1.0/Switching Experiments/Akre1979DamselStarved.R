### extracting data from Aker1979DamselStarved

### set working directory

setwd('/media/kyle/KEC_External/Documents/Data/SwitchingAndIndVar/DataFromGraphs/')

### load packages

library(dplyr)


# Convert c to relative preference

# write a function to convert and then give relative preference and predicted maximum variance ratio

ConvertC <- function(c) {
  RelPref <- c/(1+c)
  MaxVar <- 1 - RelPref
  print(as.data.frame(cbind(RelPref, MaxVar)))
}

ConvertC(0.24)

# load data

damsel.data <- read.csv('Akre1979DamselStarved.csv')

# Need to calculate the variance/sd at each ratio of prey

# first, make a vector of the relative abundances
# and an associated factor

damsel.data$Real.Relative.Density <- rep(c(3/66, 15/66, 27/66, 39/66, 51/66, 63/66), each = 5)

damsel.data$Factor.Relative.Density <- rep(c('3/66', '15/66', '27/66', '39/66', '51/66', '63/66'), each = 5)

# need to calculate variance/sd for each of the factors

St.Dev <- function(x) {
  SS <- sum((x - mean(x))^2)
  SD <- sqrt(SS/length(x))  
  SD
}


damsel.StandardDeviation <- damsel.data %>% group_by(Factor.Relative.Density) %>% summarise(StandardDeviation = St.Dev(Proportion.in.Diet), Relative.Density = mean(Real.Relative.Density))

damsel.variance <- damsel.StandardDeviation %>% mutate(Variance = StandardDeviation^2)

# fit a spline to variance data

damsel.spline <- smooth.spline(x = damsel.variance$Relative.Density, y = damsel.variance$Variance)


plot(damsel.variance$Relative.Density, damsel.variance$Variance)
lines(damsel.spline$x, damsel.spline$y)

# now predict new values using the spline

damsel.prediction <- predict(damsel.spline, x = seq(from = 0, to = 1, 0.01))

plot(damsel.prediction$x, damsel.prediction$y, type = "l")

# can now get the empirical maximum variance

damsel.prediction$x[which(damsel.prediction$y == max(damsel.prediction$y))]

### Last thing we want to do is to calculate the expected variance at each ratio and compare that to the 
### observed variance at each ratio

### predicted variance is p(1-p)/n where n is the number of feedings

# first write function to calculate p for each relative density

calcp <- function(rel.pref, rel.density) {
  rel.pref*rel.density/(rel.pref*rel.density + ((1 - rel.pref) * (1 - rel.density)))
}

est.p <- calcp(rel.pref = 0.19, rel.density = damsel.variance$Relative.Density)

# next, make vector of n's from figure 6

n <- c(47, 40, 46, 43, 23, 21)

# calculate predicted variances

pred.var <- est.p*(1-est.p)/n

# examine relationship between predicted and observed

plot(x = pred.var, y = damsel.variance$Variance)
abline(a = 0, b = 1)

### generate confidence intervals for predicted variance

GenCI_EventNsVector <- function(p.vec, n.samples, n.events, n.sim) {
  
  quant.data <- matrix(nrow = length(p.vec), ncol = 3)
  
  for(j in 1:length(p.vec)){
    mat <- matrix(nrow = n.sim, ncol = n.samples)
    for(i in 1:n.sim){
      mat[i,] <- rbinom(n = n.samples, size = n.events[j], prob = p.vec[j])
    }
    prop.mat <- mat/n.events[j]
    var.vec <- apply(prop.mat, 1, var)
    quant.data[j, ] <- c(p.vec[j], quantile(var.vec, probs = c(0.025, 0.975)))
  }
  quant.data
}

GenCI_EventNsVector(p.vec = est.p, n.samples = 5, n.events = n, n.sim = 10000)

# output gives p's and CI's

### Estimate IS and create predicted 95% CI's for IS

damsel.data.PS <- damsel.data %>% mutate(Proportion.Alternative = 1 - Proportion.in.Diet)

damsel.data.PS <- damsel.data.PS %>% group_by(Factor.Relative.Density) %>% mutate(Diff.Prop = abs(Proportion.in.Diet - mean(Proportion.in.Diet)), 
                                                                                  Diff.Alt = abs(Proportion.Alternative - mean(Proportion.Alternative))) 

damsel.data.PS <- damsel.data.PS %>% mutate(PS = 1 - Diff.Prop)

damsel.data.IS <- damsel.data.PS %>% group_by(Factor.Relative.Density) %>% summarise(IS = mean(PS))  

# use predicted mean proportions to generate 95% CI's for IS at each relative abundance

# modify confidence interval function above to get CI's for IS instead of variance

GenCI_IS <- function(p.vec, n.samples, n.events, n.sim) {
  
  quant.data <- matrix(nrow = length(p.vec), ncol = 4)
  
  for(j in 1:length(p.vec)){
    mat <- matrix(nrow = n.sim, ncol = n.samples)
    PS <- matrix(nrow = n.sim, ncol = n.samples)
    mean.prop <- vector(length = n.sim)
    IS <- vector(length = n.sim)
    
    for(i in 1:n.sim){
      mat[i,] <- rbinom(n = n.samples, size = n.events[j], prob = p.vec[j])
      prop.mat <- mat/n.events[j]
      mean.prop[i] <- mean(prop.mat[i,])
      PS[i,] <- 1 - abs(prop.mat[i,] - mean.prop[i])
      IS[i] <- mean(PS[i,])
    }
    
    
    quant.data[j, ] <- c(p.vec[j], quantile(IS, probs = c(0.025, 0.5, 0.975)))
  }
  quant.data
}

GenCI_IS(p.vec = est.p, n.samples = 5, n.events = n, n.sim = 10000)


