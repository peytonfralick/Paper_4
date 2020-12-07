### get data from Vantornhout2006Mites

### set working directory

setwd('D:/Documents/Data/SwitchingAndIndVar/DataFromGraphs/')

### load libraries

library(dplyr)

### load data

mite.data <- read.csv('Vantornhout2006Mites.csv')

colnames(mite.data)[1] <- 'Experiment'

### split up data into different experiments

FvT.data <- mite.data %>% filter(Experiment == 'FvT')

FvTeggs.data <- mite.data %>% filter(Experiment == 'FvTeggs')

Instars.data <- mite.data %>% filter(Experiment == 'Instars')

TvT.data <- mite.data %>% filter(Experiment == 'TvT')

### calculate variances from Standard Errors
# Var = (SE * sqrt(n))^2
# n = 20

FvT.data$Variance <- (FvT.data$StandardError/2 * sqrt(20))^2

FvTeggs.data$Variance <- (FvTeggs.data$StandardError/2 * sqrt(20))^2

Instars.data$Variance <- (Instars.data$StandardError/2 * sqrt(20))^2

TvT.data$Variance <- (TvT.data$StandardError/2 * sqrt(20))^2

### Relative preference of mites -- 
# FvT -- 0.2
# FvTeggs -- 0.356
# Instars -- 0.84
# TvT -- 0.466

### try to estimate the magnitude of the variance for each experiment

### calculate p first

calcp <- function(rel.pref, rel.density) {
  rel.pref*rel.density/(rel.pref*rel.density + ((1 - rel.pref) * (1 - rel.density)))
}

### FvT p 

FvT.p <- calcp(rel.pref = 0.2, rel.density = FvT.data$RelativeDensity)

# FvTeggs p 

FvTeggs.p <- calcp(rel.pref = 0.356, rel.density = FvTeggs.data$RelativeDensity)

# Instars p 

Instars.p <- calcp(rel.pref = 0.84, rel.density = Instars.data$RelativeDensity)

# TvT p

TvT.p <- calcp(rel.pref = 0.466, rel.density = TvT.data$RelativeDensity)

### calculate expected variances 
# var = p(1-p)/n  need n's for each of the experiments

#FvT

FvT.n <- c(11,12,9)

# FvTeggs

FvTeggs.n <- c(13,9,10)

# Instars

Instars.n <- c(4,3,5)

# TvT

TvT.n <- c(4,5,5)

### calculate variances

# FvT

FvT.exp.var <- FvT.p*(1-FvT.p)/FvT.n

plot(x = FvT.exp.var, y = FvT.data$Variance)
abline(a = 0, b = 1)

# FvTeggs

FvTeggs.exp.var <- FvTeggs.p*(1 - FvTeggs.p)/FvTeggs.n
plot(x = FvTeggs.exp.var, y = FvTeggs.data$Variance)
abline(a = 0, b = 1)

# Instars

Instars.exp.var <- Instars.p*(1 - Instars.p)/Instars.n

plot(x = Instars.exp.var, y = Instars.data$Variance)
abline(a = 0, b = 1)

# TvT

TvT.exp.var <- TvT.p*(1 - TvT.p)/TvT.n

plot(x = TvT.exp.var, y = TvT.data$Variance)
abline(a = 0, b = 1)


### simulate confidence intervals for variance estimates

GenCI.VarEvents <- function(p.vec, n.samples, n.events, n.sim) {
  
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

### function will use data to generate a CI, so we can see whether the observed variance is 
### within the confidence interval.

# FvT

GenCI.VarEvents(p.vec = FvT.p, n.samples = 20, n.events = FvT.n, n.sim = 10000)

#FvTEggs

GenCI.VarEvents(p.vec = FvTeggs.p, n.samples = 20, n.events = FvTeggs.n, n.sim = 10000)

# Instars

GenCI.VarEvents(p.vec = Instars.p, n.samples = 20, n.events = Instars.n, n.sim = 10000)

#TvT

GenCI.VarEvents(p.vec = TvT.p, n.samples = 20, n.events = TvT.n, n.sim = 10000)

















