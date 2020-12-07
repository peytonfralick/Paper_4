### get data from Dinis2016Beetles

### set working directory

setwd('/media/kyle/KEC_External/Documents/Data/SwitchingAndIndVar/DataFromGraphs/')

### load libraries

library(dplyr)

### load data

beetle.data <- read.csv('Dinis2016Beetles.csv')

colnames(beetle.data)[1] <- 'Species'

### seperate data into different species

Calathus.data <- beetle.data %>% filter(Species == 'Calathus')

Ptero.data <- beetle.data %>% filter(Species == 'Pterosichus')

# Calathus relative preference = 0.97

# Ptero relative preference = 0.4

### Find where maximum variance occurs

# convert standard errors to variances
# variance = (SE * sqrt(n))^2. n in this study is 25

# Calathus

Calathus.data$Variance <- (Calathus.data$StandardError * sqrt(25))^2

# Ptero

Ptero.data$Variance <- (Ptero.data$StandardError * sqrt(25))^2

### now want to find where maximum occurs for each

Calathus.smooth <- smooth.spline(x = Calathus.data$RelativeDensity, y = Calathus.data$Variance)

plot(x = Calathus.data$RelativeDensity, y = Calathus.data$Variance)
lines(x = Calathus.data$RelativeDensity, y = Calathus.smooth$y)

Calathus.predict <- predict(Calathus.smooth, x = seq(0,1,by = 0.01))

plot(x = seq(0,1, by = 0.01), y = Calathus.predict$y)

Calathus.predict$x[which(Calathus.predict$y == max(Calathus.predict$y))]

### predicts that maximum variance occurs when relative density is zero. Instead
### use the lowest observed relative abundance

### Ptero

Ptero.smooth <- smooth.spline(x = Ptero.data$RelativeDensity, y = Ptero.data$Variance)

plot(x = Ptero.data$RelativeDensity, y = Ptero.data$Variance)
lines(x = Ptero.data$RelativeDensity, y = Ptero.smooth$y)

Ptero.predict <- predict(Ptero.smooth, x = seq(0,1, by = 0.01))

plot(x = seq(0,1, by = 0.01), y = Ptero.predict$y)

Ptero.predict$x[which(Ptero.predict$y == max(Ptero.predict$y))]

##### now calculate expected variances 

### Calathus

# first calculate p

calcp <- function(rel.pref, rel.density) {
  rel.pref*rel.density/(rel.pref*rel.density + ((1 - rel.pref) * (1 - rel.density)))
}

Calathus.p <- calcp(rel.pref = 0.97, rel.density = Calathus.data$RelativeDensity)

# calculate variance p(1-p)/n where n here is from reported mean number of prey eaten

Calathus.est.var <- Calathus.p*(1-Calathus.p)/4

plot(x = Calathus.est.var, y = Calathus.data$Variance)
abline(a = 0, b = 1)

### Ptero

Ptero.p <- calcp(rel.pref = 0.49, rel.density = Ptero.data$RelativeDensity)

# Calculate variance. n here is 12

Ptero.est.var <- Ptero.p*(1-Ptero.p)/12

plot(x = Ptero.est.var, y = Ptero.data$Variance)
abline(a = 0, b = 1)

### Get confidence intervals for each of the variance estimates

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

### Calathus 

GenCI(p.vec = Calathus.p, n.samples = 25, n.events = 3, n.sim = 10000)

### Ptero

GenCI(p.vec = Ptero.p, n.samples = 25, n.events = 12, n.sim = 10000)
















