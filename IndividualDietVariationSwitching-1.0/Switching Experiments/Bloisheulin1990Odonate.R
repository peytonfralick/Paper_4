### get data from Bloisheulin1990Odonates

### set working directory

setwd('/media/kyle/KEC_External/Documents/Data/SwitchingAndIndVar/DataFromGraphs/')

### load libraries

library(dplyr)

### load data

Odonate.data <- read.csv('Bloisheulin1990Odonate.csv')

colnames(Odonate.data)[1] <- 'Group'

### split into two groups

A.data <- Odonate.data %>% filter(Group == 'A')

B.data <- Odonate.data %>% filter(Group == 'B')

### calculate variance for each observation
# reported as confidence interval -- CI = 1.96*SE
# SE = sqrt(variance)/sqrt(n)
# variance = (CI/2 * sqrt(n))^2

# A

A.data$Variance <- ((A.data$CI/2)*sqrt(10))^2

# B

B.data$Variance <- ((B.data$CI/2)*sqrt(10))^2

# what are the relative preferences?

# A = 0.52, B = 0.68

### fit spline to find maximum variance

# A

A.smooth <- smooth.spline(x = A.data$RelativeDensity, y = A.data$Variance)

plot(x = A.data$RelativeDensity, y = A.data$Variance)
lines(x = A.data$RelativeDensity, y = A.smooth$y)

A.predict <- predict(A.smooth, x = seq(0, 1, by = 0.01))

plot(x = seq(0, 1, by = 0.01), y = A.predict$y)

A.predict$x[which(A.predict$y == max(A.predict$y))]

# B

B.smooth <- smooth.spline(x = B.data$RelativeDensity, y = B.data$Variance)

plot(x = B.data$RelativeDensity, y = B.data$Variance)
lines(B.data$RelativeDensity, y = B.smooth$y)

B.predict <- predict(B.smooth, x = seq(0, 0.75, by = 0.01))

plot(x = seq(0, 0.75, by = 0.01), y = B.predict$y)

B.predict$x[which(B.predict$y == max(B.predict$y))]

### try to estimate magnitude of variance

### calculate p's for A and B

calcp <- function(rel.pref, rel.density) {
  rel.pref*rel.density/(rel.pref*rel.density + ((1 - rel.pref) * (1 - rel.density)))
}

# A

A.p <- calcp(rel.pref = 0.52, rel.density = A.data$RelativeDensity)

# B

B.p <- calcp(rel.pref = 0.68, rel.density = B.data$RelativeDensity)

### calculate variance -- var = p(1-p)/n. n = 8.

A.est.var <- A.p*(1 - A.p)/8

# examine relationship between observed and expected variances

plot(x = A.est.var, y = A.data$Variance)
abline(a = 0, b = 1)

# variance for group b

B.est.var <- B.p*(1 - B.p)/8

plot(x = B.est.var, y = B.data$Variance)
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

# Group A

GenCI(p.vec = A.p, n.samples = 10, n.events = 8, n.sim = 10000)

#Group B

GenCI(p.vec = B.p, n.samples = 10, n.events = 8, n.sim = 10000)

















































