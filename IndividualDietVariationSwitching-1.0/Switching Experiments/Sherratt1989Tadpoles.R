### get data from Sherratt1989Tadpoles

### set working directory

setwd('/media/kyle/KEC_External/Documents/Data/SwitchingAndIndVar/DataFromGraphs')

### load libraries

library(dplyr)

### load data

Tadpole.data <- read.csv('Sherratt1989Tadpoles.csv')

colnames(Tadpole.data)[1] <- 'NTreatment'

### split data into different experiments

ten.data <- Tadpole.data %>% filter(NTreatment == 10)

twenty.data <- Tadpole.data %>% filter(NTreatment == 20)

forty.data <- Tadpole.data %>% filter(NTreatment == 40)

### manipulate data to have one row for each density

St.Dev <- function(x) {
  SS <- sum((x - mean(x))^2)
  SD <- sqrt(SS/length(x))  
  SD
}


ten.data <- ten.data %>% group_by(RelativeDensity) %>% summarise(MeanProportion = mean(Proportion), StandardDeviation = St.Dev(Proportion))

twenty.data <- twenty.data %>% group_by(RelativeDensity) %>% summarise(MeanProportion = mean(Proportion), StandardDeviation = St.Dev(Proportion))

forty.data <- forty.data %>% group_by(RelativeDensity) %>% summarise(MeanProportion = mean(Proportion), StandardDeviation = St.Dev(Proportion))

### where maximum variances occur, transform to variance

ten.data$Variance <- (ten.data$StandardDeviation)^2

twenty.data$Variance <- (twenty.data$StandardDeviation)^2

forty.data$Variance <- (forty.data$StandardDeviation)^2

### what are the relative preferences ... can use the proportions at 0.5 -- 0.4, 0.313, 0.255

ConvertC <- function(c) {
  RelPref <- c/(1+c)
  MaxVar <- 1 - RelPref
  print(as.data.frame(cbind(RelPref, MaxVar)))
}

### or could use the estimates from linear regression

# ten prey

lm(MeanProportion ~ 0 + RelativeDensity, ten.data)

# c = 0.7799

ConvertC(c = 0.7799)

# twenty prey

lm(MeanProportion ~ 0 + RelativeDensity, twenty.data)

# c = 0.7457

ConvertC(c = 0.7457)

# forty prey

lm(MeanProportion ~ 0 + RelativeDensity, forty.data)

# c = 0.5964

ConvertC(c = 0.5964)

### try to predict the variances, need p's 

calcp <- function(rel.pref, rel.density) {
  rel.pref*rel.density/(rel.pref*rel.density + ((1 - rel.pref) * (1 - rel.density)))
}

### ten prey

ten.p <- calcp(rel.pref = 0.4, rel.density = ten.data$RelativeDensity)

### twenty prey

twenty.p <-  calcp(rel.pref = 0.31, rel.density = twenty.data$RelativeDensity)

### forty prey

forty.p <- calcp(rel.pref = 0.25, rel.density = forty.data$RelativeDensity)

### estimate variances as p(1-p)/n. N's are ~ 50

### ten prey

ten.est.var <- ten.p*(1-ten.p)/30

plot(x = ten.est.var, y = ten.data$Variance)
abline(a = 0, b = 1)

### twenty prey

twenty.est.var <- twenty.p*(1 - twenty.p)/30

plot(x = twenty.est.var, y = twenty.data$Variance)
abline(a = 0, b = 1)

### forty prey

forty.est.var <- forty.p*(1 - forty.p)/30

plot(x = forty.est.var, y = forty.data$Variance)
abline(a = 0, b =1)

### calculate confidence intervals for each of the variance estimates

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


### ten prey

GenCI(p.vec = ten.p, n.samples = 6, n.events = 30, n.sim = 10000)

### twenty prey

GenCI(p.vec = twenty.p, n.samples = 6, n.events = 30, n.sim = 10000)

### forty prey

GenCI(p.vec = forty.p, n.samples = 6, n.events = 30, n.sim = 10000)

### Estimate IS and create predicted 95% CI's for IS

### Calculate IS for each relative abundance

# ten prey

ten.data.ind <- Tadpole.data %>% filter(NTreatment == 10)

ten.data.PS <- ten.data.ind %>% mutate(PropAlt = 1 - Proportion)

ten.data.PS <- ten.data.PS %>% group_by(RelativeDensity) %>% mutate(Diff.Prop = abs(Proportion - mean(Proportion)), 
                                                                    Diff.Alt = abs(PropAlt - mean(PropAlt)))

ten.data.PS <- ten.data.PS %>% mutate(PS = 1 - Diff.Prop)

ten.data.IS <- ten.data.PS %>% group_by(RelativeDensity) %>% summarise(IS = mean(PS))

# twenty prey

twenty.data.ind <- Tadpole.data %>% filter(NTreatment == 20)

twenty.data.PS <- twenty.data.ind %>% mutate(PropAlt = 1 - Proportion)

twenty.data.PS <- twenty.data.PS %>% group_by(RelativeDensity) %>% mutate(Diff.Prop = abs(Proportion - mean(Proportion)), 
                                                                          Diff.Alt = abs(PropAlt - mean(PropAlt)))

twenty.data.PS <- twenty.data.PS %>% mutate(PS = 1 - Diff.Prop)

twenty.data.IS <- twenty.data.PS %>% group_by(RelativeDensity) %>% summarise( IS = mean(PS))

# forty prey

forty.data.ind <- Tadpole.data %>% filter(NTreatment == 40)

forty.data.PS <- forty.data.ind %>% mutate(PropAlt = 1 - Proportion)

forty.data.PS <- forty.data.PS %>% group_by(RelativeDensity) %>% mutate(Diff.Prop = abs(Proportion - mean(Proportion)), 
                                                                        Diff.Alt = abs(PropAlt - mean(PropAlt)))

forty.data.PS <- forty.data.PS %>% mutate(PS = 1 - Diff.Prop)

forty.data.IS <- forty.data.PS %>% group_by(RelativeDensity) %>% summarise(IS = mean(PS))

### modify the CI function above to generate confidence intervals for IS

GenCI_IS <- function(p.vec, n.samples, n.events, n.sim) {
  
  quant.data <- matrix(nrow = length(p.vec), ncol = 4)
  
  for(j in 1:length(p.vec)){
    mat <- matrix(nrow = n.sim, ncol = n.samples)
    PS <- matrix(nrow = n.sim, ncol = n.samples)
    mean.prop <- vector(length = n.sim)
    IS <- vector(length = n.sim)
    
    for(i in 1:n.sim){
      mat[i,] <- rbinom(n = n.samples, size = n.events, prob = p.vec[j])
      prop.mat <- mat/n.events
      mean.prop[i] <- mean(prop.mat[i])
      PS[i,] <- 1 - abs(prop.mat[i,] - mean.prop[i])
      IS[i] <- mean(PS[i,])
    }
    
    quant.data[j, ] <- c(p.vec[j], quantile(IS, probs = c(0.025, 0.5, 0.975)))
  }
  quant.data
}

# ten prey

GenCI_IS(p.vec = ten.p, n.samples = 6, n.events = 30, n.sim = 10000)

# twenty prey

GenCI_IS(p.vec = twenty.p, n.samples = 6, n.events = 30, n.sim = 10000)

# forty prey

GenCI_IS(p.vec = forty.p, n.samples = 6, n.events = 30, n.sim = 10000)
