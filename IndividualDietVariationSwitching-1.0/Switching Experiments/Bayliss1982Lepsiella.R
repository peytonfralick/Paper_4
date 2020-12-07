### Getting data from Bayliss on Lepsiella

### set working directory

setwd('/media/kyle/KEC_External/Documents/Data/SwitchingAndIndVar/DataFromGraphs/')

### load libraries

library(dplyr)

### load data

whelk.data <- read.csv('Bayliss1982Lepsiella.csv')

### first entry weird again

colnames(whelk.data)[1] <- 'Trained'

### calculate relative preferences

### for training data sets, no c value is given. Will have to do it myself.

whelk.data <- whelk.data %>% mutate(ProportionBalanus = BalanusEaten/(BalanusEaten + ElminiusEaten))

### add a column of numeric ratios

whelk.data$NumRatio <- ifelse(whelk.data$Ratio == '30/6', 30/36, ifelse(whelk.data$Ratio == '18/18', 0.5, 6/36))

### find relative preferences using observations at 1:1 density

#Balanus

Balanus.data <- whelk.data %>% filter(Trained == 'Balanus')

Balanus.data %>% group_by(Ratio) %>% summarise(MeanProportion = mean(ProportionBalanus))

# rel.pref = 0.91

# Control

Control.data <- whelk.data %>% filter(Trained == 'Control')

Control.data %>% group_by((Ratio)) %>% summarise(MeanProportion = mean(ProportionBalanus))

# rel.pref = 0.56

# Elminius

Elminius.data <- whelk.data %>% filter(Trained == 'Elminius')

Elminius.data %>% group_by(Ratio) %>% summarise(MeanProportion = mean(ProportionBalanus)) 

# rel.pref = 0.5

### Calculate variances for each of the experiments

Variance <- function(x) {
  SS <- sum((x - mean(x))^2)
  Var <- SS/length(x)  
  Var
}

# balanus trained

Balanus.var <- Balanus.data %>% group_by(Ratio) %>% summarise(Variance = Variance(ProportionBalanus),  NumRatio = mean(NumRatio))

# control

Control.var <- Control.data %>% group_by(Ratio) %>% summarise(Variance = Variance(ProportionBalanus), NumRatio = mean(NumRatio))

# elminius trained

Elminius.var <- Elminius.data %>% group_by(Ratio) %>% summarise(Variance = Variance(ProportionBalanus), NumRatio = mean(NumRatio))

### Now want to estimate the magnitudes of variances -- p(1-p)/n

# write function to calculate p's for each relative abundance

calcp <- function(rel.pref, rel.density) {
  rel.pref*rel.density/(rel.pref*rel.density + ((1 - rel.pref) * (1 - rel.density)))
}

### Balanus variances

Balanus.est.p <- calcp(rel.pref = 0.91, rel.density = Balanus.var$NumRatio)

# calculate n's as mean total number of prey eaten at each relative abundance

Balanus.n <- Balanus.data %>% mutate(TotalFeed = BalanusEaten + ElminiusEaten) %>% group_by(Ratio) %>% summarise(n = mean(TotalFeed))

Balanus.pred.var <- Balanus.est.p*(1 - Balanus.est.p)/Balanus.n$n

# examine relationship between predicted and observed variance

plot(x = Balanus.pred.var, y = Balanus.var$Variance)
abline(a = 0, b = 1)

### Control variances

Control.est.p <- calcp(rel.pref = 0.56, rel.density = Control.var$NumRatio)

# calculate mean n for each relative abundance

Control.n <- Control.data %>% mutate(TotalFeed = BalanusEaten + ElminiusEaten) %>% group_by(Ratio) %>% summarise(n = mean(TotalFeed))

Control.pred.var <- Control.est.p * (1 - Control.est.p)/Control.n$n

# examine relationship between predicted and observed variances

plot(x = Control.pred.var, y = Control.var$Variance)
abline(a = 0, b = 1)

### Elminius variances

Elminius.est.p <- calcp(rel.pref = 0.5, rel.density = Elminius.var$NumRatio)

# calculate mean n for each relative abundance

Elminius.n <- Elminius.data %>% mutate(TotalFeed = BalanusEaten + ElminiusEaten) %>% group_by(Ratio) %>% summarise(n = mean(TotalFeed)) 

Elminius.pred.var <- Elminius.est.p*(1-Elminius.est.p)/Elminius.n$n

# examine relationship between predicted and observed variances

plot(x = Elminius.pred.var, y = Elminius.var$Variance)
abline(a = 0, b = 1)

### generate confidence intervals for predicted variance

# function to calculate 95% confidence intervals

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

# Balanus

GenCI_EventNsVector(p.vec = Balanus.est.p, n.samples = 6, n.events = round(Balanus.n$n), n.sim = 10000)

# Control

GenCI_EventNsVector(p.vec = Control.est.p, n.samples = 6, n.events = round(Control.n$n), n.sim = 10000)

#Elminius

GenCI_EventNsVector(p.vec = Elminius.est.p, n.samples = 6, n.events = round(Elminius.n$n), n.sim = 10000)

### estimate IS and and create predicted 95% CI's for IS

### calculate IS for each relative density

# Balanus

Balanus.data.PS <- Balanus.data %>% mutate(AltProp = 1 - ProportionBalanus)

Balanus.data.PS <- Balanus.data.PS %>% group_by(Ratio) %>% mutate(Diff.Prop = abs(ProportionBalanus - mean(ProportionBalanus)), 
                                                                  Diff.Alt = abs(AltProp - mean(AltProp)))

Balanus.data.PS <- Balanus.data.PS %>% mutate(PS = 1 - Diff.Prop)

Balanus.data.IS <- Balanus.data.PS %>% group_by(Ratio) %>% summarise(IS = mean(PS))

# modify above confidence interval function to get CI's for IS

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

### generate CI's for balanus

GenCI_IS(p.vec = Balanus.est.p, n.samples = 6, n.events = round(Balanus.n$n), n.sim = 10000)

# Control

# calculate IS for control

Control.data.PS <- Control.data %>% mutate(AltProp = 1 - ProportionBalanus)

Control.data.PS <- Control.data.PS %>% group_by(Ratio) %>% mutate(Diff.Prop = abs(ProportionBalanus - mean(ProportionBalanus)), 
                                                                  Diff.Alt = abs(AltProp - mean(AltProp)))

Control.data.PS <- Control.data.PS %>% mutate(PS = 1 - Diff.Prop)

Control.data.IS <- Control.data.PS %>% group_by(Ratio) %>% summarise(IS = mean(PS))

### generate CI's for control

GenCI_IS(p.vec = Control.est.p, n.samples = 6, n.events = round(Control.n$n), n.sim = 10000)

# Elminius

Elminius.data.PS <- Elminius.data %>% mutate(AltProp = 1 - ProportionBalanus)

Elminius.data.PS <- Elminius.data.PS %>% group_by(Ratio) %>% mutate(Diff.Prop = abs(ProportionBalanus - mean(ProportionBalanus)), 
                                                                    Diff.Alt = abs(AltProp - mean(AltProp)))

Elminius.data.PS <- Elminius.data.PS %>% mutate(PS = 1 - Diff.Prop)

Elminius.data.IS <- Elminius.data.PS %>% group_by(Ratio) %>% summarise(IS = mean(PS))

### generate CI's for Elminius

GenCI_IS(p.vec = Elminius.est.p, n.samples = 6, n.events = round(Elminius.n$n), n.sim = 10000)











