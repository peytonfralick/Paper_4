### get data from Murdoch1973Ladybirds

### set working directory

setwd('/media/kyle/KEC_External/Documents/Data/SwitchingAndIndVar/DataFromGraphs/')

### load libraries

library(dplyr)

### load data

ladybird.data <- read.csv('Murdoch1973Ladybirds.csv') 

colnames(ladybird.data)[1] <- 'Trained'

### split up data by experiment

Acyrth.data <- ladybird.data %>% filter(Trained == 'Acyrth')

Aphis.data <- ladybird.data %>% filter(Trained == 'Aphis')

Control.data <- ladybird.data %>% filter(Trained == 'Control')

### write function to calculate standard deviation 

St.Dev <- function(x) {
  SS <- sum((x - mean(x))^2)
  SD <- sqrt(SS/length(x))  
  SD
}


Acyrth.data <- Acyrth.data %>% group_by(RelativeDensity) %>%mutate(Proportion = AcyrthEaten/(AcyrthEaten + AphisEaten))%>% summarise(MeanProportion = mean(Proportion), StandardDeviation = St.Dev(Proportion), n = n())

Aphis.data <- Aphis.data %>% group_by(RelativeDensity) %>% mutate(Proportion = AcyrthEaten/(AcyrthEaten + AphisEaten)) %>% summarise(MeanProportion = mean(Proportion), StandardDeviation = St.Dev(Proportion), n = n())

Control.data <- Control.data %>% group_by(RelativeDensity) %>% mutate(Proportion = AcyrthEaten/(AcyrthEaten + AphisEaten)) %>% summarise(MeanProportion = mean(Proportion), StandardDeviation = St.Dev(Proportion), n = n())

### relative preferences 

# Acyrth is 0.38

# Aphis is 0.41

# Control is 0.41

### need to estimate where the greatest variance occurs for each experiment

# first calculate variances from standard deviations var = (sd)^2

Acyrth.data$Variance <- (Acyrth.data$StandardDeviation)^2

Aphis.data$Variance <- (Aphis.data$StandardDeviation)^2

Control.data$Variance <- (Control.data$StandardDeviation)^2

### fit splines for each of the experiments

# Acyrth

Acyrth.smooth <- smooth.spline(x = Acyrth.data$RelativeDensity, y = Acyrth.data$Variance)

plot(x = Acyrth.data$RelativeDensity, y = Acyrth.data$Variance)
lines(x = Acyrth.data$RelativeDensity, y = Acyrth.smooth$y)

Acyrth.predict <- predict(Acyrth.smooth, x = seq(0, 1, by = 0.01))

plot(x = seq(0, 1, by = 0.01), y = Acyrth.predict$y)

Acyrth.predict$x[which(Acyrth.predict$y == max(Acyrth.predict$y))]

# Aphis

Aphis.smooth <- smooth.spline(x = Aphis.data$RelativeDensity, y = Aphis.data$Variance, df = 3)

plot(x = Aphis.data$RelativeDensity, y = Aphis.data$Variance)
lines(x = Aphis.data$RelativeDensity, y = Aphis.smooth$y)

Aphis.predict <- predict(Aphis.smooth, x = seq(0, 1, by = 0.01))

plot(x = seq(0, 1, by = 0.01), y = Aphis.predict$y)

Aphis.predict$x[which(Aphis.predict$y == max(Aphis.predict$y))]

# Control

Control.smooth <- smooth.spline(x = Control.data$RelativeDensity, y = Control.data$Variance)

plot(x = Control.data$RelativeDensity, y = Control.data$Variance)
lines(x = Control.data$RelativeDensity, y = Control.smooth$y)

Control.predict <- predict(Control.smooth, x = seq(0, 1, by = 0.01))

plot(x = seq(0, 1, by = 0.01), y = Control.predict$y)

Control.predict$x[which(Control.predict$y == max(Control.predict$y))]

### want to predict the magnitude of the variances 

### need p's for each experiment

calcp <- function(rel.pref, rel.density) {
  rel.pref*rel.density/(rel.pref*rel.density + ((1 - rel.pref) * (1 - rel.density)))
}

# Acyrth

Acyrth.p <- calcp(rel.pref = 0.38, rel.density = Acyrth.data$RelativeDensity)

# Aphis

Aphis.p <- calcp(rel.pref = 0.41, rel.density = Aphis.data$RelativeDensity)

# Control

Control.p <- calcp(rel.pref = 0.41, rel.density = Control.data$RelativeDensity)

### need to find n's for each of the experiments. Have a bit more detailed information for this experiment relative to others

n.data <- read.csv('Murdoch1973Ladybirds.csv')

colnames(n.data)[1] <- 'Trained'

n.data <- n.data %>% mutate(Total = AcyrthEaten + AphisEaten)

Acyrth.n <- n.data %>% filter(Trained == 'Acyrth') %>% group_by(RelativeDensity) %>% summarise(n = mean(Total))

Aphis.n <- n.data %>% filter(Trained == 'Aphis') %>% group_by(RelativeDensity) %>% summarise(n = mean(Total))

Control.n <- n.data %>% filter(Trained == 'Control') %>% group_by(RelativeDensity) %>% summarise(n = mean(Total))

### have all of the pieces to estimate the variance

# Acyrth

Acyrth.est.var <- Acyrth.p*(1 - Acyrth.p)/Acyrth.n$n

plot(x = Acyrth.est.var, y = Acyrth.data$Variance)
abline(a = 0, b = 1)

# Aphis

Aphis.est.var <- Aphis.p*(1-Aphis.p)/Aphis.n$n

plot(x = Aphis.est.var, y = Aphis.data$Variance)
abline(a = 0, b = 1)

#Control

Control.est.var <- Control.p*(1 - Control.p)/Control.n$n

plot(x = Control.est.var, y = Control.data$Variance)
abline(a = 0, b = 1)

### calculate confidence intervals for variance estimates
### need to modify function to deal with variable 
### number of events and variable number of samples.

GenCI.VarSampleVarEvent <- function(p.vec, n.samples, n.events, n.sim) {
  
  quant.data <- matrix(nrow = length(p.vec), ncol = 3)
  
  for(j in 1:length(p.vec)){
    mat <- matrix(nrow = n.sim, ncol = n.samples[j])
    for(i in 1:n.sim){
      mat[i,] <- rbinom(n = n.samples[j], size = n.events[j], prob = p.vec[j])
    }
    prop.mat <- mat/n.events[j]
    var.vec <- apply(prop.mat, 1, var)
    quant.data[j, ] <- c(p.vec[j], quantile(var.vec, probs = c(0.025, 0.975)))
  }
  quant.data
}

### Acyrth

# get number of samples

Acyrth.n.samples <- ladybird.data %>% filter(Trained == 'Acyrth') %>% group_by(RelativeDensity) %>% summarise(n = n())

GenCI.VarSampleVarEvent(p.vec = Acyrth.p, n.samples = Acyrth.n.samples$n, n.events = round(Acyrth.n$n), n.sim = 10000)

### Aphis

Aphis.n.samples <- ladybird.data %>% filter(Trained == 'Aphis') %>% group_by(RelativeDensity) %>% summarise(n = n())

GenCI.VarSampleVarEvent(p.vec = Aphis.p, n.samples = Aphis.n.samples$n, n.events = round(Aphis.n$n), n.sim = 10000)

### Control

Control.n.samples <- ladybird.data %>% filter(Trained == 'Control') %>% group_by(RelativeDensity) %>% summarise(n = n())

GenCI.VarSampleVarEvent(p.vec = Control.p, n.samples = Control.n.samples$n, n.events = round(Control.n$n), n.sim = 10000)

### Estimate IS and create 95% CI's for IS

# Estimate IS at each of the relative abundances

# Acyrth

Acyrth.data.ind <- ladybird.data %>% filter(Trained == 'Acyrth')

Acyrth.data.ind <- Acyrth.data.ind %>% mutate(Proportion = AcyrthEaten/(AcyrthEaten + AphisEaten))

Acyrth.data.PS <- Acyrth.data.ind %>% mutate(Prop.Alt = 1 - Proportion)

Acyrth.data.PS <- Acyrth.data.PS %>% group_by(RelativeDensity) %>% mutate(Diff.Prop = abs(Proportion - mean(Proportion)), 
                                                                          Diff.Alt = abs(Prop.Alt - mean(Prop.Alt)))

Acyrth.data.PS <- Acyrth.data.PS %>% mutate(PS = 1 - Diff.Prop)

Acyrth.data.IS <- Acyrth.data.PS %>% group_by(RelativeDensity) %>% summarise(IS = mean(PS))

# Aphis

Aphis.data.ind <- ladybird.data %>% filter(Trained == 'Aphis')

Aphis.data.ind <- Aphis.data.ind %>% mutate(Proportion = AcyrthEaten/(AcyrthEaten + AphisEaten))

Aphis.data.PS <- Aphis.data.ind %>% mutate(Prop.Alt = 1 - Proportion)

Aphis.data.PS <- Aphis.data.PS %>% group_by(RelativeDensity) %>% mutate(Diff.Prop = abs(Proportion - mean(Proportion)), 
                                                                        Diff.Alt = abs(Prop.Alt - mean(Prop.Alt)))

Aphis.data.PS <- Aphis.data.PS %>% mutate(PS = 1 - Diff.Prop)

Aphis.data.IS <- Aphis.data.PS %>% group_by(RelativeDensity) %>% summarise(IS = mean(PS))

# Control

Control.data.ind <- ladybird.data %>% filter(Trained == 'Control')

Control.data.ind <- Control.data.ind %>% mutate(Proportion = AcyrthEaten/(AcyrthEaten + AphisEaten))

Control.data.PS <- Control.data.ind %>% mutate(Prop.Alt = 1 - Proportion)

Control.data.PS <- Control.data.PS %>% group_by(RelativeDensity) %>% mutate(Diff.Prop = abs(Proportion - mean(Proportion)), 
                                                                            Diff.Alt = abs(Prop.Alt - mean(Prop.Alt)))

Control.data.PS <- Control.data.PS %>% mutate(PS = 1 - Diff.Prop)

Control.data.IS <- Control.data.PS %>% group_by(RelativeDensity) %>% summarise(IS = mean(PS))

### modify CI function above to calculate 95% CI's for IS

GenCI_IS <- function(p.vec, n.samples, n.events, n.sim) {
  
  quant.data <- matrix(nrow = length(p.vec), ncol = 4)
  
  for(j in 1:length(p.vec)){
    mat <- matrix(nrow = n.sim, ncol = n.samples[j])
    PS <- matrix(nrow = n.sim, ncol = n.samples[j])
    mean.prop <- vector(length = n.sim)
    IS <- vector(length = n.sim)
    
    for(i in 1:n.sim){
      mat[i,] <- rbinom(n = n.samples[j], size = n.events[j], prob = p.vec[j])
      prop.mat <- mat/n.events[j]
      mean.prop[i] <- mean(prop.mat[i,])
      PS[i,] <- 1 - abs(prop.mat[i,] - mean.prop[i])
      IS[i] <- mean(PS[i,])
    }
    
    
    quant.data[j, ] <- c(p.vec[j], quantile(IS, probs = c(0.025, 0.5, 0.975)))
  }
  quant.data
}

### Acyrth

GenCI_IS(p.vec = Acyrth.p, n.samples = Acyrth.n.samples$n, n.events = round(Acyrth.n$n), n.sim = 10000)

### Aphis

GenCI_IS(p.vec = Aphis.p, n.samples = Aphis.n.samples$n, n.events = round(Aphis.n$n), n.sim = 10000)

### Control

GenCI_IS(p.vec = Control.p, n.samples = Control.n.samples$n, n.events = round(Control.n$n), n.sim = 10000)
