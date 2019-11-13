# Code for the Project in Stochastic Modelling
# Lars Sæle & Helge Bergo
# November 2019

## Problem 1 - Modelling the common cold

# c)


# Rates
lambda <- 1/100
mu <- 1/7

# Times
years <- 1000
days <- 365
time <- years * days

# Initial conditions
x <- vector('numeric',time) 
x[1] <- 1 # S = 1, I = 2
s <- vector('numeric',time)
totalTime = vector('numeric',time)
coldTime <- 0

# Simulation over 1000 years
for(i in 1:time){
  state <- x[i]
  if(state == 1){
    s[i] = rexp(n = 1, rate = lambda)
    x[i+1] <- 2
    totalTime[i+1] = totalTime[i] + s[i]
  }else{
    s[i] = rexp(n = 1, rate = mu)
    x[i+1] <- 1
    totalTime[i+1] = totalTime[i] + s[i]
    coldTime = coldTime + s[i]
  }
}


plot(NULL, NULL, xlim = c(0, 5*365), ylim = c(0.8, 2.2), xlab = "Time (days)", 
     ylab = "State", cex.lab = 1.5, cex.axis = 1.5)
for(i in 1:100){
  lines(totalTime[i:(i+1)], rep(x[i], 2), lwd = 3)
}

coldMeanfraction = coldTime/tail(totalTime,1)

