# Code for the Project in Stochastic Modelling
# Lars Sæle & Helge Bergo
# November 2019

## Problem 1 - Modelling the common cold


# c.1)

# Rates
lambda <- 1/100
mu <- 1/7

# Times
years <- 5
days <- 365
time <- years * days

# Initial conditions
x <- vector('numeric',time) 
x[1] <- 1 # S = 1, I = 2
s <- vector('numeric',time)
#s = c(0) # Soujourn times
totalTime = vector('numeric',time)

# Simulation
for(i in 1:time){
  state <- x[i]
  if(state == 1){
    s[i] = rexp(n = 1, rate = lambda)
    #s = c(s, tail(s,1)+rexp(n = 1, rate = lambda))
    x[i+1] <- 2
    totalTime[i+1] = totalTime[i] + s[i]
  }else{
    s[i] = rexp(n = 1, rate = mu)
    #s = c(s, tail(s,1) + rexp(n = 1, rate = mu))
    x[i+1] <- 1
    totalTime[i+1] = totalTime[i] + s[i]
  }
}


plot(NULL, NULL, xlim = c(0, time), ylim = c(0.8, 2.2), xlab = "Time (days)", 
     ylab = "State", cex.lab = 1.5, cex.axis = 1.5)
for(i in 1:(length(x)-1)){
  lines(totalTime[i:(i+1)], rep(x[i], 2), lwd = 3)
}




# ----------------------------

# c.2)


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
coldTime <- vector('numeric',time)
coldT <- 0
totTime <- 0

# Simulation
for(i in 1:time){
  state <- x[i]
  if(state == 1){
    s[i] = rexp(n = 1, rate = lambda)
    #s = c(s, tail(s,1)+rexp(n = 1, rate = lambda))
    x[i+1] <- 2
    totalTime[i+1] = totalTime[i] + s[i]
    totTime = totTime + s[i]
  }else{
    s[i] = rexp(n = 1, rate = mu)
    #s = c(s, tail(s,1) + rexp(n = 1, rate = mu))
    x[i+1] <- 1
    totalTime[i+1] = totalTime[i] + s[i]
    coldTime[i] = s[i]
    coldT = coldT + s[i]
    totTime = totTime + s[i]
  }
}
 
#coldFraction = mean(coldTime)
coldFraction = sum(coldTime)/1000/365*2
coldFraction2 = mean(coldTime)*2
leif = coldT/tail(totalTime,1)


