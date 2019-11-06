# Code for the Project in Stochastic Modelling
# Lars Sæle & Helge Bergo
# September - October 2019

## Problem 1 - Modelling the outbreak of measles

## 1. d)
beta <- 0.05
gamma <- 0.20
N <- 1000

S <- vector('numeric', N)
R <- vector('numeric', N)

for(i in 1:N){
  S[i] <- 1
  R[i] <- 1
  
  p <- runif(1)  
  while(p > beta){
    S[i] <- S[i] + 1
    p <- runif(1)
  }
  
  p <- runif(1)  
  while(p > gamma){
    R[i] <- R[i] + 1
    p <- runif(1)
  }
}

prob_S <- mean(S)
prob_R <- mean(R)

print(prob_S)
print(prob_R)


# --------------------------

## 1. e)

pop <- 1000; # total population
n <- 200; # number of time steps
gamma <- 0.2;
S_n <- vector('numeric', n);
I_n <- vector('numeric', n);
R_n <- vector('numeric', n);
S_n[1] <- 950;
I_n[1] <- 50;


for(i in 2:n){
    beta <- 0.5*I_n[i-1]/pop
    I <- rbinom(1,S_n[i-1], beta)
    R <- rbinom(1,I_n[i-1], gamma)
    S_n[i] <- S_n[i-1] - I;
    I_n[i] <- I_n[i-1] + I - R;
    R_n[i] <- R_n[i-1] + R; 
}
time <- c(1:n);

plot(time,S_n,type='l',col='blue',lwd = 2,
     main='Population size over time',xlab='Time',
     ylab='Population size')
lines(time,I_n,type='l',col='red',lwd = 2)
lines(time,R_n,type='l',col='green',lwd = 2)
legend(x=0.8*n, y=0.7*pop, legend=c("S","I","R"), 
       col=c('blue','red','green'),lty=1, cex=1.2,lwd = 2)
      
# --------------------------

## 1. f)

pop <- 1000; # total population
n <- 200; # number of time steps
N <- 1000;
gamma <- 0.2;
S_n <- vector('numeric', n);
I_n <- vector('numeric', n);
R_n <- vector('numeric', n);
S_n[1] <- 950;
I_n[1] <- 50;

I_max <- vector('numeric', n);
time_max <- vector('numeric', n);

for(k in 1:N){
  
  for(i in 2:n){
    beta <- 0.5*I_n[i-1]/pop
    I <- rbinom(1,S_n[i-1], beta)
    R <- rbinom(1,I_n[i-1], gamma)
    S_n[i] <- S_n[i-1] - I;
    I_n[i] <- I_n[i-1] + I - R;
    R_n[i] <- R_n[i-1] + R; 
    
  }
  
  I_max[k] <- max(I_n);
  time_max[k] <- which.max(I_n);
  
}

print(mean(I_max))
print(mean(time_max))

