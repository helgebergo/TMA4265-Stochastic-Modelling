# Code for the Project in Stochastic Modelling
# Lars Sæle & Helge Bergo
# September - October 2019

## Problem 2 - Insurance claims

# Simulate 1000 realizations of the Poission process

N <- 10000000
lambda <- 1.5 
t <- 59 # days
mu <- lambda * t 
limit <- 100 # value to be checked agains

claims <- numeric(N) # number of claims for each realization

beta = 10;
E <- numeric(N);     # expected total claim amount 
Var <- numeric(N);   # variance of total claim amount   

for(i in 1:N){
  claims[i] = rpois(1,mu)  # Number of claims for each realization
  E[i] = claims[i]/beta;    
  Var[i] = 2*claims[i]/(beta^2); 
}

result = mean(claims > limit)

cat(sprintf("Probability: %0.6f\n", result));
cat(sprintf("Expectation: %0.5f\n", mean(E)));
cat(sprintf("Variance: %0.5f\n", mean(Var)));



# Plot of 10 realizations
realizations <- 10;

library(ggplot2) # import ggplot 2, to use several plots in one
plot = ggplot() 

for (k in 1:realizations) {
  # For each realization, the number of claims are calculated
  claims = rpois(1,mu)
  u = numeric(claims)
  # The claims are distributed using a normal distribution from 0 to t
  for (i in 1:claims) {
    u[i] = runif(1,0,t)
  }
  
  events = c(1:claims) # the accumulated list of claims (y-axis)
  w = sort(u); # the "timing" of each claim (x-axis)
  
  df = data.frame(w, events); 
  plot = plot + geom_line(aes_string(x = w, y = events, color= shQuote(k)))
}

plot + labs(x = "Time", y = "Claims") + theme(legend.position = "none")


