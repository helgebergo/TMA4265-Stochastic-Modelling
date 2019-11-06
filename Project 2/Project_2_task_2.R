# Code for the Project in Stochastic Modelling
# Lars Sæle & Helge Bergo
# November 2019

## Problem 2 - Calibrating Climate Models

# Given values
theta <- seq(0.25, 0.50, by=0.005) # parameter grid of theta
mu <- 0.5 # mean vector to n
mu_conditional <- 0
sigma <- 0.5^2 # variance squared
phi_m <- 15 # correlation function parameter

# Conditional values 
theta_conditional = c(0.3, 0.35, 0.39, 0.41, 0.45)
times_conditional = c(0.5, 0.32, 0.40, 0.35, 0.60)

# Build distance matrices H
ones <- as.matrix(rep(1.0, length(theta)));
H_A <- abs(theta %*% t(ones) - ones %*% t(theta)); 

ones <- as.matrix(rep(1.0, length(theta_conditional)));
H_B <- abs(theta_conditional %*% t(ones) - ones %*% t(theta_conditional)); 

H_AB <- matrix(0, nrow = 51, ncol = 5)
for(i in 1:51){
  for(j in 1:5){
    H_AB[i,j] <- abs(theta[i]-theta_conditional[j])
  }
}


# Build covariance matrices Sigma
Sigma_A  <- sigma * (1 + phi_m*H_A)  * exp(-phi_m*H_A)
Sigma_B  <- sigma * (1 + phi_m*H_B)  * exp(-phi_m*H_B)
Sigma_AB <- sigma * (1 + phi_m*H_AB) * exp(-phi_m*H_AB)


# Calculate the answer
E <- mu + Sigma_AB %*% solve(Sigma_B) %*% (theta_conditional - mu_conditional)


# Calculating 90% confidence interval
Var = Sigma_A - Sigma_AB %*% solve(Sigma_B) %*% t(Sigma_AB)
upper = c()
lower = c()
z = 1.64
for(i in 1:51){
  upper[i] <- E[i] + z*sqrt(Var[i,i])
  lower[i] <- E[i] - z*sqrt(Var[i,i])
}

require(plotrix) # Solution with plotrix
plotCI(theta, E, ui=upper, li=lower)
points(theta,E)

df <- data.frame(theta,E,upper, lower)
require(ggplot2) # Solution with ggplot
ggplot(df, aes(x = theta, y = E)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymax = upper, ymin = lower))



# Plotting 
df <- data.frame(theta,E,theta_conditional,times_conditional,upper,lower)
plot(theta,E,type='l')
par(new=TRUE)
plot(theta_conditional,times_conditional,col='red')

require(ggplot2) # Solution with ggplot
ggplot(df, aes(x = theta, y = E)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymax = upper, ymin = lower)) +
  geom_point(aes(x = theta_conditional, t = times_conditional))


