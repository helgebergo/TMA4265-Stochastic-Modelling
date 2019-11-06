# Code for the Project in Stochastic Modelling
# Lars Sæle & Helge Bergo
# November 2019

## Problem 2 - Calibrating Climate Models

# Given values
theta <- seq(0.25, 0.50, by=0.005) # parameter grid of theta
mu <- 0.5 # mean vector to n
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
Sigma_A  <- sigma * (1 + phi_m*H_A) %*% exp(-phi_m*H_A)
Sigma_B  <- sigma * (1 + phi_m*H_B) %*% exp(-phi_m*H_B)
Sigma_AB <- sigma * (1 + phi_m*H_AB) %*% exp(-phi_m*H_AB)


#factorize SIGMA = LL'
#draw n std norm z 
#return x = mu + L*z


