# Code for the Project in Stochastic Modelling
# Lars Sæle & Helge Bergo
# November 2019

## Problem 2 - Calibrating Climate Models
#rm(list = ls()) # clear environment and variables

# Given values
theta <- seq(0.25, 0.50, by=0.005) 
mu <- 0.5 # mean vector to n
sigma <- 0.5 # variance squared
phi_m <- 15 # Matern type correlation function parameter

# Conditional values 
theta_conditional = c(0.3, 0.35, 0.39, 0.41, 0.45) # mu_B
Y_conditional = c(0.5, 0.32, 0.40, 0.35, 0.60) # x_B
# theta_conditional = c(0.3, 0.35, 0.39, 0.41, 0.45, 0.33) # mu_B with extra conditional parameter
# Y_conditional = c(0.5, 0.32, 0.40, 0.35, 0.60, 0.40) # x_B
#theta <- setdiff(theta, theta_conditional); # remove the steps where we have a conditional value

# Build distance matrices H
ones <- as.matrix(rep(1.0, length(theta)));
H_A <- abs(theta %*% t(ones) - ones %*% t(theta)); 

ones <- as.matrix(rep(1.0, length(theta_conditional)));
H_B <- abs(theta_conditional %*% t(ones) - ones %*% t(theta_conditional)); 

H_AB <- matrix(0, nrow = length(theta), ncol = length(theta_conditional))
for(i in 1:length(theta)){
  for(j in 1:length(theta_conditional)){
    H_AB[i,j] <- abs(theta[i]-theta_conditional[j])
  }
}


# Build covariance matrices Sigma
Sigma_A  <- sigma * (1 + phi_m * H_A)  * exp(-phi_m * H_A)
Sigma_B  <- sigma * (1 + phi_m * H_B)  * exp(-phi_m * H_B)
Sigma_AB <- sigma * (1 + phi_m * H_AB) * exp(-phi_m * H_AB)



# Calculate the expected value
Y <- mu + Sigma_AB %*% solve(Sigma_B) %*% (Y_conditional - theta_conditional)


# Calculating 90% confidence interval
Var = Sigma_A - Sigma_AB %*% solve(Sigma_B) %*% t(Sigma_AB)
z = 1.64
upper <- Y + z*sqrt(diag(Var))
lower <- Y - z*sqrt(diag(Var))

# Plotting
require(plotrix) # Solution with plotrix
plotCI(theta, Y, ui=upper, li=lower)

 
# df <- data.frame(theta,Y,upper, lower)
# require(ggplot2) # Solution with ggplot
# ggplot(df, aes(x = theta, y = Y)) +
#   geom_point(size = 2) +
#   geom_errorbar(aes(ymax = upper, ymin = lower))


# Plot with error and conditional values
plot(NULL,NULL, xlim = c(0.25,0.5), ylim = c(0.1, 1.5), 
     xlab = "Theta", ylab = "y(theta)", cex.lab = 1.5)
lines(theta,Y,col="black")
lines(theta,upper,lty=2)
lines(theta,lower,lty=2)
points(theta_conditional,Y_conditional,col = "red", pch = 19)
legend(0.35,1.2,legend = c("Model","Bounds", "Reality"),
       col = c("black","black","red"),cex = 0.8, lty = c(1,2,NA), pch = c(NA,NA,19))



## 2.b)

# Calculating the treshold 
#z <- vector('numeric',51)
treshold = 0.3
p <- pnorm(treshold, mean = Y, sd = diag(Var), lower.tail = TRUE)

plot(NULL, NULL, xlim = c(0.25,0.50), ylim = c(0, 0.002), xlab = "Theta", 
     ylab = ('p(x[A] < 0.3 | x[B])'), cex.lab = 1.5, cex.axis = 1.5)
points(theta,p)

## 2.c)



# -----------

### Test plotting

# # Test to see how much we miss the conditional values
# plot(NULL,NULL, xlim = c(0.25,0.5), ylim = c(0.3, 0.8), 
#      xlab = "Theta", ylab = "y(theta)", cex.lab = 1.5)
# lines(theta,Y,col="blue")
# points(theta_conditional,Y_conditional,col = "red")
# legend(0.35,0.8,legend = c("Model", "Reality"),
#        col = c("blue","red"),cex = 0.8,lty=1:3)


# # 
# # # H_AB
# plot(NULL, NULL, xlim = c(0,50), ylim = c(0, 0.2), 
# xlab = "Index", ylab = "Correlation", cex.lab = 1.5, cex.axis = 1.5)
# for(i in 1:length(theta_conditional)){
#   lines(H_AB[,i])
# }
# 
# # H_A
# plot(NULL, NULL, xlim = c(0,50), ylim = c(0, 0.2), 
# xlab = "Index", ylab = "Correlation", cex.lab = 1.5, cex.axis = 1.5)
# for(i in 1:length(theta)){
#   lines(H_A[,i])
# }
# # H_B
# plot(NULL, NULL, xlim = c(0,50), ylim = c(0, 0.2), 
# xlab = "Index", ylab = "Correlation", cex.lab = 1.5, cex.axis = 1.5)
# for(i in 1:length(theta_conditional)){
#   lines(H_B[,i])
# }
# 
# # Sigma_AB
# plot(NULL, NULL, xlim = c(0,50), ylim = c(0.05, 0.25), 
# xlab = "Index", ylab = "Correlation", cex.lab = 1.5, cex.axis = 1.5)
# for(i in 1:length(theta_conditional)){
#   lines(Sigma_AB[,i])
# }

