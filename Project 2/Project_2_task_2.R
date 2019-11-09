# Code for the Project in Stochastic Modelling
# Lars Sæle & Helge Bergo
# November 2019

## Problem 2 - Calibrating Climate Models
#rm(list = ls()) # clear environment and variables
#library(Matrix)

# Given values
theta <- seq(0.25, 0.50, by=0.005) # mu_A
sigma <- 0.5^2 # variance squared
phi_m <- 15 # Matern type correlation function parameter

# Conditional values 
theta_conditional = c(0.3, 0.35, 0.39, 0.41, 0.45) # mu_B
Y_conditional = c(0.5, 0.32, 0.40, 0.35, 0.60) # x_B

# Extra measurements from c):
theta_conditional <- c(theta_conditional, 0.33) # mu_B with extra conditional parameter
Y_conditional <- c(Y_conditional, 0.40) # x_B
#theta <- setdiff(theta, theta_conditional); # remove the steps where we have a conditional value


# Build distance matrices H
ones_A <- as.matrix(rep(1.0, length(theta)));
H_A <- abs(theta %*% t(ones_A) - ones_A %*% t(theta)); 

ones_B <- as.matrix(rep(1.0, length(theta_conditional)));
H_B <- abs(theta_conditional %*% t(ones_B) - ones_B %*% t(theta_conditional)); 

H_AB <- matrix(0, nrow = length(theta), ncol = length(theta_conditional))
H_AB <- abs(theta %*% t(ones_B) - ones_A %*% t(theta_conditional)); 


# Build covariance matrices Sigma
Sigma_A  <- sigma * (1 + phi_m * H_A)  * exp(-phi_m * H_A)
Sigma_B  <- sigma * (1 + phi_m * H_B)  * exp(-phi_m * H_B)
Sigma_AB <- sigma * (1 + phi_m * H_AB) * exp(-phi_m * H_AB)


# Calculate the expected value
Y <- theta + Sigma_AB %*% solve(Sigma_B) %*% (Y_conditional - theta_conditional) 


# Calculating 90% confidence interval
Var <- Sigma_A - Sigma_AB %*% solve(Sigma_B) %*% t(Sigma_AB)
diag(Var)[diag(Var) < 0] <- 0 # Remove negative values from floating point errors
z <- 1.64
upper <- Y + z*sqrt(diag(Var))
lower <- Y - z*sqrt(diag(Var))


# Plot with error and conditional values
#pdf("plot5.pdf") 
plot(NULL,NULL, xlim = c(0.25,0.5), ylim = c(0.2, 1.0), 
     xlab = expression(paste(theta)), ylab = expression(paste(E(Y(theta)))), cex.lab = 1.5)
par(mar=c(5,6,4,1)+.1)
lines(theta,Y,col="black")
lines(theta,upper,lty=2)
lines(theta,lower,lty=2)
points(theta_conditional,Y_conditional,col = "red", pch = 19)
legend(0.35,0.9,legend = c("Model","Bounds", "Measurements"),
       col = c("black","black","red"),cex = 0.8, lty = c(1,2,NA), pch = c(NA,NA,19))
#dev.off()


## 2.b)

# Calculating the treshold 
treshold = 0.30
p <- pnorm((treshold - Y)/sqrt(diag(Var)), lower.tail = TRUE)


#pdf("treshold5.pdf") 
plot(NULL, NULL, xlim = c(0.25,0.50), ylim = c(0.0, 0.25), xlab = expression(paste(theta)),
     ylab = (expression(paste("P(x"[A],"<0.3|x"[B],")"))), cex.lab = 1.5, cex.axis = 1.5)
par(mar=c(5,6,4,1)+.1)
lines(theta,p)
#dev.off()

best_guess <- theta[match(max(p[5:50]),p)]
print("Max value: ")
print(best_guess)


