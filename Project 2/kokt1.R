# Time period
days = 365; # Ignoring leap years
years = 5;

# Rates
lambda = 1/100;
mu = 1/7;
rates = c(lambda, mu);

# State (1: Susceptible, 2: Infected)
xVals = c(1); # Start as susceptible
tTimes = c(0);

for (i in 1:(days*years)) {
  currState = xVals[i];
  sjTime = rexp(1, rate = rates[currState]);
  nextState = currState %% 2 + 1; # Switch state
  
  xVals = c(xVals, nextState);
  tTimes = c(tTimes, tail(tTimes,1)+sjTime);
}

# Plotting

plot(NULL, NULL, xlim = c(0, days*years), ylim = c(0.8, 2.2), xlab = "Time (days)", ylab = "State", cex.lab = 1.5, cex.axis = 1.5)
for(i in 1:(length(xVals)-1)){
  lines(tTimes[i:(i+1)], rep(xVals[i], 2), lwd = 4)
}
lines(tail(tTimes, 1) + c(0,1), c(1,1), lwd = 4)



## c.2)

# Time period
years = 1000;

# State
xVals = c(1);
totTime = 0;      # Total time the simulation runs
infectedTime = 0; # Time the individual is infected out of the total time

for (i in 1:(days*years)) {
  currState = xVals[i];
  sjTime = rexp(1, rate = rates[currState]);
  nextState = currState %% 2 + 1; # Switch state
  
  xVals = c(xVals, nextState);
  
  if (currState == 2) { # If infected; update infected time
    infectedTime = infectedTime + sjTime;
  }
  
  totTime = totTime + sjTime;
}

fracInfected = infectedTime / totTime;
print(fracInfected)

