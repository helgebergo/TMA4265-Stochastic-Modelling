### Calculate probability directly
lambda = 1.5;
days = 59;
mu = lambda*days;
claims = 100;

prob = ppois(claims, mu, lower.tail = FALSE)

### Simulate realizations and estimate probability/expectation/variance in 2a/2b
N = 1000;                  #Repeat simulation N times to get a representative estimate
numSim = 1000;             #Number of realizations
probs = numeric(N);        #Vector of estimated probabilites for each simulation in 2a)
expectations = numeric(N); #Vector of estimated expected values for each simulation in 2b)
variances = numeric(N);    #Vector of estimated variances for each simulation in 2b)
beta = 10;

for (i in 1:N){
  sim = rpois(numSim, mu);
  probs[i] = sum(sim > claims)/numSim;
  expectations[i] = sim/beta;     #Formula for expectation after law of total expectation
  variances[i] = 2*sim/(beta^2);  #Formula for variance after law of total variance
}

cat(sprintf("Probability std: %f\n", sd(probs)));
cat(sprintf("Probability mean: %f\n", mean(probs)));

cat(sprintf("Expectation std: %f\n", sd(expectations)));
cat(sprintf("Expectation mean: %f\n", mean(expectations)));

cat(sprintf("Variance std: %f\n", sd(variances)));
cat(sprintf("Variance mean: %f\n", mean(variances)));

### Figure with 10 realizations of X(t)
library(ggplot2) #Must install package "ggplot2" to plot

t = 59;
realizations = 10;

plot = ggplot(); #Initialize plot to be iteratively updated

for (realization in 1:realizations) {
  n = rpois(1, lambda*t);
  u = numeric(n);
  for (i in 1:n) {
    u[i] = runif(1,0,t);
  }
  w = sort(u);
  
  w = c(0, w, t);         #Add start and end time for plotting purposes
  events = c(1, 1:n, n);  #Duplicate first and last event for plotting purposes
  
  df = data.frame(w, events); #Create dataframe to use ggplot
  plot = plot + geom_line(aes_string(x = w, y = events, color= shQuote(realization)));
}

plot + labs(x = "Time", y = "Events") + theme(legend.position = "none")

