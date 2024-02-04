## Code_4_Appendix_B.R
# Read data
data(cars)
X <- cars$speed
X <- X-mean(X)
Y <- cars$dist
# Perform JAGS computation
N <- length(X) # number of data points

# Model
mod <- textConnection("model {
   mu ~ dnorm(0, 0.01); 
   tau.pro ~ dgamma(0.001,0.001); 
   sd.pro <- 1/sqrt(tau.pro);
   beta ~ dnorm(0,0.01)
   
   for(i in 1:N) {
      predY[i] <- mu + beta*X[i]; 
      Y[i] ~ dnorm(predY[i], tau.pro);
   }
}")
# Init and definitions
jags.data = list("Y"=Y,"N"=N,"X"=X)
jags.params=c("sd.pro","mu","beta")
jags.inits <- list("tau.pro" = 1, "mu" = 0, "beta" = 0)
# Run
model <- jags.model(mod, data=jags.data, inits=jags.inits, n.chains=1)
update(model, n.iter=10000)
samples <- coda.samples(model, jags.params,10000)
out <- do.call(rbind.data.frame, samples)
names(out)
# Comparison with R linear fitting
plot(X,Y)
abline(lm(Y~X),col="red")
abline(mean(out$mu),mean(out$beta),col="blue")
# Coefficients
coefficients(lm(Y~X))
c(mean(out$mu),mean(out$beta))