## Code_2_Appendix_B.R
N <- 1000
lambda <- 2
# using R
Pois.R <- rpois(N,lambda)

# USING RJAGS
library(rjags)
# Model
model <- textConnection("model {
      Y ~ dpois(lambda);
}")
# Init
jags.inits <- NULL
# Data
jags.data <- list("lambda"=lambda)
# Run
model <- jags.model(model, data=jags.data, inits=jags.inits, n.chains=1)
update(model, n.iter=N)
out <- do.call(rbind.data.frame, samples)
Pois.J <- out$Y
# Compare R and JAGS results
par(mfrow=c(1,2))
hist(Pois.R,breaks=max(Pois.R),ylim=c(0,500),col="red",main="",xlab="Poisson")
summary(Pois.R)
hist(Pois.J,breaks=max(Pois.J),col="blue",xlab="Poisson",ylim=c(0,500),main="")