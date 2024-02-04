## Code_1_Appendix_B.R
# Use rjags library
library(rjags)
# Export JAGS model
model <- textConnection("model {
  Y ~ dbin(0.5,10)
  Pr <- step(3.1 - Y)
}")
# Init
jags.inits <- NULL
# Data
jags.data <- list()
# Perform Bayesian analysis using JAGS
model <- jags.model(model, data=jags.data, inits=jags.inits, n.chains=1)
update(model, n.iter=10000)
# Parameters
jags.params <- c("Y","Pr") # parameters to be monitored
samples <- coda.samples(model, jags.params,10000)
# Analysis of the simulation results
plot(samples)
summary(samples)
# extract P and Y
out <- do.call(rbind.data.frame, samples)
names(out)
# Probability
summary(out$Pr)
sum(out$Pr==0)/length(out$Pr)
sum(out$Pr==1)/length(out$Pr)
# Flipping results
summary(out$Y)
hist(out$Y,breaks=c(0,1,2,3,4,5,6,7,8,9,10),main="",freq=FALSE,xlab="Number of heads")
