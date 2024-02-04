## Code_Append_C_3.R
model {
  for (k in 1:N) {
    Y[k] ~ dnorm(mu[k], tau)
    mu[k] <- alpha + beta*(x[k] - x.bar)
  }
  x.bar <- mean(x)
  alpha ~ dnorm(0.0, 1.0E-4)
  beta ~ dnorm(0.0, 1.0E-4) 
  sigma <- 1.0/sqrt(tau)
  tau ~ dgamma(1.0E-3, 1.0E-3)
}
