# This file contains code adapted from
# https://github.com/daniele-falco/software_comparison/tree/main/linear_models

##############3
## "L" (Lasso) cases

# Nimble 
library(nimble) 

#MODEL: linear model, lasso prior


# Put their data simulation code into a function
sim_lm_data_L_case <- function(N = 100, K = 4, seed = 1234, number_zero_coefs = 0) {
  set.seed(seed)
  #SIMULATION OF DATA
  # N=100 #number of observations
  # K<-16  #length of our parameters 16 or also 120 
  
  # Matrix X of covariates
  X <- matrix(nrow=N, ncol=K)
  X[,1] <- rep(1,N) #intercept
  for(i in 2:K){
    X[,i] <- rnorm(n = N, mean = 0 ,sd = 1) #covariates 
  }

  ## Beraha et al's code commented out cases for p, leaving one at a time.
  ## Here I will accomplish that using if-then-else
  if(K == 16) {
    #when p=16
    beta=c(1,-1,-2,3,2,0.5,0.05,-0.005,0.001,1,5,-0.02,2.5,1.5,-0.5,0.01)
  }
  if(K == 30) {
    #when p=30
    if(number_zero_coefs == 2)
      beta=c(1,-1,-2,3,2,0.5,0.05,-0.005,0.001,1,5,-0.02,2.5,1.5,-0.5,
             0.01,1.2,-1.5,-2.9,3.2,1.4,0.05,0.5,-0.009,0.001,2.5,4.5,-0.09,0,0)
  
    if(number_zero_coefs == 15)
      beta=c(1,-1,-2,3,2,0.5,0.05,-0.005,0.001,1,5,-0.02,2.5,1.5,-0.5, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    
    if(number_zero_coefs == 28) 
      beta=c( 2,-2,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  }
  if(K == 100) {
    #when p=100
    if(number_zero_coefs == 2)
      beta=c(1,-1,-2,3,2,0.5,0.05,-0.005,0.001,1,5,-0.02,2.5,1.5,-0.5,
             0.01,1.2,-1.5,-2.9,3.2,1.4,0.05,0.5,-0.009,0.001,2.5,4.5,-0.09,1,2,3,-1,-2,-3,0.5,1.5,2.5,-1,-1.5,-2.5,-3,0.2,0.6,0.7,-0.5,-0.9,1,-1,-2,3,2,0.5,0.05,-0.005,0.001,1,5,-0.02,2.5,1.5,-0.5,
             0.01,1.2,-1.5,-2.9,3.2,1.4,0.05,0.5,-0.009,0.001,2.5,4.5,-0.09,1,2,3,-1,-2,-3,0.5,1.5,2.5,-1,-1.5,-2.5,-3,0.2,0.6,0.7,-0.5,-0.9,-0.2,0.3,0.9,1.5,1.8,2.4,0,0)
    if(number_zero_coefs == 50)
      beta=c(1,-1,-2,3,2,0.5,0.05,-0.005,0.001,1,5,-0.02,2.5,1.5,-0.5,
             0.01,1.2,-1.5,-2.9,3.2,1.4,0.05,0.5,-0.009,0.001,2.5,4.5,-0.09,1,2,3,-1,-2,-3,0.5,1.5,2.5,-1,-1.5,-2.5,-3,0.2,0.6,0.7,-0.5,-0.9,1,-1,-2,3,
             0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0)
    if(number_zero_coefs == 98)
      beta=c( 2,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  }

  ## The next three lines are in Beraha et al's code (commented-out in their code as well as here).
  ## I am not curently seeing repeated simulations reported in their arXiv paper,
  ## but will try to look more.  At this moment I will leave these commented-out.
  #for repeated simulations with p=30
  #beta <- runif(15, min = -7, max = 7)
  #beta=c(beta,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)


  sigmasquare <- 5
  # for repeated simulations
  #sigmasquare <- runif(1, min = 2, max = 10)
  tau <- sigmasquare^(-1) #Since Jags uses precision rather than variance
  
  # Y: response
  Y <- rep(0,N)
  for (i in 1:N){
    Y[i] <- rnorm(n=1, mean=X[i,]%*%beta,sd=sigmasquare^0.5)
  }


  # Data
  N <- length(Y)
  p <- dim(X)[2]
  
  # Initial values
  b0 <- rep(0,K)
  sigma2 <- 1
  lambda2=1
  sigma0 <- 1 
  nu0 <- 0.0001 

  return(list(constants = list(n=N,k=p),
              data = list(X=X,Y=Y),
              inits = list(beta=b0, sigmasq=sigma2, lambda2=lambda2),
              dimensions = list(a = c(N,p))))
}

#NIMBLE


linearCodeL <- nimbleCode({
  for(i in 1:n){
    for(j in 1:k){
      a[i,j] <- X[i,j]*beta[j]
    }
    mu[i] <- sum(a[i,])
    # likelihood
    Y[i] ~ dnorm(mu[i],sd=sigma)
  }
  ## priors
  for(i in 1:k){    
    b0[i] <- 0
  }
  
  nu0 <- 0.001
  sigma0 <- 1
  
  sigmasq ~ dinvgamma(nu0*0.5,nu0*sigma0*sigma0*0.5)
  tau <- pow(sigmasq,-1)  ## residual std dev
  sigma <- pow(sigmasq,0.5)
  for (j in 1:k){
    beta[j] ~ ddexp(location=0, scale=1.0 / (pow(lambda2,0.5)))
  }
  lambda2 ~ dexp(0.1)
})

linearCodeL_better <- nimbleCode({
  mu[1:n] <- (X[1:n,1:k] %*% beta[1:k])[1:n,1] # alternate
  for(i in 1:n){
    # likelihood
    Y[i] ~ dnorm(mu[i],sd=sigma)
  }
  ## priors
  for(i in 1:k){    
    b0[i] <- 0
  }
  
  nu0 <- 0.001
  sigma0 <- 1
  
  sigmasq ~ dinvgamma(nu0*0.5,nu0*sigma0*sigma0*0.5)
  tau <- pow(sigmasq,-1)  ## residual std dev
  sigma <- pow(sigmasq,0.5)
  for (j in 1:k){
    beta[j] ~ ddexp(location=0, scale=1.0 / (pow(lambda2,0.5)))
  }
  lambda2 ~ dexp(0.1)
})

run_nimble <- function(code, case) {
  t1<- system.time(
  {
     linear <- nimbleModel(code = code, 
                            dimensions = case$dimensions,
                            name = "linear",
                            constants = case$constants, 
                            data = case$data,
                            inits = case$inits)
     Clinear <- compileNimble(linear)
     linearConf <- configureMCMC(linear, print = TRUE)
     linearConf$addMonitors(c("tau", "sigmasq", "beta"))
     linearMCMC <- buildMCMC(linearConf)
     ClinearMCMC <- compileNimble(linearMCMC, project = linear)
     t2<- system.time(
     {
       Psamples <- runMCMC(ClinearMCMC, niter=20000, nburnin=10000, thin=2,
                           nchains=1,samplesAsCodaMCMC = TRUE,summary=TRUE)
     })
  })
  list(t1 = t1, t2 = t2, Psamples = Psamples, niter = 20000, nburnin = 10000)
}

cases <- list()
cases[["N=100, p=16, z=0"]]    <- sim_lm_data_L_case(N=100,  K=16,  number_zero_coefs = 0,  seed = 431234)
cases[["N=1000, p=16, z=0"]]   <- sim_lm_data_L_case(N=1000, K=16,  number_zero_coefs = 0,  seed = 431234)
cases[["N=1000, p=30, z=2"]]   <- sim_lm_data_L_case(N=1000, K=30,  number_zero_coefs = 2,  seed = 431234)
cases[["N=1000, p=30, z=15"]]  <- sim_lm_data_L_case(N=1000, K=30,  number_zero_coefs = 15, seed = 431234)
cases[["N=1000, p=30, z=28"]]  <- sim_lm_data_L_case(N=1000, K=30,  number_zero_coefs = 28, seed = 431234)
cases[["N=1000, p=100, z=2"]]  <- sim_lm_data_L_case(N=1000, K=100, number_zero_coefs = 2,  seed = 431234)
cases[["N=1000, p=100, z=50"]] <- sim_lm_data_L_case(N=1000, K=100, number_zero_coefs = 50, seed = 431234)
cases[["N=1000, p=100, z=98"]] <- sim_lm_data_L_case(N=1000, K=100, number_zero_coefs = 98, seed = 431234)

# It is recommented to generate results in multiple R sessions to reduce
# later-slower effects from many runs in memory.
# For example, in each session run the above code following by one
# of the blocks below:

# Block 1:
results_L_BFG <- list()
for(case_name in names(cases)) {
  results_L_BFG[[ case_name ]] <- run_nimble(linearCodeL, cases[[case_name]])
}
save(results_L_BFG, file = "results_L_BFG.RData")

# Block 2:
results_L_better <- list()
for(case_name in names(cases)) {
  results_L_better[[ case_name ]] <- run_nimble(linearCodeL_better, cases[[case_name]])
}
save(results_L_better, file = "results_L_better.RData")
