# This file contains code adapted from
# https://github.com/daniele-falco/software_comparison/tree/main/linear_models
#
# In this file we fix the way the model is written to allow nimble to 
# detect conjugacy.

library(nimble)
# library(extraDistr) # This package from Beraha et al's code conflicts with nimble's dinvgamma

# Put their data simulation code into a function
sim_lm_data_conj_case <- function(N = 100, K = 4, seed = 1234) {
  set.seed(seed)
  
  ####################SIMULATION OF DATA
  # N=100 #number of observations
  # K<-4  #length of our parameters
  
  # Matrix X of covariates with binary covariates
  X <- matrix(nrow=N, ncol=K)
  X[,1] <- rep(1,N) #intercept
  
  if(K==4){
    p_b=c(0.1,0.5,0.8)}
  if(K==16){
    p_b=c(0.1,0.2,0.3,0.4,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.6,0.7,0.8,0.9)}
  if(K==50){
    p_b=runif((K-1),min=0.05,max=0.95)}
  ## Beraha et al's code appears to be missing the K=30 case shown in the arXiv paper
  ## We will assume it is generated in the same way as K=50
  if(K==30){
    p_b=runif((K-1),min=0.05,max=0.95)}
  
  
  ##### Beraha et al's code has one block of code to generate
  ##### Bernoulli-distributed covariates, followed by another 
  ##### block of code that over-writes these with normally
  ##### distributed covariates.  
  ##### We assume the first one is for the "LM-C - Bin" case.
  ##### That is the case run in this file.
 for(i in 2:K){
   X[,i] <- extraDistr::rbern(N,p_b[i-1]) #covariates 
 }
  
  ## # Matrix X of covariates with continuous covariates
  ## X <- matrix(nrow=N, ncol=K)
  ## X[,1] <- rep(1,N) #intercept
  ## for(i in 2:K){
  ##   X[,i] <- rnorm(n = N, mean = 0 ,sd = 1) #covariates 
  ## }
  
  # TRUE values
  beta <- runif(K, min = -7, max = 7)
  sigmasquare <- runif(1, min = 2, max = 10)
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
  sigma0 <- 1 
  nu0 <- 0.0001 
  B0 <- diag(p)

  return(list(constants = list(n=N,k=p),
              data = list(X=X,Y=Y,B0=B0),
              inits = list(beta=b0, sigmasq=sigma0),
              dimensions = list(a = c(N,p))))
}

cases <- list()
cases[["N=100, p=4"]]   <- sim_lm_data_conj_case(N=100,  K=4,  seed = 1)
cases[["N=1000, p=4"]]  <- sim_lm_data_conj_case(N=1000, K=4,  seed = 1)
cases[["N=100, p=16"]]  <- sim_lm_data_conj_case(N=100,  K=16, seed = 1)
cases[["N=1000, p=16"]] <- sim_lm_data_conj_case(N=1000, K=16, seed = 1)
cases[["N=30, p=50"]]   <- sim_lm_data_conj_case(N=30,   K=50, seed = 1)

# Beraha et al's code
linearCodeZ <- nimbleCode({
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
  cova[1:k,1:k] <- sigmasq*B0[1:k,1:k]
  beta[1:k] ~ dmnorm(b0[1:k],cov=cova[1:k,1:k])   
  
})

# Simpler, better code
linearCodeZconj <- nimbleCode({
  mu[1:n] <- (X[1:n,1:k] %*% beta[1:k])[1:n,1] 
  for(i in 1:n) {
    # mu[i] <- inprod(X[i,], beta[1:k]) # alternate
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
  for(i in 1:k)
    beta[i] ~ dnorm(b0[i], var = sigmasq)
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
          Psamples <- runMCMC(ClinearMCMC, niter=11000, nburnin=1000, thin=2,
                              nchains=1,samplesAsCodaMCMC = TRUE,summary=TRUE)
        })
    })
  list(t1 = t1, t2 = t2, Psamples = Psamples, niter = 11000, nburnin = 1000)
}

# Run the above code to define functions and nimbleCode.
# Next, it is recommended to do each of the below blocks in separate R sessions.
# In each R session the above code needs to be run.
# Splitting these result runs into separate R sessions should reduce
# the chance of later-slower effects due to accumulated memory use or
# any other reasons.

# Block 1:
# Here is how we understand BFG to run this case.
# It is labeled "conj" but did not actually achieve conjugate sampling.
results_C_bin_BFG <- list()
for(case_name in names(cases)) {
  results_C_bin_BFG[[ case_name ]] <- run_nimble(linearCodeZ, cases[[case_name]])
}
save(results_C_bin_BFG, file = "results_C_bin_BFG.RData")

# Block 2:
# Here are results from making the model code simpler and better.
# This actually achieves conjugacy for the beta[i]s
results_C_bin_better <- list()
for(case_name in names(cases)) {
  results_C_bin_better[[ case_name ]] <- run_nimble(linearCodeZconj, cases[[case_name]])
}
save(results_C_bin_better, file = "results_C_bin_better.RData")
