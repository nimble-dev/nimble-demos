# This file contains code adapted from
# https://github.com/daniele-falco/software_comparison/tree/main/linear_models

################
## "WI" cases

library(nimble)

# Put their data simulation code into a function
sim_lm_data_WI_case <- function(N = 100, K = 4, seed = 1234) {
  set.seed(seed)

 # N=100 #number of observations
 # K<-4  #length of our parameters

  # Matrix X of covariates with continuous covariates
  X <- matrix(nrow=N, ncol=K)
  X[,1] <- rep(1,N) #intercept
  for(i in 2:K){
    X[,i] <- rnorm(n = N, mean = 0 ,sd = 1) #covariates 
  }
  
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
  sigma2 <- 1
  
  return(list(constants = list(n=N,k=p),
              data = list(X=X,Y=Y),
              inits = list(beta=b0, sigma=sigma2),
              dimensions = list(a = c(N,p))))
}

linearCodeWI <- nimbleCode({
  for(i in 1:n){
    for(j in 1:k){
      a[i,j] <- X[i,j]*beta[j]
    }
    mu[i] <- sum(a[i,])
    # likelihood
    Y[i] ~ dnorm(mu[i],sd=sigma)
  }
  ## priors
  for(j in 1:k){
    beta[j] ~ dnorm(0,sd=100)
  }  
  sigma ~ T(dt(mu = 0, tau =2.5, df=1), 0, )
})

linearCodeWI_better <- nimbleCode({
  mu[1:n] <- (X[1:n,1:k] %*% beta[1:k])[1:n,1] # alternate
  for(i in 1:n){
    # likelihood
    Y[i] ~ dnorm(mu[i],sd=sigma)
  }
  ## priors
  for(j in 1:k){
    beta[j] ~ dnorm(0,sd=100)
  } 
  sigma ~ T(dt(mu = 0, tau =2.5, df=1), 0, )
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
      linearConf$addMonitors(c("sigma", "beta"))
      linearMCMC <- buildMCMC(linearConf)
      ClinearMCMC <- compileNimble(linearMCMC, project = linear)
      t2<- system.time(
        {
          Psamples <- runMCMC(ClinearMCMC, niter=15000, nburnin=5000, thin=2,
                              nchains=1,samplesAsCodaMCMC = TRUE,summary=TRUE)
        })
    })
  list(t1 = t1, t2 = t2, Psamples = Psamples, niter = 15000, nburnin = 5000)
}

cases <- list()
cases[["N=100, p=4"]]   <- sim_lm_data_WI_case(N=100,  K=4,  seed = 431234) #431234 is in BFG's code
cases[["N=1000, p=4"]]  <- sim_lm_data_WI_case(N=1000, K=4,  seed = 431234)
cases[["N=100, p=16"]]  <- sim_lm_data_WI_case(N=100,  K=16, seed = 431234)
cases[["N=1000, p=16"]] <- sim_lm_data_WI_case(N=1000, K=16,  seed = 431234)

# It is recommented to generate results in multiple R sessions to reduce
# later-slower effects from many runs in memory.
# For example, in each session run the above code following by one
# of the blocks below:

# Block 1: 
results_WI_BFG <- list()
for(case_name in names(cases)) {
  results_WI_BFG[[ case_name ]] <- run_nimble(linearCodeWI, cases[[case_name]])
}
save(results_WI_BFG, file = "results_WI_BFG.RData")

# Block 2:
results_WI_better <- list()
for(case_name in names(cases)) {
  results_WI_better[[ case_name ]] <- run_nimble(linearCodeWI_better, cases[[case_name]])
}
save(results_WI_better, file = "results_WI_better.RData")

