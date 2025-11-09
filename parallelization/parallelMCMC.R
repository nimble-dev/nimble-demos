## This provides a wrapper function to run NIMBLE's MCMC in parallel over multiple chains.
## The interface is similar to that of `runMCMC` and `nimbleMCMC`.
## `mcmcConfig` is an optional argument that allows one to specify an MCMC configuration. It
## should take in a model and return an MCMC configuration object (e.g., by calling `configureMCMC`).
## `sourceFiles` is an optional argument that supports `mcmcConfig`.
## It should be a character vector of filenames to `source()` that contain objects
## used in the MCMC configuration.
## `ncores` by default uses the number of cores detected on the machine but can be set to a different value.

nimbleParallelMCMC <- function (code, constants = list(), data = list(), inits = list(), dimensions = list(),
    mcmcConfig = NULL, monitors = NULL, monitors2 = NULL, thin = 1, thin2 = 1, niter = 10000, nburnin = 0, nchains = 1, 
    check = FALSE, calculate = FALSE, setSeed = FALSE,
    samples = TRUE, samplesAsCodaMCMC = FALSE, summary = FALSE, 
    WAIC = FALSE, perChainWAIC = FALSE, sourceFiles = NULL, ncores = parallel::detectCores()) {
    ## `mcmcConfig` should be a function that takes a model (and `...`) and outputs an MCMC configuration.
    ## `sourceFiles` should be a character vector of files to `source()` that contain objects
    ## used in the MCMC configuration.
    
    if (missing(code) || inherits(code, "modelBaseClass"))
        stop("nimbleParallelMCMC: must provide model code; cannot be run with actual model, unlike `nimbleMCMC`")
    if (!samples && !summary && !WAIC) 
        stop("no output specified, use samples = TRUE, summary = TRUE, or WAIC = TRUE")

    if(nchains < 2)
        stop("nimbleParallelMCMC requires at least two chains")

    
    if(length(inits)) {
        if(!is.function(inits) && !is.list(inits))
            stop('inits must be a function, a list of initial values, or a list (of length nchains) of lists of inital values')
        if(is.list(inits) && (length(inits) > 0) && is.list(inits[[1]]) && (length(inits) != nchains))
            stop('inits must be a function, a list of initial values, or a list (of length nchains) of lists of inital values')
    }

    hasMonitors2 <- !is.null(monitors2)
    
    if(ncores < nchains)
        cat("Fewer cores available than number of MCMC chains; at most ", ncores, " chains will run in parallel at a given time.\n")

    if(ncores > nchains)
        ncores <- nchains

    if(length(setSeed) > 1)
        stop("`setSeed` should be a single numeric value or `FALSE`")
    if(!setSeed) iseed <- NULL else iseed <- setSeed
    thisCluster <- parallel::makeCluster(ncores)
    ## Set independent random number streams.
    parallel::clusterSetRNGStream(cl = thisCluster, iseed = iseed)

    runOneMCMC <- function(chainNum, code, constants, data, inits, dimensions, mcmcConfig,
                           monitors, monitors2, thin, thin2, niter, nburnin, check, calculate, WAIC, sourceFiles) {
        library(nimble)

        tmp <- sapply(sourceFiles, source)
        
        if(is.function(inits))
            inits <- inits()
        if(is.list(inits) && length(inits) > 0 && is.list(inits[[1]])) {
            inits <- inits[[chainNum]]
        } 
        
        myModel <- nimbleModel(code = code,
                               constants = constants,
                               data = data,
                               dimensions = dimensions,
                               inits = inits,
                               check = check,
                               calculate = calculate)

        
        if(WAIC) {
            neededMonitors <- myModel$getParents(myModel$getNodeNames(dataOnly = TRUE), stochOnly = TRUE)
            if(!all(neededMonitors %in% monitors)) {
                monitors2 <- unique(c(monitors2, neededMonitors))
            }
        }

        if(is.null(mcmcConfig)) {
            conf <- configureMCMC(myModel, monitors2 = monitors2)
        } else conf <- mcmcConfig(myModel, monitors2 = monitors2)
        
        ## Plug in stuff given we want control of it here rather than via `mcmcConfig`.
        conf$enableWAIC <- WAIC
        conf$thin <- thin
        conf$thin2 <- thin2
        
        myMCMC <- buildMCMC(conf)

        CmyModel <- compileNimble(myModel)
        CmyMCMC <- compileNimble(myMCMC)
        CmyMCMC$run(niter, nburnin = nburnin, thin = thin, thin2 = thin2)

        return(list(smp1 = as.matrix(CmyMCMC$mvSamples), smp2 = as.matrix(CmyMCMC$mvSamples2)))
    }

    samplesList  <- vector('list', nchains); names(samplesList)  <- paste0('chain', 1:nchains)
    samplesList2 <- vector('list', nchains); names(samplesList2) <- paste0('chain', 1:nchains)
    
    output <- parLapply(cl = thisCluster, seq_len(nchains),
                        fun = runOneMCMC, code = code, constants = constants, data = data, inits = inits,
                        dimensions = dimensions, mcmcConfig = mcmcConfig, monitors = monitors, monitors2 = monitors2,
                        thin = thin, thin2 = thin2, niter = niter, nburnin = nburnin,
                        check = check, calculate = calculate, WAIC = WAIC, sourceFiles = sourceFiles)


    for(i in seq_len(nchains)) {
        samplesList[[i]] <- output[[i]][[1]]
        samplesList2[[i]] <- output[[i]][[2]]
    }
                                        
    if(WAIC) {
        ## We need to compiled model to be able to run offline WAIC.
        myModel <- nimbleModel(code = code,
                               constants = constants,
                               data = data,
                               dimensions = dimensions,
                               check = FALSE,
                               calculate = FALSE)
        
        CmyModel <- compileNimble(myModel)

        neededMonitors <- myModel$getParents(myModel$getNodeNames(dataOnly = TRUE), stochOnly = TRUE)
        if(all(neededMonitors %in% myModel$getVarNames(nodes = colnames(samplesList[[1]]))))
           samplesListToUse <- samplesList else samplesListToUse <- samplesList2
 
        if(perChainWAIC) {
            perChainWAICvalue <- list(); length(perChainWAICvalue) <- nchains
            for(i in seq_along(samplesListToUse)) {
                perChainWAICvalue[[i]] <- calculateWAIC(samplesListToUse[[i]] , myModel)
            }
        }

        samplesPerChain <- dim(samplesListToUse[[1]])[1]
        posteriorSamplesMatrix <- matrix(0, nrow = samplesPerChain*nchains, ncol = dim(samplesListToUse[[1]])[2])
        for(i in seq_along(samplesListToUse)) {
            posteriorSamplesMatrix[((i-1)*samplesPerChain + 1):(i*samplesPerChain),] <- samplesListToUse[[i]][,]
        }
        colnames(posteriorSamplesMatrix) <- colnames(samplesListToUse[[1]])
        WAICvalue <- calculateWAIC(posteriorSamplesMatrix, myModel)
    }
    if(samplesAsCodaMCMC) {
        if(!require(coda))
            stop("`coda` must be installed when `samplesAsCodaMCMC = TRUE`")
        samplesList <- as.mcmc.list(lapply(samplesList, as.mcmc))
        if(hasMonitors2)
            samplesList2 <- as.mcmc.list(lapply(samplesList2, as.mcmc))
    }
    if(summary) {
        summaryObject <- lapply(samplesList, samplesSummary)
        names(summaryObject) <- paste0('chain', 1:nchains)
        summaryObject$all.chains <- samplesSummary(do.call('rbind', samplesList))
        if(hasMonitors2) {
            summaryObject2 <- lapply(samplesList2, samplesSummary)
            summaryObject2$all.chains <- samplesSummary(do.call('rbind', samplesList2))
            summaryObject <- mapply(rbind, summaryObject, summaryObject2, SIMPLIFY = FALSE)      ## combine summaries
        }
    }

    parallel::stopCluster(thisCluster)
    
    retList <- list()
    if(samples) {
        retList$samples <- samplesList
        if(hasMonitors2)   retList$samples2 <- samplesList2
    }
    if(summary)   retList$summary <- summaryObject
    if(WAIC)      retList$WAIC    <- WAICvalue
    if(perChainWAIC) retList$perChainWAIC <- perChainWAICvalue
    if(length(retList) == 1) retList <- retList[[1]]
    return(retList)

}
