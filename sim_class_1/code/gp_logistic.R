library(pgdraw)
library(MASS)

corr <- function(phi, D) exp(- phi * D)

updateOmegas <- function(xtrain, betas, ftrain) {
    pgdraw(1, abs(drop(xtrain %*% betas + ftrain)))    
}

updateTestGP <- function (resids, sig, phi, omegas, kernDistMatFull) {
    library(FastGP)
    
    ntrain <- length(resids)
    nsample <- nrow(kernDistMatFull)
    ntest <- nsample - ntrain

    gpCov <- sig * corr(phi, kernDistMatFull)
    trainTrainCov <- gpCov[1:ntrain, 1:ntrain]
    testTestCov <- gpCov[(ntrain + 1):nsample, (ntrain + 1):nsample]
    testTrainCov <- gpCov[(ntrain + 1):nsample, 1:ntrain]

    invTrainTrainCov <- chol2inv(chol(trainTrainCov + diag(1 / omegas)))

    testPostMean <- (testTrainCov %*% invTrainTrainCov) %*% resids
    testPostCov <- testTestCov - tcrossprod(testTrainCov %*% invTrainTrainCov, testTrainCov)
    
    gpSamp <- drop(rcpp_rmvnorm_stable(1, testPostCov, testPostMean))

    gpSamp
}

updateTrainGP <- function (resids, sig, phi, omegas, kernDistMatTrain) {
    library(FastGP)
    
    ntrain <- nrow(kernDistMatTrain)

    gpCov <- sig * corr(phi, kernDistMatTrain)
    
    trainPostCov <- gpCov - gpCov %*% chol2inv(chol(gpCov + diag(1 / omegas))) %*% gpCov
    trainPostMean <- trainPostCov %*% (resids * omegas)
        
    gpSamp <- drop(rcpp_rmvnorm_stable(1, trainPostCov, trainPostMean))

    gpSamp
}

sampleBetas <- function(ztilde, xtilde) {
  ndim <- ncol(xtilde)
  
  xtx <- crossprod(xtilde)
  cholXtx <- tryCatch(chol(xtx), error = function(e) "error") ## take care of rank deficiency
  while (identical(cholXtx, "error")) {
    xtx <- xtx + diag(1e-8, ndim)
    cholXtx <- tryCatch(chol(xtx), error = function(e) "error")
  }
  
  xtxInv <- chol2inv(cholXtx)
  xtz <- crossprod(xtilde, ztilde)
  
  betaHat <- drop(xtxInv %*% xtz)
  
  drop(rcpp_rmvnorm_stable(1, xtxInv, betaHat))
}

transformPars <- function (pars, aphi = 0, bphi = 5) {
    
    thetas <- numeric(2)
    thetas[1] <- log((pars[1] - aphi) / (bphi - pars[1])) ## phi
    thetas[2] <- log(pars[2]) ## sig
    
    thetas
}

logLikPars <- function (resids, thetas, omegas, kernDistMatTrain, aphi = 0, bphi = 2) {

    ## transformed parameters
    phi <- aphi + (bphi - aphi) / (1 + exp(-thetas[1]))
    sig <- exp(thetas[2])

    tmp <- sig * corr(phi, kernDistMatTrain) + diag(1 / omegas)
    tmp1 <- chol(tmp)
    tmp2 <- chol2inv(tmp1)

    lpost <- - 0.5 * sum(resids * (tmp2 %*% resids)) - sum(log(diag(tmp1)))
    
    lpost
}

sliceSample <- function(resids, init, omegas, kernDistMatTrain, aphi, bphi, priorMu, priorSd) {

    thetas <- transformPars(init, aphi, bphi)

    ## Initialize
    nu <- as.numeric(priorMu + priorSd %*% rnorm(length(init)))
    u <- runif(1)

    ### Update variables using slice sampling
    logLik <- logLikPars(resids, thetas, omegas, kernDistMatTrain, aphi, bphi)
    logy <- logLik + log(u)    ## slice threshold
    phase <- runif(1, 0, 2 * pi)                              ## draw an initial proposal for mixture element
    phase_min <- phase - 2*pi                                     ## lower bound for theta to start with
    phase_max <- phase                                          ## upper bound for theta to start with
    thetas_prime <- (thetas - priorMu) * cos(phase) + (nu - priorMu) * sin(phase) + priorMu  ## generate a proposal for transformed variables
    while(1) {
        logLik_prime <- logLikPars(resids, thetas_prime, omegas, kernDistMatTrain, aphi, bphi)

        if(logLik_prime > logy){
            break
        }
        ## if loop does not break, shrink the interval
        if(phase<0){
            phase_min <- phase
        }else{
            phase_max <- phase
        }
        phase <- runif(1, phase_min, phase_max)
        thetas_prime <- (thetas - priorMu) * cos(phase) + (nu-priorMu) * sin(phase) + priorMu
    }

    res <- list()
    res$phi <- aphi + (bphi - aphi) / (1 + exp(-thetas_prime[1]))
    res$sig <- exp(thetas_prime[2])

    res
}

fitGPClass <- function (y, x, xtest, kernDistMatFull,
                        priorMu, priorSd,
                        aphi = 0, bphi = 5,
                        niter = 3000, nburn = 1000, nthin = 2) {
    library(FastGP)
    
    ntrain <- nrow(x)
    ndim <- ncol(x)
    cts <- 0

    sampBetas <- list()
    sampSig <- list()
    sampPhi <- list()    
    sampGP <- list()
    sampProb <- list()

    sFtrain <- rep(0, ntrain)
    sBetas <- rep(0, ndim)
    sSig <- NA
    sPhi <- NA
    sOmegas <- rep(NA, ntrain)
    
    sSigPrev <- 1
    sPhiPrev <- 1

    kernDistMatTrain <- kernDistMatFull[1:ntrain, 1:ntrain]
    
    startTime <- proc.time()
    for (its in 1:niter) {
        sOmegas <- updateOmegas(x, sBetas, sFtrain)
        
        zz <- (y - 0.5) / sOmegas        
        rmat <- chol(sSigPrev * exp(- sPhiPrev * kernDistMatTrain) + diag(1 / sOmegas))
        ztilde <- forwardsolve(rmat, x = zz, upper.tri = TRUE, transpose = TRUE)
        xtilde <- forwardsolve(rmat, x = x, upper.tri = TRUE, transpose = TRUE)

        sBetas <- sampleBetas(ztilde, xtilde) 
        
        resids <- zz - drop(x %*% sBetas)
        inits0 <- c(sPhiPrev, sSigPrev)
        sampPhiSig <- sliceSample(resids, inits0, sOmegas, kernDistMatTrain, aphi, bphi, priorMu, priorSd)
        sPhi <- sampPhiSig$phi
        sSig <- sampPhiSig$sig

        sFtrain <- updateTrainGP(resids, sSig, sPhi, sOmegas, kernDistMatTrain)     

        
        if (its %% 10 == 0) cat("iter: ", its,  "  sSig: ", round(sSig, 2), "  sPhi: ", sPhi, "\n")

        if (its > nburn & its %% nthin == 0) {
            cts <- cts + 1
            sampBetas[[cts]] <- sBetas
            sampSig[[cts]] <- sSig
            sampPhi[[cts]] <- sPhi
            sampGP[[cts]] <- updateTestGP(resids, sSig, sPhi, sOmegas, kernDistMatFull)     
            sampProb[[cts]] <- 1 / (1 + exp(- (as.numeric(xtest %*% sBetas) + sampGP[[cts]])))
        }
        sSigPrev <- sSig
        sPhiPrev <- sPhi
    }
    endTime <- proc.time()

    list(
        beta = do.call(rbind, sampBetas),
        sig = unlist(sampSig),
        phi = unlist(sampPhi),
        gp = do.call(rbind, sampGP),
        prob = do.call(rbind, sampProb),
        time = endTime[3] - startTime[3]
    )
}
