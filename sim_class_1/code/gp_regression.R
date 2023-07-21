library(FastGP)
library(matrixStats)

corr <- function(phi, D) exp(- phi * D)

sampleBetaSigma <- function(ytilde, xtilde) {
    ndim <- ncol(xtilde)
    nsample <- nrow(xtilde)

    xtxInv <- chol2inv(chol(crossprod(xtilde)))
    xty <- crossprod(xtilde, ytilde)

    betaHat <- drop(xtxInv %*% xty)
    yHat <- drop(xtilde %*% betaHat)

    rss <- sum((ytilde - yHat) * (ytilde - yHat))

    sig <- rss / rchisq(1, nsample - ndim)
    betas <- drop(rcpp_rmvnorm_stable(1, sig * xtxInv, betaHat))

    list(betas = betas, sig = sig)
}

transformPars <- function (pars, aphi = 0, bphi = 5) {

    thetas <- numeric(2)
    thetas[1] <- log((pars[1] - aphi) / (bphi - pars[1]))
    thetas[2] <- log(pars[2])

    thetas
}

logLikPars <- function (resids, sig, thetas, kernDistMat, aphi = 0, bphi = 2) {

    nsample <- length(resids)
    ## transformed parameters
    phi <- aphi + (bphi - aphi) / (1 + exp(-thetas[1]))
    alpha <- exp(thetas[2])

    tmp <- sig * (corr(phi, kernDistMat) + alpha * diag(1, nsample))
    tmp1 <- chol(tmp)
    tmp2 <- chol2inv(tmp1)

    lpost <- - 0.5 * sum(resids * (tmp2 %*% resids)) - sum(log(diag(tmp1)))

    lpost
}

sliceSample <- function(resids, sig, init, kernDistMat, aphi, bphi, priorMu, priorSd) {

    thetas <- transformPars(init, aphi, bphi)

    ## Initialize
    nu <- as.numeric(priorMu + priorSd %*% rnorm(length(init)))
    u <- runif(1)

    ### Update variables using slice sampling
    logLik <- logLikPars(resids, sig, thetas, kernDistMat, aphi, bphi)
    logy <- logLik + log(u)    ## slice threshold
    phase <- runif(1, 0, 2 * pi)                              ## draw an initial proposal for mixture element
    phase_min <- phase - 2*pi                                     ## lower bound for theta to start with
    phase_max <- phase                                          ## upper bound for theta to start with
    thetas_prime <- (thetas - priorMu) * cos(phase) + (nu - priorMu) * sin(phase) + priorMu  ## generate a proposal for transformed variables
    while(1) {
        logLik_prime <- logLikPars(resids, sig, thetas_prime, kernDistMat, aphi, bphi)

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
    res$alpha <- exp(thetas_prime[2])

    res
}

updateGP <- function (resids, sig, tau, phi, kernDistMatFull) {
    library(FastGP)

    ntrain <- length(resids)
    nsample <- nrow(kernDistMatFull)
    ntest <- nsample - ntrain

    gpCov <- sig * corr(phi, kernDistMatFull)

    trainTrainCov <- gpCov[1:ntrain, 1:ntrain]
    testTestCov <- gpCov[(ntrain + 1):nsample, (ntrain + 1):nsample]
    testTrainCov <- gpCov[(ntrain + 1):nsample, 1:ntrain]

    invTrainTrainCov <- chol2inv(chol(trainTrainCov + tau * diag(1, ntrain)))

    testPostMean <- (testTrainCov %*% invTrainTrainCov) %*% resids
    testPostCov <- testTestCov - tcrossprod(testTrainCov %*% invTrainTrainCov, testTrainCov)

    gpSamp <- drop(rcpp_rmvnorm_stable(1, testPostCov, testPostMean))

    gpSamp
}

fitGPR <- function (y, x, xtest, kernDistMatFull, priorMu, priorSd,
                    aphi = 0, bphi = 5,
                    niter = 3000, nburn = 1000, nthin = 2) {

    ntrain <- nrow(x)
    ndim <- ncol(x)
    cts <- 0

    sampBetas <- list()
    sampTau <- list()
    sampSig <- list()
    sampPhi <- list()
    sampAlpha <- list()
    sampGP <- list()
    sampY <- list()

    sBetas <- rep(NA, ndim)
    sTau <- NA
    sSig <- NA
    sPhi <- NA

    sTauPrev <- 1
    sSigPrev <- 1
    sPhiPrev <- 1
    sAlphaPrev <- sTauPrev / sSigPrev

    kernDistMat <- kernDistMatFull[1:ntrain, 1:ntrain]

    startTime <- proc.time()
    for (its in 1:niter) {
        rmat <- chol(exp(- sPhiPrev * kernDistMat) + sAlphaPrev * diag(1, ntrain))
        ytilde <- forwardsolve(rmat, x = y, upper.tri = TRUE, transpose = TRUE)
        xtilde <- forwardsolve(rmat, x = x, upper.tri = TRUE, transpose = TRUE)

        betasSig <- sampleBetaSigma(ytilde, xtilde)
        sBetas <- betasSig$betas
        sSig <- betasSig$sig

        resids <- as.numeric(y - x %*% sBetas)
        inits0 <- c(sPhiPrev, sAlphaPrev)
        sampPhiAlpha <- sliceSample(resids, sSig, inits0, kernDistMat, aphi, bphi, priorMu, priorSd)
        sPhi <- sampPhiAlpha$phi
        sAlpha <- sampPhiAlpha$alpha
        sTau <- sAlpha * sSig

        cat("iter: ", its,  "  sTau: ", round(sTau, 2), "  sSig: ", round(sSig, 2), "  sPhi: ", sPhi, " betas: ",
            round((sBetas), 2), "\n")
        ## if (its %% 10 == 0)

        if (its > nburn & its %% nthin == 0) {
            cts <- cts + 1
            sampBetas[[cts]] <- sBetas
            sampTau[[cts]] <- sTau
            sampSig[[cts]] <- sSig
            sampPhi[[cts]] <- sPhi
            sampAlpha[[cts]] <- sAlpha
            sampGP[[cts]] <- updateGP(resids, sSig, sTau, sPhi, kernDistMatFull)
            sampY[[cts]] <- as.numeric(xtest %*% sBetas) + sampGP[[cts]] + rnorm(nrow(xtest), mean = 0, sd = sqrt(sTau))
        }
        sTauPrev <- sTau
        sSigPrev <- sSig
        sPhiPrev <- sPhi
        sAlphaPrev <- sAlpha
    }
    endTime <- proc.time()

    list(
        beta = sampBetas,
        tau = unlist(sampTau),
        sig = unlist(sampSig),
        phi = unlist(sampPhi),
        alpha = unlist(sampAlpha),
        gp = sampGP,
        pred = sampY,
        time = endTime[3] - startTime[3]
    )
}
