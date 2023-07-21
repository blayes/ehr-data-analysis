setwd("~/icd_biom/code")
rm(list = ls())

## sample size, p, and the parametric linear predictor

ndim <- 2
ntrain <- 1000
ntest <- 100
nsample <- ntrain + ntest
betas <- rep(c(0.1, 0.2), length = ndim)

###
## There are four chonditions and each is represented using 2 letters
## divied them into five equal sized groups to define chronic condtions
popChrIcd <- list("chr1" = c("A", "B", "AA", "BB", "BA", "AB"),
                  "chr2" = c("C", "D", "CC", "DD", "CD", "DC"),
                  "chr3" = c("E", "F", "EE", "FF", "EF", "FE"),
                  "chr4" = c("G", "H", "GG", "HH", "GH", "HG")
                  )

## group 1 (y = 0) has chronic conditions 1, 3, 5 and the other group
## (y = 1) has chronic conditions 2 and 4.

set.seed(12345)

simData <- list()
for (cc in 1:10) {
    
    dat <- readRDS("../data/sim_data_lin_corr.rds")
    cvdat <- dat[[cc]]

    xtrain <- cvdat$xtrain
    xtest <- cvdat$xtest
    xmat <- rbind(xtrain, xtest)
    linpred <- drop(xmat %*% betas)
    
    icdList <- c(cvdat$icdtrain, cvdat$icdtest)

    icdMat <- matrix(0, nsample, length(unlist(popChrIcd)))
    colnames(icdMat) <- as.vector(unlist(popChrIcd))

    for (ii in 1:nsample) {
        idx <- which(colnames(icdMat) %in% unlist(icdList[[ii]]))
        icdMat[ii, idx] <- 1    
    }

    beta0 <- c(rep(-1, 6), rep(1, 6), rep(-1, 6), rep(1, 6))

    icdPred <- drop(icdMat %*% (beta0))

    y <- linpred + icdPred + rnorm(nsample, 0, 0.1)

    ## chr1 & chr3 => y = 0 and chr2 & chr4 => y = 1
    ## cbind(probs, as.numeric(sapply(icdList, function(x) names(x)[1] == "chr2")), y)

    trainId <- 1:ntrain
    testId <- (ntrain + 1):nsample

    ytrain <- y[trainId]
    cvdat$ptrain <- NULL
    
    ytest <- y[-trainId]
    cvdat$ptest <- NULL

    cvdat$ytrain <- ytrain
    cvdat$ytest <- ytest
    
    simData[[cc]] <- cvdat
}

saveRDS(simData, "../data/sim_data_lin_corr_reg.rds")

