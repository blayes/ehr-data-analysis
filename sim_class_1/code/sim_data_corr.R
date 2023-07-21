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
    xmat <- cbind(1, rnorm(nsample))
    colnames(xmat) <- c("Int", "Slope")

    linpred <- drop(xmat %*% betas)
    
    icdList <- vector("list", nsample)
    for (ii in 1:nsample) {
        if (runif(1) > 0.5) {        
            icdList[[ii]] <- list("chr1" = "",
                                  "chr2" = "",
                                  "chr3" = "",
                                  "chr4" = ""                              
                                  )
            chrCondNo <- drop(rmultinom(1, 6, c(0.5, 0.5)))
            icdList[[ii]]$chr2 <- sample(popChrIcd$chr2, chrCondNo[1])
            icdList[[ii]]$chr4 <- sample(popChrIcd$chr4, chrCondNo[2])
            chrCondNo <- drop(rmultinom(1, 3, c(0.5, 0.5)))
            icdList[[ii]]$chr1 <- sample(popChrIcd$chr1, chrCondNo[1])
            icdList[[ii]]$chr3 <- sample(popChrIcd$chr3, chrCondNo[2])        
        } else {
            icdList[[ii]] <- list("chr1" = "",
                                  "chr2" = "",
                                  "chr3" = "",
                                  "chr4" = ""                              
                                  )        
            chrCondNo <- drop(rmultinom(1, 6, c(0.5, 0.5)))
            icdList[[ii]]$chr1 <- sample(popChrIcd$chr1, chrCondNo[1])
            icdList[[ii]]$chr3 <- sample(popChrIcd$chr3, chrCondNo[2])
            chrCondNo <- drop(rmultinom(1, 3, c(0.5, 0.5)))
            icdList[[ii]]$chr2 <- sample(popChrIcd$chr2, chrCondNo[1])
            icdList[[ii]]$chr4 <- sample(popChrIcd$chr4, chrCondNo[2])        
        }        
    }

    icdMat <- matrix(0, nsample, length(unlist(popChrIcd)))
    colnames(icdMat) <- as.vector(unlist(popChrIcd))

    for (ii in 1:nsample) {
        idx <- which(colnames(icdMat) %in% unlist(icdList[[ii]]))
        icdMat[ii, idx] <- 1    
    }

    beta0 <- c(rep(-1, 6), rep(1, 6), rep(-1, 6), rep(1, 6))

    icdPred <- drop(icdMat %*% (beta0))

    probs <- 1 / (1 + exp(- (linpred + icdPred)))

    y <- rbinom(n = length(probs), size = 1, prob = probs)

    ## chr1 & chr3 => y = 0 and chr2 & chr4 => y = 1
    ## cbind(probs, as.numeric(sapply(icdList, function(x) names(x)[1] == "chr2")), y)

    trainId <- 1:ntrain
    testId <- (ntrain + 1):nsample

    icd <- sapply(icdList, function(x) paste(unlist(x), collapse = " "))

    ytrain <- y[trainId]
    xtrain <- xmat[trainId, ]
    xtrain.icd <- cbind.data.frame(xtrain, icd = icd[trainId])
    xtrain.feat <- cbind(xmat[trainId, ], icdMat[trainId, ])
    ptrain <- probs[trainId]

    ytest <- y[-trainId]
    xtest <- xmat[-trainId, ]
    xtest.icd <- cbind.data.frame(xtest, icd = icd[-trainId])
    xtest.feat <- cbind(xmat[-trainId, ], icdMat[-trainId, ])
    ptest <- probs[-trainId]

    simData[[cc]] <- list(
        icdtrain = icdList[trainId],
        ytrain = y[trainId],
        xtrain = xmat[trainId, ],
        xtrain.icd = cbind.data.frame(xtrain, icd = icd[trainId]),
        xtrain.feat = cbind(xmat[trainId, ], icdMat[trainId, ]),
        ptrain = probs[trainId],
        icdtest = icdList[-trainId],
        ytest = y[-trainId],
        xtest = xmat[-trainId, ],
        xtest.icd = cbind.data.frame(xtest, icd = icd[-trainId]),
        xtest.feat = cbind(xmat[-trainId, ], icdMat[-trainId, ]),
        ptest = probs[-trainId]
    )
}

saveRDS(simData, "../data/sim_data_lin_corr.rds")

