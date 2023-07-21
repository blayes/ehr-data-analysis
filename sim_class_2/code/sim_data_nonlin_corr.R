setwd("~/icd_biom/sim_class_2/code")

rm(list = ls())

library(kernlab)

ntrain <- 1000
ntest <- 100
nsample <- ntrain + ntest

## There are four chonditions and each is represented using 2 letters
## divied them into five equal sized groups to define chronic condtions
popChrIcd <- list("chr1" = c("A", "B", "AA", "BB", "BA", "AB"),
                  "chr2" = c("C", "D", "CC", "DD", "CD", "DC"),
                  "chr3" = c("E", "F", "EE", "FF", "EF", "FE"),
                  "chr4" = c("G", "H", "GG", "HH", "GH", "HG")
                  )

## group 1 (y = 0) has chronic conditions 1, 3 and the other group
## (y = 1) has chronic conditions 2 and 4.

set.seed(12345)

simData <- list()
for (cc in 1:10) {

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
            chrCondNo <- drop(rmultinom(1, 6, c(0.5, 0.5)))
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
            chrCondNo <- drop(rmultinom(1, 6, c(0.5, 0.5)))
            icdList[[ii]]$chr2 <- sample(popChrIcd$chr2, chrCondNo[1])
            icdList[[ii]]$chr4 <- sample(popChrIcd$chr4, chrCondNo[2])        
        }        
    }
    icd <- sapply(icdList, function(x) paste(unlist(x), collapse = " "))
    
    icdMat <- matrix(0, nsample, length(unlist(popChrIcd)))
    colnames(icdMat) <- as.vector(unlist(popChrIcd))

    for (ii in 1:nsample) {
        idx <- which(colnames(icdMat) %in% unlist(icdList[[ii]]))
        icdMat[ii, idx] <- 1
    }
    

    icdFeat <- model.matrix(  ~ .^2, data = as.data.frame(icdMat))
    intCol <- c("BA:AB", "EF:FE", "CD:DC", "GH:HG")
    icdInt <- icdFeat[ , which(colnames(icdFeat) %in% intCol)]

    icdPred <- drop(icdInt %*% c(-1, 2, -1, 2))
    
    probs <- 1 / (1 + exp(- (icdPred)))

    y <- as.numeric(probs > 0.5)
    
    trainId <- 1:ntrain
    testId <- (ntrain + 1):nsample

    ytrain <- y[trainId]
    xtrain <- matrix(1, ntrain, 1)
    xtrain.icd <- cbind.data.frame(xtrain, icd = icd[trainId])
    xtrain.feat <- cbind(xtrain, icdMat[trainId, ])
    ptrain <- probs[trainId]
    
    ytest <- y[-trainId]
    xtest <- matrix(1, ntest, 1) 
    xtest.icd <- cbind.data.frame(xtest, icd = icd[-trainId])
    xtest.feat <- cbind(xtest, icdMat[-trainId, ])
    ptest <- probs[-trainId]
    
    simData[[cc]] <- list(
        icdtrain = icdList[trainId],
        ytrain = y[trainId],
        xtrain = xtrain,
        xtrain.icd = cbind.data.frame(xtrain, icd = icd[trainId]),
        xtrain.feat = cbind(xtrain, icdMat[trainId, ]),
        ptrain = probs[trainId],
        icdtest = icdList[-trainId],
        ytest = y[-trainId],
        xtest = xtest,
        xtest.icd = cbind.data.frame(xtest, icd = icd[-trainId]),
        xtest.feat = cbind(xtest, icdMat[-trainId, ]),
        ptest = probs[-trainId]
    )
}

saveRDS(simData, "../data/sim_data_nonlin_corr.rds")

