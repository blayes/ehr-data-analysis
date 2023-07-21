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

set.seed(12345)

simData <- list()
for (cc in 1:10) {

    dat <- readRDS("../data/sim_data_nonlin_corr.rds")
    cvdat <- dat[[cc]]

    icdList <- c(cvdat$icdtrain, cvdat$icdtest)

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

    y <- icdPred + rnorm(nsample, 0, 0.1)

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

saveRDS(simData, "../data/sim_data_nonlin_corr_reg.rds")
