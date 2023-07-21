library(FastGP)
library(glmnet)
library(kernlab)
library(ranger)

cmdArgs <- commandArgs(trailingOnly = TRUE)
ii <- as.numeric(cmdArgs)

xx <- readRDS("../data/data.rds")
lvls <- unique(xx$cancer)
pair <- list()

for (i in 1:5) {
    for (j in (i+1):6) {
        name <- paste(lvls[i], "&", lvls[j], sep = " ")
        pair[[name]] <- c(i, j)
    }
}

Idx <- list()
for (ll in 1:length(lvls)) {
    Idx[[ll]] <- which(lvls[ll] == xx$cancer)
}

i <- pair[[ii]][1]
j <- pair[[ii]][2]
idx <- c(unlist(Idx[[i]]), unlist(Idx[[j]]))
ICD <- xx$icd[idx]
kernDistMat <- xx$kernDistMat[idx, idx]
SkMat <- xx$SkMat[idx, idx]


ndim <- 2
ntrain <- 1000
ntest <- 100
nsample <- ntrain + ntest
betas <- rep(c(0.1, 0.2), length = ndim)

set.seed(12345)

simData <- list()
for (cc in 1:10) {
    
    dat <- readRDS(paste0("../data/sim_data_lin_corr_phecode", ii, ".rds"))
    cvdat <- dat[[cc]]
    
    xtrain <- cvdat$xtrain
    xtest <- cvdat$xtest
    xmat <- rbind(xtrain, xtest)
    linpred <- drop(xmat %*% betas)
    ff <- cvdat$ff
    
                                        # y <- linpred + ff + rnorm(nsample, 0, 0.1)
    y <- 3*tan(cvdat$x) + 3*(cvdat$label) + linpred
    
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
fname <- paste0("../data/sim_data_lin_corr_reg_", ii, ".rds")
saveRDS(simData, fname)
