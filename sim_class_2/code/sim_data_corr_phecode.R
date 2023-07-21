library(mvtnorm)
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

allsubstr <- function(x) {
    subs <- list()
    for (n in 1:nchar(x)) {
        subs[[n]] <- substring(x, 1:(nchar(x) - n + 1), n:nchar(x))
    }
    unlist(subs)
}

i <- pair[[ii]][1]
j <- pair[[ii]][2]
idx <- c(unlist(Idx[[i]]), unlist(Idx[[j]]))
ICD <- xx$icd[idx]
kernDistMat <- xx$kernDistMat[idx, idx]
SkMat <- xx$SkMat[idx, idx]

Label <- rep(1, length(ICD)) 
Label[c(unlist(Idx[[j]]))] <- -1

ndim <- 2
ntrain <- 1000
ntest <- 100
nsample <- ntrain + ntest
betas <- rep(c(0.1, 0.2), length = ndim)

set.seed(12345)

ICDtoPhe = read.csv("../Phecode/Phecode_map_v1_2_icd10cm_beta.csv")

simData <- list()
for (cc in 1:10) {
    xmat <- cbind(1, rnorm(nsample))
    colnames(xmat) <- c("Int", "Slope")
    
    linpred <- drop(xmat %*% betas)
    
    sampleIdx <- sample(1:length(ICD), nsample, replace = TRUE)
    
    icd <- ICD[sampleIdx]
    label <- Label[sampleIdx]
    
    icdsp = strsplit(icd, " ")
    ## transfer icd-10-cm to PheCode 
    phesp = list()
    for (i in 1:length(icdsp)) {
        phesp[[i]] = ICDtoPhe$phecode[match(unlist(icdsp[[i]]), ICDtoPhe$icd10cm)]
    }
    
    phe = icd
    for (i in 1:length(icd)) {
        phe[i] = paste(phesp[[ii]], collapse = " ")
    }
    
    
    ncol = length(unique(unlist(phesp)))
    ## Set up an empty matrix
    pheMat = matrix(0, nrow = nsample, ncol = ncol)  
    colnames(pheMat) = unique(unlist(phesp))
    ## Fill it in
    for (k in 1:nsample) {
        index <- which(colnames(pheMat) %in% unlist(phesp[[k]]))
        pheMat[k, index] <- 1    
    }
    
    icdspu <- unique(unlist(icdsp))
    icdsubl <- list()
    for (i in 1:length(icdspu)) {
        icdsubl[[i]] <- allsubstr(icdspu[i])
    }
    icdsub <- unique(unlist(icdsubl))
    ncol1 <- length(icdsub)
    ## Set up an empty matrix
    icdMatsub <- matrix(0, nrow = nsample, ncol = ncol1)  
    colnames(icdMatsub) <- icdsub
    ## Fill it in
    for (k in 1:nsample) {
        for (j in 1:length(icdsp[[k]])) {
            index <- which(colnames(icdMatsub) %in% allsubstr(unlist(icdsp[[k]])[j]))
            icdMatsub[k, index] <- icdMatsub[k, index] + 1
        }
    }
    icdMat <- icdMatsub
    
    skMat <- SkMat[sampleIdx, sampleIdx]
    cholKern <- rcppeigen_get_chol_stable(skMat)
    diagKern <- rcppeigen_get_chol_diag(skMat)
    diagKern[which(diagKern < 0)] <- 0
    
    ff <- drop(cholKern %*% (sqrt(diagKern) * rnorm(nsample)))
    
    x <- apply(icdMatsub, 1, function(x) sum(x^2))
    probs <- 1 / (1 + exp(- (3*tan(x) + 3*label + linpred)))
    
    y <- rbinom(n = length(probs), size = 1, prob = probs)
    
    trainId <- 1:ntrain
    testId <- (ntrain + 1):nsample
    
    ytrain <- y[trainId]
    xtrain <- xmat[trainId, ]
    xtrain.phe <- cbind.data.frame(xtrain, phe = phe[trainId])
    xtrain.feat <- cbind(xmat[trainId, ], pheMat[trainId, ])
    ptrain <- probs[trainId]
    
    ytest <- y[-trainId]
    xtest <- xmat[-trainId, ]
    xtest.phe <- cbind.data.frame(xtest, phe = phe[-trainId])
    xtest.feat <- cbind(xmat[-trainId, ], pheMat[-trainId, ])
    ptest <- probs[-trainId]
    
    simData[[cc]] <- list(
        phetrain = phe[trainId],
        ytrain = y[trainId],
        xtrain = xmat[trainId, ],
        xtrain.phe = cbind.data.frame(xtrain, phe = phe[trainId]),
        xtrain.feat = cbind(xmat[trainId, ], pheMat[trainId, ]),
        ptrain = probs[trainId],
        phetest = phe[-trainId],
        ytest = y[-trainId],
        xtest = xmat[-trainId, ],
        xtest.phe = cbind.data.frame(xtest, phe = phe[-trainId]),
        xtest.feat = cbind(xmat[-trainId, ], pheMat[-trainId, ]),
        ptest = probs[-trainId],
        KernMat = kernDistMat[sampleIdx, sampleIdx],
        ff = ff,
        x = x,
        label = label
    )
}

fname <- paste0("../data/sim_data_lin_corr_phecode", ii, ".rds")
saveRDS(simData, fname)

