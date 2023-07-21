cmdArgs <- commandArgs(trailingOnly = TRUE)

mtd <- as.numeric(cmdArgs[1])
id <- as.numeric(cmdArgs[2])


if (mtd == 1) {
    library(matrixStats)
    library(kernlab)
    library(glmnet)

    dat <- readRDS("../data/sim_data_lin_corr.rds")
    cvdat <- dat[[id]]

    ytrain <- cvdat$ytrain
    xtrain <- cvdat$xtrain
    xtrain.icd <- cvdat$xtrain.icd
    xtrain.feat <- cvdat$xtrain.feat
    ptrain <- cvdat$ptrain

    ytest <- cvdat$ytest
    xtest <- cvdat$xtest
    xtest.icd <- cvdat$xtest.icd
    xtest.feat <- cvdat$xtest.feat
    ptest <- cvdat$ptest

    ## logistic regression
    dftrain = cbind.data.frame(y = ytrain, xtrain.feat)
    dftest = cbind.data.frame(y = ytest, xtest.feat)
    fit.lg <- glm(y ~ 0 + ., family = "binomial", data = dftrain)
    pred.lg <- predict(fit.lg, dftest, type = "response")

    ## ridge regression
    fit.r <- cv.glmnet(xtrain.feat, ytrain, alpha = 0, family = "binomial")
    pred.r <- predict(fit.r, newx = xtest.feat, s = "lambda.min", type = "response")

    ## lasso regression
    fit.lasso <- cv.glmnet(xtrain.feat, ytrain, alpha = 1, family = "binomial")
    pred.lasso <- predict(fit.lasso, newx = xtest.feat, s = "lambda.min", type = "response")

    ## spectrum and boundrange kernel
    sks <- stringdot(type= "spectrum")
    skb <- stringdot(type= "boundrange")

    fit.ksvp.sp <- ksvm(as.character(xtrain.icd[ , 3]), as.factor(ytrain), kernel = sks, scale = c(), prob.model = TRUE)
    fit.ksvp.bd <- ksvm(as.character(xtrain.icd[ , 3]), as.factor(ytrain), kernel = skb, scale = c(), prob.model = TRUE)

    fit.gp.sp <- gausspr(as.character(xtrain.icd[ , 3]), as.factor(ytrain), kernel = sks, scale = c(), prob.model = TRUE)
    fit.gp.bd <- gausspr(as.character(xtrain.icd[ , 3]), as.factor(ytrain), kernel = skb, scale = c(), prob.model = TRUE)

    pred.ksvp.sp <- predict(fit.ksvp.sp, as.character(xtest.icd[ , 3]), type = "probabilities")[ , 2]
    pred.ksvp.bd <- predict(fit.ksvp.bd, as.character(xtest.icd[ , 3]), type = "probabilities")[ , 2]

    pred.gp.sp <- predict(fit.gp.sp, as.character(xtest.icd[ , 3]), type = "probabilities")[ , 2]
    pred.gp.bd <- predict(fit.gp.bd, as.character(xtest.icd[ , 3]), type = "probabilities")[ , 2]

    res <- list(
        logis = list(fit = fit.lg, pred = pred.lg),
        ridge = list(fit = fit.r, pred = pred.r),
        lasso = list(fit = fit.lasso, pred = pred.lasso),
        svm.sp = list(fit = fit.ksvp.sp, pred = pred.ksvp.sp),
        svm.bd = list(fit = fit.ksvp.bd, pred = pred.ksvp.bd),
        gp.sp = list(fit = fit.gp.sp, pred = pred.gp.sp),
        gp.bd = list(fit = fit.gp.bd, pred = pred.gp.bd)
    )

    fname <- paste0("../data/comp_sim1_", id, ".rds")
    saveRDS(res, fname)
} else if (mtd == 2) {

    estimateCodeFreq <- function(dx) {
        unq_codes = as.character(dx)
        freq_codes = aggregate(dx, list(unq_codes), length)
        colnames(freq_codes) = c("codes", "freq")
        cl_list =  vector("list", nrow(freq_codes))
        names(cl_list) = freq_codes$codes
        for (i in 1:nrow(freq_codes)) {
            len = nchar(freq_codes$codes[i])
            cl <- data.frame(matrix(ncol = 2, nrow = len))
            cl$X2 = freq_codes$freq[i]
            for(j in 1:len){
                cl$X1[j] = substr(freq_codes$codes[i], 1, j)
            }
            cl_list[[i]] = cl
        }

        sc_freq = do.call(rbind, cl_list)
        sc_freq = aggregate(sc_freq$X2, by = list(sc_freq$X1), sum)
        colnames(sc_freq) = c("subtrees", "freq")
        sc_freq
    }

    popChrIcd <- list("chr1" = c("A", "B", "AA", "BB", "BA", "AB"),
                      "chr2" = c("C", "D", "CC", "DD", "CD", "DC"),
                      "chr3" = c("E", "F", "EE", "FF", "EF", "FE"),
                      "chr4" = c("G", "H", "GG", "HH", "GH", "HG")
                      )

    code_freq = estimateCodeFreq(unlist(popChrIcd))

    ##
    kernFunction = function (code1, code2, code_freq) {
        if ((length(code1) > 0) & (length(code2) > 0)) {
            f1 = estimateCodeFreq(code1)
            f2 = estimateCodeFreq(code2)
            mfs = merge(f1, f2, by = "subtrees")
            if(nrow(mfs) == 0) {
                return(0)
            } else{
                mfsOverall = merge(mfs, code_freq, by = "subtrees")
                return(sum(mfsOverall$freq.x * mfsOverall$freq.y
                           * sapply(mfsOverall$subtrees, nchar)
                           ))
            }
        } else {
            return(0)
        }
    }

    kernelMatrix = function (p_code1, p_code2, code_freq) {
        pl1 = length(p_code1)
        pl2 = length(p_code2)

        K = matrix(NA, ncol = pl1, nrow = pl2)
        for (i in 1:pl1) {
            cat("i = ", i, "\n")
            for (j in i:pl2) {
                K [i, j] = kernFunction(p_code1[[i]]$chr1, p_code2[[j]]$chr1, code_freq) + kernFunction(p_code1[[i]]$chr2, p_code2[[j]]$chr2, code_freq) +
                    kernFunction(p_code1[[i]]$chr3, p_code2[[j]]$chr3, code_freq) + kernFunction(p_code1[[i]]$chr4, p_code2[[j]]$chr4, code_freq)
            }
        }
        K[lower.tri(K)] <- t(K)[lower.tri(K)]
        return(K)
    }

    dat <- readRDS("../data/sim_data_lin_corr.rds")
    cvdat <- dat[[id]]

    ytrain <- cvdat$ytrain
    xtrain <- cvdat$xtrain

    xtest <- cvdat$xtest

    icdList <- c(cvdat$icdtrain, cvdat$icdtest)

    simMat <- kernelMatrix (icdList, icdList, code_freq)
    psim <- simMat / sqrt(outer(diag(simMat), diag(simMat), "*"))
    dsim <- outer(diag(psim), diag(psim), "+") - 2 * psim

    kernDistMat <- dsim ## this is the distance matrix computed from the similarity matrix

    saveRDS(kernDistMat, paste0("../data/kernDistMat_class_", id, ".rds"))

    source("gp_logistic.R")

    res <-
        fitGPClass (ytrain, xtrain, xtest, kernDistMat,
                    c(0, 0), diag(1, 2),
                    aphi = 0, bphi = 4,
                    niter = 10000, nburn = 5000, nthin = 5)

    fname <- paste0("../data/gp_sim1_", id, ".rds")

    saveRDS(list(fit = res, sim = simMat), fname)
} else if (mtd == 3) {
    library(matrixStats)
    library(kernlab)
    library(glmnet)

    dat <- readRDS("../data/sim_data_lin_corr_reg.rds")
    cvdat <- dat[[id]]

    ytrain <- cvdat$ytrain
    xtrain <- cvdat$xtrain
    xtrain.icd <- cvdat$xtrain.icd
    xtrain.feat <- cvdat$xtrain.feat

    ytest <- cvdat$ytest
    xtest <- cvdat$xtest
    xtest.icd <- cvdat$xtest.icd
    xtest.feat <- cvdat$xtest.feat

    ## lin regression
    dftrain = cbind.data.frame(y = ytrain, xtrain.feat)
    dftest = cbind.data.frame(y = ytest, xtest.feat)
    fit.lm <- lm(y ~ 0 + ., data = dftrain)
    pred.lm <- predict(fit.lm, dftest)

    ## ridge regression
    fit.r <- cv.glmnet(xtrain.feat, ytrain, alpha = 0)
    pred.r <- predict(fit.r, newx = xtest.feat, s = "lambda.min", type = "response")

    ## lasso regression
    fit.lasso <- cv.glmnet(xtrain.feat, ytrain, alpha = 1)
    pred.lasso <- predict(fit.lasso, newx = xtest.feat, s = "lambda.min", type = "response")

    ## spectrum and boundrange kernel
    sks <- stringdot(type= "spectrum")
    skb <- stringdot(type= "boundrange")

    fit.ksvp.sp <- ksvm(as.character(xtrain.icd[ , 3]), ytrain, kernel = sks, scale = c())
    fit.ksvp.bd <- ksvm(as.character(xtrain.icd[ , 3]), ytrain, kernel = skb, scale = c())

    fit.gp.sp <- gausspr(as.character(xtrain.icd[ , 3]), ytrain, kernel = sks, scale = c(), variance.model = TRUE)
    fit.gp.bd <- gausspr(as.character(xtrain.icd[ , 3]), ytrain, kernel = skb, scale = c(), variance.model = TRUE)

    pred.ksvp.sp <- predict(fit.ksvp.sp, as.character(xtest.icd[ , 3]), type = "response")
    pred.ksvp.bd <- predict(fit.ksvp.bd, as.character(xtest.icd[ , 3]), type = "response")

    pred.gp.sp <- predict(fit.gp.sp, as.character(xtest.icd[ , 3]), type = "response")
    pred.gp.bd <- predict(fit.gp.bd, as.character(xtest.icd[ , 3]), type = "response")

    res <- list(
        lin = list(fit = fit.lm, pred = pred.lm),
        ridge = list(fit = fit.r, pred = pred.r),
        lasso = list(fit = fit.lasso, pred = pred.lasso),
        svm.sp = list(fit = fit.ksvp.sp, pred = pred.ksvp.sp),
        svm.bd = list(fit = fit.ksvp.bd, pred = pred.ksvp.bd),
        gp.sp = list(fit = fit.gp.sp, pred = pred.gp.sp),
        gp.bd = list(fit = fit.gp.bd, pred = pred.gp.bd)
    )

    fname <- paste0("../data/comp_sim1_reg_", id, ".rds")
    saveRDS(res, fname)
} else if (mtd == 4) {
    estimateCodeFreq <- function(dx) {
        unq_codes = as.character(dx)
        freq_codes = aggregate(dx, list(unq_codes), length)
        colnames(freq_codes) = c("codes", "freq")
        cl_list =  vector("list", nrow(freq_codes))
        names(cl_list) = freq_codes$codes
        for (i in 1:nrow(freq_codes)) {
            len = nchar(freq_codes$codes[i])
            cl <- data.frame(matrix(ncol = 2, nrow = len))
            cl$X2 = freq_codes$freq[i]
            for(j in 1:len){
                cl$X1[j] = substr(freq_codes$codes[i], 1, j)
            }
            cl_list[[i]] = cl
        }

        sc_freq = do.call(rbind, cl_list)
        sc_freq = aggregate(sc_freq$X2, by = list(sc_freq$X1), sum)
        colnames(sc_freq) = c("subtrees", "freq")
        sc_freq
    }

    popChrIcd <- list("chr1" = c("A", "B", "AA", "BB", "BA", "AB"),
                      "chr2" = c("C", "D", "CC", "DD", "CD", "DC"),
                      "chr3" = c("E", "F", "EE", "FF", "EF", "FE"),
                      "chr4" = c("G", "H", "GG", "HH", "GH", "HG")
                      )

    code_freq = estimateCodeFreq(unlist(popChrIcd))

    ##
    kernFunction = function (code1, code2, code_freq) {
        if ((length(code1) > 0) & (length(code2) > 0)) {
            f1 = estimateCodeFreq(code1)
            f2 = estimateCodeFreq(code2)
            mfs = merge(f1, f2, by = "subtrees")
            if(nrow(mfs) == 0) {
                return(0)
            } else{
                mfsOverall = merge(mfs, code_freq, by = "subtrees")
                return(sum(mfsOverall$freq.x * mfsOverall$freq.y
                           * sapply(mfsOverall$subtrees, nchar)
                           ))
            }
        } else {
            return(0)
        }
    }

    kernelMatrix = function (p_code1, p_code2, code_freq) {
        pl1 = length(p_code1)
        pl2 = length(p_code2)

        K = matrix(NA, ncol = pl1, nrow = pl2)
        for (i in 1:pl1) {
            cat("i = ", i, "\n")
            for (j in i:pl2) {
                K [i, j] = kernFunction(p_code1[[i]]$chr1, p_code2[[j]]$chr1, code_freq) + kernFunction(p_code1[[i]]$chr2, p_code2[[j]]$chr2, code_freq) +
                    kernFunction(p_code1[[i]]$chr3, p_code2[[j]]$chr3, code_freq) + kernFunction(p_code1[[i]]$chr4, p_code2[[j]]$chr4, code_freq)
            }
        }
        K[lower.tri(K)] <- t(K)[lower.tri(K)]
        return(K)
    }

    dat <- readRDS("../data/sim_data_lin_corr_reg.rds")
    cvdat <- dat[[id]]

    ytrain <- cvdat$ytrain
    xtrain <- cvdat$xtrain

    xtest <- cvdat$xtest

    icdList <- c(cvdat$icdtrain, cvdat$icdtest)

    ## simMat <- kernelMatrix (icdList, icdList, code_freq)
    ## psim <- simMat / sqrt(outer(diag(simMat), diag(simMat), "*"))
    ## dsim <- outer(diag(psim), diag(psim), "+") - 2 * psim

    ## kernDistMat <- dsim ## this is the distance matrix computed from the similarity matrix

    ## saveRDS(kernDistMat, paste0("../data/kernDistMat_reg_", id, ".rds"))
    kernDistMat <- readRDS(paste0("../data/kernDistMat_reg_", id, ".rds"))

    source("gp_regression.R")

    res <- fitGPR (ytrain, xtrain, xtest, kernDistMat,
                   c(0, 0), diag(1, 2),
                   aphi = 0, bphi = 4,
                   niter = 10000, nburn = 5000, nthin = 5)

    fname <- paste0("../data/gp_sim1_reg_", id, ".rds")

    saveRDS(list(fit = res, sim = kernDistMat), fname)
} else {
    print("peace")
}
