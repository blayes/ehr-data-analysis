library(glmnet)
library(matrixStats)
library(xtable)
library(ranger)
library(abind)

summ <- function (preds, ytest) {
    library(pROC)

    roc.p <- roc(ytest, preds)
    auc.p <- as.numeric(auc(roc.p))
    p <- coords(roc.p, "best", transpose = TRUE)[1]
    spec <- coords(roc.p, "best", transpose = TRUE)[2]
    sens <- coords(roc.p, "best", transpose = TRUE)[3]
    rate <- mean(as.numeric(preds > p) == ytest)

    list(p = p,
         auc = auc.p,
         spec = spec,
         sens = sens,
         rate = rate
         )
}

summArray <- array()
for (ii in 1:15) {
    test <- readRDS(paste0("../data/sim_data_lin_corr_phecode", ii, ".rds"))
    
    logClass <- matrix(NA, 10, 4)
    lassoClass <- matrix(NA, 10, 4)
    ridClass <- matrix(NA, 10, 4)
    logcClass <- matrix(NA, 10, 4)
    lassocClass <- matrix(NA, 10, 4)
    ridcClass <- matrix(NA, 10, 4)
    
    colnames(logClass) <- colnames(lassoClass) <- colnames(ridClass) <- 
      colnames(logcClass) <- colnames(lassocClass) <- colnames(ridcClass) <- 
      colnames(gpSpClass) <- colnames(gpBdClass) <-
        colnames(svmSpClass) <- colnames(svmBdClass) <- colnames(rfClass) <- 
        c("auc", "er", "sens", "spec")
    for (cc in 1:10) {
        res1 <- readRDS(paste0("../data/comp_sim1_", ii, "_", cc, ".rds"))
        logClass[cc, "auc"] <- summ(res1$logis$pred, test[[cc]]$ytest)$auc
        logClass[cc, "er"] <- summ(res1$logis$pred, test[[cc]]$ytest)$rate
        logClass[cc, "sens"] <- summ(res1$logis$pred, test[[cc]]$ytest)$sens
        logClass[cc, "spec"] <- summ(res1$logis$pred, test[[cc]]$ytest)$spec

        ridClass[cc, "auc"] <- summ(as.numeric(res1$ridge$pred), test[[cc]]$ytest)$auc
        ridClass[cc, "er"] <- summ(as.numeric(res1$ridge$pred), test[[cc]]$ytest)$rate
        ridClass[cc, "sens"] <- summ(as.numeric(res1$ridge$pred), test[[cc]]$ytest)$sens
        ridClass[cc, "spec"] <- summ(as.numeric(res1$ridge$pred), test[[cc]]$ytest)$spec

        lassoClass[cc, "auc"] <- summ(as.numeric(res1$lasso$pred), test[[cc]]$ytest)$auc
        lassoClass[cc, "er"] <- summ(as.numeric(res1$lasso$pred), test[[cc]]$ytest)$rate
        lassoClass[cc, "sens"] <- summ(as.numeric(res1$lasso$pred), test[[cc]]$ytest)$sens
        lassoClass[cc, "spec"] <- summ(as.numeric(res1$lasso$pred), test[[cc]]$ytest)$spec

        logcClass[cc, "auc"] <- summ(res1$logis.c$pred, test[[cc]]$ytest)$auc
        logcClass[cc, "er"] <- summ(res1$logis.c$pred, test[[cc]]$ytest)$rate
        logcClass[cc, "sens"] <- summ(res1$logis.c$pred, test[[cc]]$ytest)$sens
        logcClass[cc, "spec"] <- summ(res1$logis.c$pred, test[[cc]]$ytest)$spec
        
        ridcClass[cc, "auc"] <- summ(as.numeric(res1$ridge.c$pred), test[[cc]]$ytest)$auc
        ridcClass[cc, "er"] <- summ(as.numeric(res1$ridge.c$pred), test[[cc]]$ytest)$rate
        ridcClass[cc, "sens"] <- summ(as.numeric(res1$ridge.c$pred), test[[cc]]$ytest)$sens
        ridcClass[cc, "spec"] <- summ(as.numeric(res1$ridge.c$pred), test[[cc]]$ytest)$spec
        
        lassocClass[cc, "auc"] <- summ(as.numeric(res1$lasso.c$pred), test[[cc]]$ytest)$auc
        lassocClass[cc, "er"] <- summ(as.numeric(res1$lasso.c$pred), test[[cc]]$ytest)$rate
        lassocClass[cc, "sens"] <- summ(as.numeric(res1$lasso.c$pred), test[[cc]]$ytest)$sens
        lassocClass[cc, "spec"] <- summ(as.numeric(res1$lasso.c$pred), test[[cc]]$ytest)$spec
    }
    
    
    resSummEst <- rbind.data.frame(colMeans(logClass),
                                   colMeans(ridClass),
¯                                   colMeans(lassoClass),
                                   colMeans(logcClass),
                                   colMeans(ridcClass),
                                   colMeans(lassocClass)
    )
    rownames(resSummEst) <- c("logistic", "ridge", "lasso",
                              "logistic.c", "ridge.c", "lasso.c"
                              )
    colnames(resSummEst) <- colnames(lassocClass)
    
    resSummSd <- rbind.data.frame(colSds(logClass),
                                  colSds(ridClass),
                                  colSds(lassoClass),
                                  colSds(logcClass),
                                  colSds(ridcClass),
                                  colSds(lassocClass)
    )
    rownames(resSummSd) <- c("logistic", "ridge", "lasso",
                             "logistic.c", "ridge.c", "lasso.c"
                             )
    colnames(resSummSd) <- colnames(lassocClass)
    
    icdGpClass <- matrix(NA, 10, 4)
    icdLogClass <- matrix(NA, 10, 4)
    icdRidClass <- matrix(NA, 10, 4)
    icdLassoClass <- matrix(NA, 10, 4)
    icdRfClass <- matrix(NA, 10, 4)
    
    icdLogcClass <- matrix(NA, 10, 4)
    icdRidcClass <- matrix(NA, 10, 4)
    icdLassocClass <- matrix(NA, 10, 4)
    
    colnames(icdGpClass) <- 
      colnames(icdLogcClass) <- colnames(icdRidcClass) <- colnames(icdLassocClass) <- colnames(icdRfClass) <-  c("auc", "er", "sens", "spec")
    
    for (cc in 1:10) {
        res2 <- readRDS(paste0("../data/gp_sim1_", ii, "_", cc, ".rds"))
        icdGpClass[cc, "auc"] <- summ(colMeans(res2$fit$prob), test[[cc]]$ytest)$auc
        icdGpClass[cc, "er"] <- summ(colMeans(res2$fit$prob), test[[cc]]$ytest)$rate
        icdGpClass[cc, "sens"] <- summ(colMeans(res2$fit$prob), test[[cc]]$ytest)$sens
        icdGpClass[cc, "spec"] <- summ(colMeans(res2$fit$prob), test[[cc]]$ytest)$spec
        
        simMat <- res2$sim
        
        ksvd <- svd(simMat, 10, 10)
        U <- ksvd$u[ , 1:10]
        D <- sqrt(ksvd$d[1:10])
        featMat <- U * matrix(D, nrow(simMat), 10, byrow = TRUE)
        colnames(featMat) <- paste0("U", 1:10)
        
        dftrain <- cbind.data.frame(y = test[[cc]]$ytrain, test[[cc]]$xtrain, featMat[1:1000, ])
        dftest <- cbind.data.frame(y = test[[cc]]$ytest, test[[cc]]$xtest, featMat[-(1:1000), ])
        
        ## logistic regression
        fit.lg <- glm(y ~ 0 + ., family = "binomial", data = dftrain)
        pred.lg <- predict(fit.lg, dftest, type = "response")
        icdLogClass[cc, "auc"] <- summ(pred.lg, test[[cc]]$ytest)$auc
        icdLogClass[cc, "er"] <- summ(pred.lg, test[[cc]]$ytest)$rate
        icdLogClass[cc, "sens"] <- summ(pred.lg, test[[cc]]$ytest)$sens
        icdLogClass[cc, "spec"] <- summ(pred.lg, test[[cc]]$ytest)$spec
        
        ## ridge regression
        fit.r <- cv.glmnet(cbind(test[[cc]]$xtrain, featMat[1:1000, ]), test[[cc]]$ytrain, alpha = 0, family = "binomial")
        pred.r <- drop(predict(fit.r, newx = cbind(test[[cc]]$xtest, featMat[-(1:1000), ]), s = "lambda.min", type = "response"))
        icdRidClass[cc, "auc"] <- summ(pred.r, test[[cc]]$ytest)$auc
        icdRidClass[cc, "er"] <- summ(pred.r, test[[cc]]$ytest)$rate
        icdRidClass[cc, "sens"] <- summ(pred.r, test[[cc]]$ytest)$sens
        icdRidClass[cc, "spec"] <- summ(pred.r, test[[cc]]$ytest)$spec
        
        ## lasso regression
        fit.lasso <- cv.glmnet(cbind(test[[cc]]$xtrain, featMat[1:1000, ]), test[[cc]]$ytrain, alpha = 1, family = "binomial")
        pred.lasso <- drop(predict(fit.lasso, newx = cbind(test[[cc]]$xtest, featMat[-(1:1000), ]), s = "lambda.min", type = "response"))
        icdLassoClass[cc, "auc"] <- summ(pred.lasso, test[[cc]]$ytest)$auc
        icdLassoClass[cc, "er"] <- summ(pred.lasso, test[[cc]]$ytest)$rate
        icdLassoClass[cc, "sens"] <- summ(pred.lasso, test[[cc]]$ytest)$sens
        icdLassoClass[cc, "spec"] <- summ(pred.lasso, test[[cc]]$ytest)$spec
        
        ## random forest
        fit.rf <- ranger(y ~ . , data = dftrain, classification = TRUE, importance = 'impurity', respect.unordered.factors = TRUE)
        pred.rf <- predict(fit.rf, dftest)$predictions
        icdRfClass[cc, "auc"] <-  summ(pred.rf, test[[cc]]$ytest)$auc
        icdRfClass[cc, "er"] <-   summ(pred.rf, test[[cc]]$ytest)$rate
        icdRfClass[cc, "sens"] <- summ(pred.rf, test[[cc]]$ytest)$sens
        icdRfClass[cc, "spec"] <- summ(pred.rf, test[[cc]]$ytest)$spec
    }
    
    resIcdEst <- rbind.data.frame(colMeans(icdGpClass), 
                                  colMeans(icdLogClass), 
                                  colMeans(icdRidClass), 
                                  colMeans(icdLassoClass),
                                  colMeans(icdRfClass))
    
    rownames(resIcdEst) <- c("icd.gp", "icd.logistic", "icd.ridge", "icd.lasso", "icd.rf")
    colnames(resIcdEst) <- colnames(lassoClass)
    
    resIcdSd <- rbind.data.frame(colSds(icdGpClass), 
                                 colSds(icdLogClass), 
                                 colSds(icdRidClass), 
                                 colSds(icdLassoClass),
                                 colSds(icdRfClass))
    
    rownames(resIcdSd) <- c("icd.gp", "icd.logistic", "icd.ridge", "icd.lasso", "icd.rf")
    colnames(resIcdSd) <- colnames(lassoClass)
    
    if (ii == 1) summArray <- rbind(resIcdEst, resSummEst)
    else summArray <- abind(summArray, rbind(resIcdEst, resSummEst), along = 3)
    
}

xtable(rowMeans(summArray, dims = 2))
