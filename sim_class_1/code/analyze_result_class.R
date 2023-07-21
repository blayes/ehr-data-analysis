library(glmnet)
library(matrixStats)
library(xtable)

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

test <- readRDS("../data/sim_data_lin_corr.rds")

logClass <- matrix(NA, 10, 4)
lassoClass <- matrix(NA, 10, 4)
ridClass <- matrix(NA, 10, 4)
gpSpClass <- matrix(NA, 10, 4)
gpBdClass <- matrix(NA, 10, 4)
svmSpClass <- matrix(NA, 10, 4)
svmBdClass <- matrix(NA, 10, 4)

colnames(logClass) <- colnames(lassoClass) <- colnames(ridClass) <- colnames(gpSpClass) <- colnames(gpBdClass) <-
    colnames(svmSpClass) <- colnames(svmBdClass) <- c("auc", "er", "sens", "spec")
for (cc in 1:10) {
    res1 <- readRDS(paste0("../data/comp_sim1_", cc, ".rds"))
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

    svmSpClass[cc, "auc"] <- summ(as.numeric(res1$svm.sp$pred), test[[cc]]$ytest)$auc
    svmSpClass[cc, "er"] <- summ(as.numeric(res1$svm.sp$pred), test[[cc]]$ytest)$rate
    svmSpClass[cc, "sens"] <- summ(as.numeric(res1$svm.sp$pred), test[[cc]]$ytest)$sens
    svmSpClass[cc, "spec"] <- summ(as.numeric(res1$svm.sp$pred), test[[cc]]$ytest)$spec

    svmBdClass[cc, "auc"] <- summ(as.numeric(res1$svm.bd$pred), test[[cc]]$ytest)$auc
    svmBdClass[cc, "er"] <- summ(as.numeric(res1$svm.bd$pred), test[[cc]]$ytest)$rate
    svmBdClass[cc, "sens"] <- summ(as.numeric(res1$svm.bd$pred), test[[cc]]$ytest)$sens
    svmBdClass[cc, "spec"] <- summ(as.numeric(res1$svm.bd$pred), test[[cc]]$ytest)$spec

    gpSpClass[cc, "auc"] <- summ(as.numeric(res1$gp.sp$pred), test[[cc]]$ytest)$auc
    gpSpClass[cc, "er"] <- summ(as.numeric(res1$gp.sp$pred), test[[cc]]$ytest)$rate
    gpSpClass[cc, "sens"] <- summ(as.numeric(res1$gp.sp$pred), test[[cc]]$ytest)$sens
    gpSpClass[cc, "spec"] <- summ(as.numeric(res1$gp.sp$pred), test[[cc]]$ytest)$spec

    gpBdClass[cc, "auc"] <- summ(as.numeric(res1$gp.bd$pred), test[[cc]]$ytest)$auc
    gpBdClass[cc, "er"] <- summ(as.numeric(res1$gp.bd$pred), test[[cc]]$ytest)$rate
    gpBdClass[cc, "sens"] <- summ(as.numeric(res1$gp.bd$pred), test[[cc]]$ytest)$sens
    gpBdClass[cc, "spec"] <- summ(as.numeric(res1$gp.bd$pred), test[[cc]]$ytest)$spec        
}


resSummEst <- rbind.data.frame(colMeans(logClass),
                               colMeans(ridClass),
                               colMeans(lassoClass),
                               colMeans(svmSpClass),
                               colMeans(svmBdClass),
                               colMeans(gpSpClass),
                               colMeans(gpBdClass)
                               )
rownames(resSummEst) <- c("logistic", "ridge", "lasso", "svm.sp", "svm.bd", "gp.sp", "gp.bd")
colnames(resSummEst) <- colnames(lassoClass)

resSummSd <- rbind.data.frame(colSds(logClass),
                              colSds(ridClass),
                              colSds(lassoClass),
                              colSds(svmSpClass),
                              colSds(svmBdClass),
                              colSds(gpSpClass),
                              colSds(gpBdClass)
                              )
rownames(resSummSd) <- c("logistic", "ridge", "lasso", "svm.sp", "svm.bd", "gp.sp", "gp.bd")
colnames(resSummSd) <- colnames(lassoClass)

icdGpClass <- matrix(NA, 10, 4)
icdLogClass <- matrix(NA, 10, 4)
icdRidClass <- matrix(NA, 10, 4)
icdLassoClass <- matrix(NA, 10, 4)

colnames(icdGpClass) <- colnames(icdLogClass) <- colnames(icdRidClass) <- colnames(icdLassoClass) <- c("auc", "er", "sens", "spec")

for (cc in 1:10) {
    res2 <- readRDS(paste0("../data/gp_sim1_", cc, ".rds"))
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

}

resIcdEst <- rbind.data.frame(colMeans(icdGpClass), 
                              colMeans(icdLogClass), 
                              colMeans(icdRidClass), 
                              colMeans(icdLassoClass))

rownames(resIcdEst) <- c("icd.gp", "icd.logistic", "icd.ridge", "icd.lasso")
colnames(resIcdEst) <- colnames(lassoClass)

resIcdSd <- rbind.data.frame(colSds(icdGpClass), 
                             colSds(icdLogClass), 
                             colSds(icdRidClass), 
                             colSds(icdLassoClass))

rownames(resIcdSd) <- c("icd.gp", "icd.logistic", "icd.ridge", "icd.lasso")
colnames(resIcdSd) <- colnames(lassoClass)

xtable(round(rbind(resIcdEst, resSummEst), 2))
