

rm(list = ls())

library(glmnet)
library(matrixStats)
library(xtable)

summArray <- array()

for (ii in 1:15) {
  test <- readRDS(paste0("../data/sim_data_lin_corr_reg_", ii, ".rds"))

  linReg <- matrix(NA, 10, 6)
  lassoReg <- matrix(NA, 10, 6)
  ridReg <- matrix(NA, 10, 6)
  
  colnames(linReg) <- colnames(lassoReg) <- colnames(ridReg) <- 
    c("mse.b", "mse.tau", "mse.pred",
      "cov.b", "cov.tau", "cov.pred")
  
  for (cc in 1:10) {
    res1 <- readRDS(paste0("../data/comp_sim1_reg_1_", cc, ".rds"))
    
    linReg[cc, "mse.pred"] <- sqrt(mean((res1$lin$pred - test[[cc]]$ytest)^2))
    linReg[cc, "mse.tau"] <- abs(mean(residuals(res1$lin$fit)^2) - 0.01)
    linReg[cc, "mse.b"] <- sqrt(sum((coef(res1$lin$fit)[1:2] - c(0.1, 0.2))^2))
    
    lassoReg[cc, "mse.pred"] <- sqrt(mean((drop(res1$lasso$pred) - test[[cc]]$ytest)^2))
    lassoReg[cc, "mse.tau"] <- abs(mean((drop(res1$lasso$pred) - test[[cc]]$ytest)^2) - 0.01)
    lassoReg[cc, "mse.b"] <- sqrt(sum((coef(res1$lasso$fit)[c(1, 3)] - c(0.1, 0.2))^2))
    
    ridReg[cc, "mse.pred"] <- sqrt(mean((drop(res1$ridge$pred) - test[[cc]]$ytest)^2))
    ridReg[cc, "mse.tau"] <- abs(mean((drop(res1$ridge$pred) - test[[cc]]$ytest)^2) - 0.01)
    ridReg[cc, "mse.b"] <- sqrt(sum((coef(res1$ridge$fit)[c(1, 3)] - c(0.1, 0.2))^2))
    
    xtest.feat = test[[cc]]$xtest.feat
    phecount.test = rowSums(xtest.feat[, -c(1:2)])
    xtest.count = cbind.data.frame(xtest.feat[, c(1:2)], lc = log(1+phecount.train))
    tmp <- predict(res1$lin$fit, xtest.count, interval = "confidence")
    linReg[cc, "cov.pred"] <- mean((test[[cc]]$ytest <= tmp[ , 3]) & (test[[cc]]$ytest >= tmp[ , 2]))
    tmp <- confint(res1$lin$fit, parm = c("Int", "Slope"))
    linReg[cc, "cov.b"] <- mean((c(0.1, 0.2) <= tmp[ , 2]) & (c(0.1, 0.2) >= tmp[ , 1]))
  }
  
  resSummEst <- rbind(colMeans(linReg^2, na.rm = TRUE),
                      colMeans(ridReg^2, na.rm = TRUE),
                      colMeans(lassoReg^2, na.rm = TRUE)
  )
  rownames(resSummEst) <- c("linear", "ridge", "lasso")
  
  resSummSd <- rbind(colSds(linReg, na.rm = TRUE),
                     colSds(ridReg, na.rm = TRUE),
                     colSds(lassoReg, na.rm = TRUE)
  )
  names(resSummSd) <- c("linear", "ridge", "lasso")

  
  icdGpReg <- matrix(NA, 10, 6)
  icdLinReg <- matrix(NA, 10, 6)
  icdRidReg <- matrix(NA, 10, 6)
  icdLassoReg <- matrix(NA, 10, 6)
  icdRfReg<- matrix(NA, 10, 6)
  
  colnames(icdLassoReg) <- colnames(icdRidReg) <- colnames(icdLinReg) <- colnames(icdGpReg) <- colnames(icdRfReg) <- 
    c("mse.b", "mse.tau", "mse.pred",
      "cov.b", "cov.tau", "cov.pred")
  
  
  for (cc in 1:10) {
    res2 <- readRDS(paste0("../data/gp_sim1_reg_", ii, "_", cc, ".rds"))
    icdGpReg[cc, "mse.pred"] <- sqrt(mean((colMeans(do.call(rbind, res2$fit$pred)) -  test[[cc]]$ytest)^2))
    icdGpReg[cc, "mse.tau"] <- (mean(abs(res2$fit$tau - 0.01)))
    icdGpReg[cc, "mse.b"] <- sqrt(mean((colMeans(do.call(rbind, res2$fit$beta)) -  c(0.1, 0.2))^2))
    
    bb <- colQuantiles(do.call(rbind, res2$fit$beta), probs = c(0.025, 0.975))
    yy <- colQuantiles(do.call(rbind, res2$fit$pred), probs = c(0.025, 0.975))
    tt <- quantile(res2$fit$tau, probs = c(0.025, 0.975))
    
    icdGpReg[cc, "cov.pred"] <- mean((test[[cc]]$ytest <= yy[, 2]) & (test[[cc]]$ytest >= yy[ , 1]))
    icdGpReg[cc, "cov.tau"] <- 0.01 <= tt[2] & 0.01 >= tt[1]
    icdGpReg[cc, "cov.b"] <- mean((c(0.1, 0.2) <= bb[, 2]) & (c(0.1, 0.2) >= bb[ , 1]))
    
    simMat <- res2$sim
    
    ksvd <- svd(simMat, 10, 10)
    U <- ksvd$u[ , 1:10]
    D <- sqrt(ksvd$d[1:10])
    featMat <- U * matrix(D, nrow(simMat), 10, byrow = TRUE)
    colnames(featMat) <- paste0("U", 1:10)
    
    dftrain <- cbind.data.frame(y = test[[cc]]$ytrain, test[[cc]]$xtrain, featMat[1:1000, ])
    dftest <- cbind.data.frame(y = test[[cc]]$ytest, test[[cc]]$xtest, featMat[-(1:1000), ])
    
    ## lin regression
    fit.lin <- lm(y ~ 0 + ., data = dftrain)
    pred.lin <- predict(fit.lin, dftest, type = "response")
    
    icdLinReg[cc, "mse.pred"] <- sqrt(mean((pred.lin - test[[cc]]$ytest)^2))
    icdLinReg[cc, "mse.tau"] <- abs(mean((pred.lin - test[[cc]]$ytest)^2) - 0.01)
    icdLinReg[cc, "mse.b"] <- mean((coef(fit.lin)[1:2] - c(0.1, 0.2))^2)
    
    tmp <- predict(fit.lin, dftest, interval = "confidence")
    icdLinReg[cc, "cov.pred"] <- mean((test[[cc]]$ytest <= tmp[ , 3]) & (test[[cc]]$ytest >= tmp[ , 2]))
    tmp <- confint(fit.lin, parm = c("Int", "Slope"))
    icdLinReg[cc, "cov.b"] <- mean((c(0.1, 0.2) <= tmp[ , 2]) & (c(0.1, 0.2) >= tmp[ , 1]))
    
    ## ridge regression
    fit.r <- cv.glmnet(cbind(test[[cc]]$xtrain, featMat[1:1000, ]), test[[cc]]$ytrain, alpha = 0)
    pred.r <- drop(predict(fit.r, newx = cbind(test[[cc]]$xtest, featMat[-(1:1000), ]), s = "lambda.min", type = "response"))
    icdRidReg[cc, "mse.pred"] <- sqrt(mean((drop(pred.r) - test[[cc]]$ytest)^2))
    icdRidReg[cc, "mse.tau"] <- abs(mean((drop(pred.r) - test[[cc]]$ytest)^2) - 0.01)
    icdRidReg[cc, "mse.b"] <- sqrt(sum((coef(fit.r)[c(1, 3)] - c(0.1, 0.2))^2))
    
    ## lasso regression
    fit.lasso <- cv.glmnet(cbind(test[[cc]]$xtrain, featMat[1:1000, ]), test[[cc]]$ytrain, alpha = 1)
    pred.lasso <- drop(predict(fit.lasso, newx = cbind(test[[cc]]$xtest, featMat[-(1:1000), ]), s = "lambda.min", type = "response"))
    
    icdLassoReg[cc, "mse.pred"] <- sqrt(mean((drop(pred.lasso) - test[[cc]]$ytest)^2))
    icdLassoReg[cc, "mse.tau"] <- abs(mean((drop(pred.lasso) - test[[cc]]$ytest)^2) - 0.01)
    icdLassoReg[cc, "mse.b"] <- sqrt(sum((coef(fit.lasso)[c(1, 3)] - c(0.1, 0.2))^2))
    
    ## random forest
    fit.rf <- ranger(y ~ . , data = dftrain, classification = TRUE, importance = 'impurity', respect.unordered.factors = TRUE)
    pred.rf <- predict(fit.rf, dftest)$predictions
    
    icdRfReg[cc, "mse.pred"] <- sqrt(mean((drop(pred.rf) - test[[cc]]$ytest)^2))
    icdRfReg[cc, "mse.tau"] <- abs(mean((drop(pred.rf) - test[[cc]]$ytest)^2) - 0.01)
    icdRfReg[cc, "mse.b"] <- sqrt(sum((coef(fit.rf)[c(1, 3)] - c(0.1, 0.2))^2))
  }
  
  resIcdEst <- rbind(colMeans(icdGpReg, na.rm = TRUE),
                     colMeans(icdLinReg, na.rm = TRUE),
                     colMeans(icdRidReg, na.rm = TRUE),
                     colMeans(icdLassoReg, na.rm = TRUE),
                     colMeans(icdRfReg, na.rm = TRUE))
  
  rownames(resIcdEst) <- c("icd.gp", "icd.linear", "icd.ridge", "icd.lasso", "icd.rf")
  
  resIcdSd <- rbind(colSds(icdGpReg, na.rm = TRUE),
                colSds(icdLinReg, na.rm = TRUE),
                colSds(icdRidReg, na.rm = TRUE),
                colSds(icdLassoReg, na.rm = TRUE),
                colSds(icdRfReg, na.rm = TRUE))
  
  rownames(resIcdSd) <- c("icd.gp", "icd.linear", "icd.ridge", "icd.lasso", "icd.rf")
  
  if (ii == 1) summArray <- rbind(resIcdEst, resSummEst)
  else summArray <- abind(summArray, rbind(resIcdEst, resSummEst), along = 3)
  
}

xtable(rowMeans(summArray, dims = 2))
