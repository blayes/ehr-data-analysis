cmdArgs <- commandArgs(trailingOnly = TRUE)

mtd <- as.numeric(cmdArgs[1])
m <- as.numeric(cmdArgs[2])
ii <- ceiling(m / 10)
cc <- m - (ii-1) * 10

if (mtd == 1) {
    library(matrixStats)
    library(kernlab)
    library(glmnet)
    library(ranger)

    dat <- readRDS(paste0("../data/sim_data_lin_corr_phecode", ii, ".rds"))
    cvdat <- dat[[cc]]

    ytrain <- cvdat$ytrain
    xtrain <- cvdat$xtrain
    xtrain.phe <- cvdat$xtrain.phe
    xtrain.feat <- cvdat$xtrain.feat
    ptrain <- cvdat$ptrain

    ytest <- cvdat$ytest
    xtest <- cvdat$xtest
    xtest.phe <- cvdat$xtest.phe
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
    pred.lasso <- predict(fit.lasso, newx = xtest.phe, s = "lambda.min", type = "response")
    
    phecount.train = rowSums(xtrain.feat[, -c(1:2)])
    phecount.test = rowSums(xtest.feat[, -c(1:2)])
    
    # lin with log(1+count) only 
    dftrain.c = cbind.data.frame(y = ytrain, xtrain.feat[, c(1:2)], lc = log(1+phecount.train))
    dftest.c = cbind.data.frame(y = ytest, xtest.feat[, c(1:2)], lc = log(1 + phecount.test))
    fit.lg.c <- glm(y ~ 0 + . ,family = "binomial", data = dftrain.c)
    pred.lg.c <- predict(fit.lg.c, dftest.c, type = "response")

    res <- list(
        logis = list(fit = fit.lg, pred = pred.lg),
        ridge = list(fit = fit.r, pred = pred.r),
        lasso = list(fit = fit.lasso, pred = pred.lasso),
        logis.c = list(fit = fit.lg.c, pred = pred.lg.c)
    )

    fname <- paste0("../data/comp_sim1_", ii, "_", cc, ".rds")
    saveRDS(res, fname)
} else if (mtd == 3) {
    library(matrixStats)
    library(kernlab)
    library(glmnet)
    library(ranger)

    dat <- readRDS(paste0("../data/sim_data_lin_corr_reg_", ii, ".rds"))
    cvdat <- dat[[cc]]

    ytrain <- cvdat$ytrain
    xtrain <- cvdat$xtrain
    xtrain.phe <- cvdat$xtrain.phe
    xtrain.feat <- cvdat$xtrain.feat

    ytest <- cvdat$ytest
    xtest <- cvdat$xtest
    xtest.phe <- cvdat$xtest.phe
    xtest.feat <- cvdat$xtest.feat

    ## lin regression
    
    phecount.train = rowSums(xtrain.feat[, -c(1:2)])
    phecount.test = rowSums(xtest.feat[, -c(1:2)])
    
    dftrain = cbind.data.frame(y = ytrain, xtrain.feat[, c(1:2)], lc = log(1+phecount.train))
    dftest = cbind.data.frame(y = ytest, xtest.feat[, c(1:2)], lc = log(1 + phecount.test))
    fit.lm <- lm(y ~ 0 + ., data = dftrain)
    pred.lm <- predict(fit.lm, dftest)



    ## ridge regression
    fit.r <- cv.glmnet(xtrain.feat, ytrain, alpha = 0)
    pred.r <- predict(fit.r, newx = xtest.feat, s = "lambda.min", type = "response")

    ## lasso regression
    fit.lasso <- cv.glmnet(xtrain.feat, ytrain, alpha = 1)
    pred.lasso <- predict(fit.lasso, newx = xtest.feat, s = "lambda.min", type = "response")
 
    res <- list(
        lin = list(fit = fit.lm, pred = pred.lm),
        ridge = list(fit = fit.r, pred = pred.r),
        lasso = list(fit = fit.lasso, pred = pred.lasso)
    )

    fname <- paste0("../data/comp_sim1_reg_", ii, "_", cc, ".rds")
    saveRDS(res, fname)
} else {
    print("peace")
}
