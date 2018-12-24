# 
# nested cross validation for reGAIN feature selection in flu vaccine
#
rm(list = ls())

library(Rinbix)
library(caret)

# load 5000 filtered genes obtained by Coeffecient of Variation
load("B_ex.cv.fltr.RData")
# load phenotype
load("bay.fc.log.resid.RData")


# library(parallel)
# library(foreach)
# library(doParallel)
# 
# # Calculate the number of cores
# no_cores <- 20
# 
# # Initiate cluster
# cl <- makeCluster(no_cores)
# registerDoParallel(cl)


baylor_features <- list()
# create outer folds
outer_folds <- caret::createFolds(bay.fc.log.resid, 10, list = FALSE)
# create inner folds
for (i in 1:10) {
 inner.features <- list()
 inner_folds <- caret::createFolds(bay.fc.log.resid[outer_folds!=i], 10, list = TRUE)
 for (j in 1:10) {
   fold_idx <- which(outer_folds != i)[-inner_folds[[j]]]
   # reGAIN
   baylor_expr <- B_ex.cv.fltr[fold_idx,]
   class <- bay.fc.log.resid[fold_idx]
   regain_matrix <- data.frame(baylor_expr, class)
   rownames(regain_matrix) <- NULL
   baylor_regain_resid <- regainInbix(regain_matrix)$reGAIN
   #baylor_regain <- regainParallel(regain_matrix, regressionFamily = "gaussian", numCores = 4)
   idx <- NULL
   regain_sorted <- unique(sort(abs(baylor_regain), TRUE))
   for (k in 1:2000){
     idx <- c(idx, which(baylor_regain==regain_sorted[k], arr.ind = TRUE))
   }
   inner.features[[j]] <- colnames(baylor_regain)[unique(as.vector(idx))[1:200]]
 }
 baylor_features[[i]] <- inner.features
}
save(baylor_regain_resid, file = "baylor_regain_resid.RData")

# Finish
#stopCluster(cl)
