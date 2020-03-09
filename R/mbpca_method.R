#' Prediction new samples using MBPCA based method
#' @param x a list of matrix/data.frame
#' @param train.var the training variable
#' @param test.var the variables to be predicted
#' @param ncomp the number of components to calculate
#' @param method the deflation method used
#' @param validation in the validation step, whether the correlation of rows or columns should
#'   be calculated
#' @param ... other parameters passed to mogsa::mbpca
#' 
pred_mbpca <- function(x, train.var, test.var, ncomp = 2, method = "globalScore", 
                       validation = c("row", "column")[1], ...) { 
  
  x <- lapply(x, as.matrix)
  xs <- lapply(x, function(x) x[, train.var])
  rs <- mogsa::mbpca(x = xs, ncomp = ncomp, method = method, ...)
  nr <- rep(names(x), times = rs@tab.dim[1, ])
  
  lmats <- lapply(names(x), function(xx) {
    rs@loading[nr == xx, ]
  })
  names(lmats) <- names(x)
  
  prd <- BiocParallel::bplapply(names(lmats), function(m) {
    # predicting matrix - m
    prdmat <- setdiff(names(lmats), m)
    rmat <- lapply(prdmat, function(mm) {
      q <- t(x[[mm]]) %*% lmats[[mm]]
      mm <- q %*%  t(lmats[[m]])
      t(mm)
    })
    names(rmat) <- prdmat
    rmat$mean <- Reduce("+", rmat)/length(rmat)
    rmatz <- lapply(rmat, function(x) t(scale(t(x))))
    rmat$meanZ <- Reduce("+", rmatz)/length(rmat)
    
    ## correlation of predicted variable (train set)
    pred.cor.train <- sapply(rmat, function(mm) {
      if (validation == "row") {
        r <- sapply(1:nrow(x[[m]]), function(i) {
          cor(x[[m]][i, train.var], mm[i, train.var])
        })
      } else {
        r <- sapply(1:length(train.var), function(i) {
          cor(x[[m]][, train.var[i]], mm[, train.var[i]])
        })
      }
      r
    })
    
    ## correlation of predicted variable (test set)
    pred.cor.test <- sapply(rmat, function(mm) {
      
      if (validation == "row") {
        r <- sapply(1:nrow(x[[m]]), function(i) {
          cor(x[[m]][i, test.var], mm[i, test.var])
        })
        # r <- sapply(1:nrow(x[[m]]), function(i) {
        #   cor(unlist(x[[m]][i, test.var]), mm[i, ])
        # })
      } else {
        r <- sapply(1:length(test.var), function(i) {
          cor(x[[m]][, test.var[i]], mm[, test.var[i]])
        })
      }
      r
      
    })
    
    list(predicted = rmat, cor.train = pred.cor.train, cor.test = pred.cor.test)
  })
  names(prd) <- names(lmats)
  prd
  }



