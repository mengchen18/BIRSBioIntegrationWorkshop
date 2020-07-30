#' Prediction new samples using MBPCA based method
#' @param x a list of matrix/data.frame
#' @param train.var the training variable
#' @param test.var the variables to be predicted
#' @param validation in the validation step, whether the correlation of rows or columns should
#'   be calculated
#' @param ... other parameters passed to BioParallel::bplapply
#' 
pred_lm <- function(x, train.var, test.var, validation = c("row", "column")[1], ...) {
  
  xs <- lapply(x, function(x) as.data.frame(x[, train.var]))
  
  # prd <- BiocParallel::bplapply(names(xs), function(m) {
  prd <- lapply(names(xs), function(m) {
    
    # predicting matrix - m
    prdmat <- setdiff(names(xs), m)
    
    rmat <- lapply(prdmat, function(mm) {
      # using matrix mm to predict
      sapply(test.var, function(tv) {
        df <- data.frame(y = x[[mm]][, tv], xs[[mm]])
        mod <- lm(y ~ ., data = df)
        predict(mod, xs[[m]])
      })
      })
    names(rmat) <- prdmat
    rmat$mean <- Reduce("+", rmat)/length(rmat)
    rmatz <- lapply(rmat, function(x) t(scale(t(x))))
    rmat$meanZ <- Reduce("+", rmatz)/length(rmat)
    
    ## correlation of predicted variable (test set)
    pred.cor.test <- sapply(rmat, function(mm) {
      if (validation == "row") {
        r <- sapply(1:nrow(x[[m]]), function(i) {
          cor(unlist(x[[m]][i, test.var]), mm[i, ])
        })
      } else {
        r <- sapply(1:ncol(mm), function(i) {
          cor(unlist(x[[m]][, test.var[i]]), mm[, i])
        })
      }
      r
    })
    list(predicted = rmat, cor.test = pred.cor.test)  
  }, ...)
  names(prd) <- names(xs)
  prd
  
}






