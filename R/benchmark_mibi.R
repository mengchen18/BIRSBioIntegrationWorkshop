library(mogsa)
library(stringr)

source("R/mbpca_method.R")
source("R/lm_method.R")
load("Dat/mibiSCE.rda")

meta <- as.data.frame(colData(mibi.sce))
features <- as.data.frame(rowData(mibi.sce))
expr <- assay(mibi.sce)
expr <- expr[features$is_protein == 1, ]
rownames(expr) <- make.names(rownames(expr))

exprlist <- lapply(unique(meta$SampleID), function(v) {
  t(expr[, meta$SampleID == v])
})
names(exprlist) <- paste0("S", str_pad(unique(meta$SampleID), width = 2, pad = "0"))


s.train <- sample(1:ncol(exprlist[[1]]), size = 20)
s.test <- setdiff(1:ncol(exprlist[[1]]), s.train)

v <- pred_mbpca(
  x = exprlist, train.var = s.train, test.var = s.test, ncomp = 3, method = "globalScore", validation = "column"
  )
meanz <- lapply(v, function(x) x$cor.test[, "meanZ"])
boxplot(meanz)



########
v2 <- pred_lm( x = exprlist, train.var = s.train, test.var = s.test, validation = "column")


meanz <- lapply(v, function(x) x$cor.test[, "meanZ"])
meanz2 <- lapply(v2, function(x) x$cor.test[, "meanZ"])

boxplot(list(meanz[[1]], meanz2[[1]]))

