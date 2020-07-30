library(mogsa)
library(stringr)
library(SummarizedExperiment)

source("R/mbpca_method.R")
source("R/lm_method.R")
load("Dat/masstagSCE.rda")


dat <- livecells.sce
meta <- as.data.frame(colData(dat), stringsAsFactors = FALSE)
i <- which(meta$Clinical.Subtype == "TN")

features <- as.data.frame(rowData(dat))
expr <- assay(dat)
rownames(expr) <- make.names(rownames(expr))

pt <- as.character(unique(meta$patient_id.x[i]))
exprlist <- lapply(pt, function(v) {
  scale(t(expr[, meta$patient_id.x == v]), center = TRUE, scale = TRUE)
})
names(exprlist) <- pt
sapply(exprlist, nrow)


######
s.train <- sample(1:ncol(exprlist[[1]]), size = 25)
s.test <- setdiff(1:ncol(exprlist[[1]]), s.train)

v <- pred_mbpca(
  x = exprlist, train.var = s.train, test.var = s.test, ncomp = 7, 
  method = "blockLoading", center = FALSE, validation = "column"
)
meanz <- lapply(v, function(x) x$cor.test[, "meanZ"])
boxplot(meanz)


########
v2 <- pred_lm( x = exprlist, train.var = s.train, test.var = s.test, validation = "column" )
meanz2 <- lapply(v2, function(x) x$cor.test[, "meanZ"])
boxplot(meanz2)

boxplot(c(meanz, meanz2))
