library(mogsa)
library(omicade4)
data("NCI60_4arrays")

source("R/mbpca_method.R")
source("R/lm_method.R")

### using 40 in training and 20 for testing

ss <- sample(1:60, size = 40)
v <- pred_mbpca(x = NCI60_4arrays, train.var = ss, test.var = setdiff(1:60, ss), ncomp = 5, method = "globalScore")
v2 <- pred_lm(x = NCI60_4arrays, train.var = ss, test.var = setdiff(1:60, ss))


t1 <- lapply(v, function(x) x$"cor.train"[, "meanZ"])
t2 <- lapply(v, function(x) x$"cor.test"[, "meanZ"])
l <- mapply(function(x,y) {
  list(MBPCA = x, OLS = y)
}, t1, t2, SIMPLIFY = FALSE)
l <- unlist(l, recursive = FALSE)
png("Fig/benchmark_NCI60_train40_test20.png", width = 12, height = 12, units = "cm", res = 300)
par(mar = c(8, 4, 1, 1))
boxplot(l, las = 2, ylab = "Correlation (predicted vs measured)")
dev.off()

####
ss <- sample(1:60, size = 20)
v <- pred_mbpca(x = NCI60_4arrays, train.var = ss, test.var = setdiff(1:60, ss), ncomp = 5, method = "globalScore")
v2 <- pred_lm(x = NCI60_4arrays, train.var = ss, test.var = setdiff(1:60, ss))


t1 <- lapply(v, function(x) x$"cor.train"[, "meanZ"])
t2 <- lapply(v, function(x) x$"cor.test"[, "meanZ"])
l <- mapply(function(x,y) {
  list(MBPCA = x, OLS = y)
}, t1, t2, SIMPLIFY = FALSE)
l <- unlist(l, recursive = FALSE)

png("Fig/benchmark_NCI60_train20_test40.png", width = 12, height = 12, units = "cm", res = 300)
par(mar = c(8, 4, 1, 1))
boxplot(l, las = 2, ylab = "Correlation (predicted vs measured)")
dev.off()

