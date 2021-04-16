context("F-stats")


set.seed(20161005)
## From limma::limFit example page:
sd <- 0.3 * sqrt(4 / rchisq(100, df = 4))
y <- matrix(rnorm(100 * 6, sd = sd), 100, 6)
rownames(y) <- paste("Gene", seq_len(100))
y[1:2, 4:6] <- y[1:2, 4:6] + 4

## Define the alternative and null models
pheno <- data.frame(group = rep(c(0, 1), each = 3), RIN = runif(6) + 8)
mod <- model.matrix(~ pheno$group + pheno$RIN)
mod0 <- model.matrix(~ pheno$RIN)

## Fit the models
library("limma")
fit <- lmFit(y, mod)
fit0 <- lmFit(y, mod0)

## Calculate the F statistics for these nested models
finfo <- getF(fit, fit0, y)

## Compute F-stats with derfinderHelper to double check
library("derfinderHelper")
fstat <- as.numeric(fstats.apply(data = y, mod = mod, mod0 = mod0, method = "regular"))
fpval <- pf(fstat, ncol(mod) - ncol(mod0), ncol(y) - ncol(mod), lower.tail = FALSE)

test_that("getF", {
    expect_equivalent(finfo$fstat, fstat)
    expect_equivalent(finfo$f_pval, fpval)
})

## Test getF with vector output from lm()

ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
weight <- c(ctl, trt)
lm.fit <- lm(weight ~ group)
lm.fit0 <- lm(weight ~ 1)

test_that("getF with vector", {
    expect_equal(getF(lm.fit, lm.fit0, weight)$f_pval, anova(lm.fit0, lm.fit)[["Pr(>F)"]][[2]])
    expect_equal(getF(lm.fit, lm.fit0, weight)$fstat, anova(lm.fit0, lm.fit)[["F"]][[2]])
})
