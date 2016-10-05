context('F-stats')


set.seed(20161005)
## From limma::limFit example page:
sd <- 0.3 * sqrt(4 / rchisq(100, df = 4))
y <- matrix(rnorm(100 * 6, sd = sd), 100, 6)
rownames(y) <- paste('Gene', seq_len(100))
y[1:2, 4:6] <- y[1:2, 4:6] + 4

## Define the alternative and null models
pheno <- data.frame(group = rep(c(0, 1), each = 3), RIN = runif(6) + 8)
mod <- model.matrix(~ pheno$group + pheno$RIN)
mod0 <- model.matrix(~ pheno$RIN)

## Fit the models
library('limma')
fit <- lmFit(y, mod)
fit0 <- lmFit(y, mod0)

## Calculate the F statistics for these nested models
finfo <- getF(fit, fit0, y)

## Compute F-stats with derfinderHelper to double check
library('derfinderHelper')
fstat <- as.numeric(fstats.apply(data = y, mod = mod, mod0 = mod0, method = 'regular'))
fpval <- pf(fstat, ncol(mod) - ncol(mod0), ncol(y) - ncol(mod), lower.tail = FALSE)

test_that('getF', {
    expect_equivalent(finfo$fstat, fstat)
    expect_equivalent(finfo$f_pval, fpval)
})
