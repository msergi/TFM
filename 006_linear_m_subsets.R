# options(echo = TRUE)
setwd('/home/isglobal.lan/smari/data/smari/data_res2')
library("parallel")

arg <- commandArgs(trailingOnly = T)
a <- paste('easub', arg, '.RData', sep = '') # object: easub
b <- paste('masub', arg, '.RData', sep = '') # object: masub
c <- paste('overlaps', arg, '.RData', sep = '') # object: overlaps
load(a)
load(b)
load(c)
load("pheno.RData")

stopifnot(colnames(masub) == colnames(easub))
stopifnot(rownames(pheno) == colnames(masub))


resCorr <- function(x) {
  cpg <- overlaps[x, 1]
  tc <- overlaps[x, 2]
  fit <- lm(easub[tc,] ~ masub[cpg,] + pheno$cohort + pheno$e3_sex +
              pheno$age_sample_years  + pheno$NK_6 + pheno$Bcell_6 +
              pheno$CD4T_6 + pheno$CD8T_6 + pheno$Gran_6 + pheno$Mono_6)
  return(c(cpg,
           tc,
           summary(fit)$coef[2,1],
           summary(fit)$coef[2,2],
           summary(fit)$coef[2,4],
           confint(fit)[2,1],
           confint(fit)[2,2]))
}

len <- 1:nrow(overlaps)
results <- mclapply(len, resCorr, mc.cores = 16)
output <- matrix(unlist(results), ncol = 7, byrow = TRUE)
filename <- paste('output', arg, '.txt', sep = '')
write.table(output, file = filename, col.names = F, row.names = F,
            quote = F)