# k3.inp karşılığı — tidySEM
library(tidySEM)
library(OpenMx)
source("fmm_k3_generate.R")          # sim, items, nFac, nPer, J

# Sınıf k için lavaan sözdizimi; sınıf 1 faktör ortalamaları 0 ([f@0]), diğerleri serbest
lav <- function(k) paste(c(
  sapply(1:nFac, function(f) {
    it <- (f - 1) * nPer + 1:nPer
    sprintf("f%d =~ NA*y%d + %s", f, it[1], paste0("l", it, "*y", it, collapse = " + "))
  }),
  sprintf("f%d ~~ 1*f%d", 1:nFac, 1:nFac),
  "f1 ~~ c12*f2", "f1 ~~ c13*f3", "f2 ~~ c23*f3",
  sprintf("y%d ~~ e%d*y%d", 1:J, 1:J, 1:J),
  sprintf("y%d ~ i%d*1", 1:J, 1:J),
  if (k == 1) sprintf("f%d ~ 0*1", 1:nFac) else sprintf("f%d ~ m%d_f%d*1", 1:nFac, k, 1:nFac)
), collapse = "\n")

set.seed(1)
fit2 <- mx_mixture(model = lapply(1:2, lav), data = sim$data)   # TECH11 için K-1
fit3 <- mx_mixture(model = lapply(1:3, lav), data = sim$data)
stopifnot(length(omxGetParameters(fit3)) == 99)                  # Mplus Npar

table_fit(fit3)
class_prob(fit3, type = c("sum.posterior", "avg.mostlikely"))
table(true = sim$class, assigned = class_prob(fit3, type = "individual")$individual[, "predicted"])
lr_lmr(fit2, fit3)                     # TECH11'in VLMR'si DEĞİL (nonnest2 Vuong testi)
BLRT(fit2, fit3, replications = 100)   # TECH14 benzeri
