## ===========================================================================
## run_local.R -- tek makinede SimDesign::runSimulation ile tum tasarim.
## Eski data_generate.R + maincode.R + loop.R + metrics.R'nin yerini alir:
## veri diske yazilmaz, her replikasyon uretilip hemen analiz edilir.
##
## Kesilirse (elektrik, cokme) ayni scripti tekrar calistirin: SimDesign
## gecici dosyadan kaldigi kosuldan devam eder (resume = TRUE).
## ===========================================================================
library(SimDesign)

## Bu dosyanin bulundugu klasor; gerekirse elle degistirin.
sim_dir <- "D:/Ranalysis/makale/EPOD_Ordered_EBIFA/publishng/JMEEP/RV1/RV_AE/simdesign"
setwd(sim_dir)
source("functions.R")

Design <- make_design()                 # 8748 kosul
FO     <- make_fixed_objects("population")
n_rep  <- 100

## iseed bir kez genSeeds() ile uretilip SABIT tutulur (HPC ile ayni deger).
iseed <- 1276149341L
seeds <- genSeeds(Design, iseed = iseed)

res <- runSimulation(
  design        = Design,
  replications  = n_rep,
  generate      = Generate,
  analyse       = Analyse,
  summarise     = Summarise,
  fixed_objects = FO,
  packages      = c("MASS", "bifactor", "fungible", "latentFactoR"),
  seed          = seeds,
  parallel      = TRUE,
  ncores        = max(1L, parallelly::availableCores(omit = 2L)),
  filename      = "bifa_sim_results",   # bitince bifa_sim_results.rds
  store_results = FALSE,                # 8748 x 100 satir RAM'e sigmaz
  ## Replikasyon duzeyinde hizalanmis yukleri saklamak isterseniz
  ## (ek analiz / reSummarise icin) asagidakini acin; ~2-3 GB disk.
  # save_results = TRUE,
  control       = list(allow_na = TRUE)  # basarisiz yontem = NA, veri yeniden cekilmez
)

res
long <- to_long(res)
saveRDS(long, "bifa_sim_results_long.rds")

## Hata / uyari kontrolu
SimExtract(res, what = "errors")
SimExtract(res, what = "warnings")
