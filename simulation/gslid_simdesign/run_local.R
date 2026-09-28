library(SimDesign)
source("functions.R")

Design <- make_design()
FO     <- make_fixed_objects("population")
n_rep  <- 100
iseed  <- 1276149341L

res <- runSimulation(
  design        = Design,
  replications  = n_rep,
  generate      = Generate,
  analyse       = Analyse,
  summarise     = Summarise,
  fixed_objects = FO,
  packages      = c("MASS", "bifactor", "fungible", "latentFactoR"),
  seed          = genSeeds(Design, iseed = iseed),
  parallel      = TRUE,
  ncores        = max(1L, parallelly::availableCores(omit = 2L)),
  filename      = "gslid_sim_results",
  store_results = FALSE,
  control       = list(allow_na = TRUE)
)

res
saveRDS(to_long(res), "gslid_sim_results_long.rds")
SimExtract(res, what = "errors")
SimExtract(res, what = "warnings")
