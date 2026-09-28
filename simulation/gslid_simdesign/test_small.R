library(SimDesign)
source("functions.R")

Design <- make_design()
FO     <- make_fixed_objects("population")

pick <- with(Design, which(
  (NGEN == 1 & NFAC == 3 & NVAR == 3 & GLOAD == "medium" & FLOAD == "medium" &
     FRHO == 0 & NOBS == 250 & NCAT %in% c(2, 6)) |
  (NGEN == 2 & NFAC == 4 & NVAR == 5 & GLOAD == "low" & FLOAD == "high" &
     GRHO == 0.5 & FRHO == 0.5 & NOBS == 1000 & LSKEW == "moderately" & NCAT == 3)))

res <- runSimulation(
  design        = Design[pick, ],
  replications  = 10,
  generate      = Generate,
  analyse       = Analyse,
  summarise     = Summarise,
  fixed_objects = FO,
  packages      = c("MASS", "bifactor", "fungible", "latentFactoR"),
  seed          = genSeeds(length(pick), iseed = 1276149341L),
  save          = FALSE,
  control       = list(allow_na = TRUE)
)

long <- to_long(res)
print(long[, c("ROW", "NGEN", "LSKEW", "NCAT", "Method", "conv_rate", "heywood_rate",
               "CC_general", "CC_specific", "AB_general", "AB_specific",
               "RB_general", "RB_specific")], digits = 3)
SimExtract(res, what = "errors")
