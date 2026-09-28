## ===========================================================================
## test_small.R -- tam calistirmadan once birkac kosulda hizli kontrol.
## conv_rate ~0 ise BiFAD cikti yapisi beklenenden farklidir (loading_slot).
## ===========================================================================
library(SimDesign)
source("functions.R")

Design <- make_design()
FO     <- make_fixed_objects("population")

## Uc uca kosullar: en kucuk/en buyuk model, her LSKEW, NCAT = 2 ve 7
pick <- with(Design, which(
  (NFAC == 2 & NVAR == 4 & GLOAD == "medium" & FLOAD == "medium" & FRHO == 0 &
     NOBS == 250 & NCAT %in% c(2, 7)) |
  (NFAC == 4 & NVAR == 8 & GLOAD == "low" & FLOAD == "high" & FRHO == 0.5 &
     NOBS == 1000 & LSKEW == "moderately" & NCAT == 3)))

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
print(long[, c("ROW", "LSKEW", "NCAT", "Method", "conv_rate", "CC_general",
               "CC_specific", "AB_general", "AB_specific", "RB_general",
               "RB_specific")], digits = 3)
SimExtract(res, what = "errors")
