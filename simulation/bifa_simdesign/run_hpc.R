library(SimDesign)
source("functions.R")

Design <- make_design()
FO     <- make_fixed_objects("population")
n_rep  <- 100
iseed  <- 1276149341L
hpc    <- hpc_settings()
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

array2row <- function(arrayID) {
  rows <- ((arrayID - 1L) * hpc$rows_per_job + 1L):(arrayID * hpc$rows_per_job)
  rows <- rows[rows <= nrow(Design)]
  rows[!file.exists(file.path(hpc$out_dir, sprintf("%s-%d.rds", hpc$filename, rows)))]
}

arrayID <- getArrayID(type = "slurm")
todo    <- array2row(arrayID)
if (!length(todo)) {
  message("arrayID ", arrayID, ": no rows to run (finished or outside the design; ",
          "required tasks: ", ceiling(nrow(Design) / hpc$rows_per_job), "). Exiting.")
  quit(save = "no", status = 0)
}
message("arrayID ", arrayID, " -> rows: ", paste(todo, collapse = ", "))

control <- list(allow_na = TRUE)
if (nzchar(Sys.getenv("MAX_TIME"))) control$max_time <- Sys.getenv("MAX_TIME")
if (nzchar(Sys.getenv("MAX_RAM")))  control$max_RAM  <- Sys.getenv("MAX_RAM")

runArraySimulation(
  design        = Design,
  replications  = n_rep,
  generate      = Generate,
  analyse       = Analyse,
  summarise     = Summarise,
  fixed_objects = FO,
  packages      = c("MASS", "bifactor", "fungible", "latentFactoR"),
  iseed         = iseed,
  arrayID       = arrayID,
  array2row     = array2row,
  filename      = hpc$filename,
  dirname       = hpc$out_dir,
  parallel      = ncores > 1L,
  ncores        = ncores,
  control       = control
)
