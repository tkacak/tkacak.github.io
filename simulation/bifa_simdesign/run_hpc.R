## ===========================================================================
## run_hpc.R -- SLURM job array ile SimDesign::runArraySimulation.
## Her array gorevi ROWS_PER_JOB kosul calistirir; her kosul ayri bir
## dosyaya yazilir: <OUT_DIR>/bifa-<satir>.rds
## Gonderim: submit_slurm.sh  |  Birlestirme: collect_hpc.R
## ===========================================================================
library(SimDesign)

script_dir <- Sys.getenv("SIM_DIR", getwd())
source(file.path(script_dir, "functions.R"))

Design <- make_design()                 # 7290 kosul
FO     <- make_fixed_objects("population")
n_rep  <- 100

## run_local.R ile ayni iseed: iki surum ayni RNG akislarini kullanir.
iseed <- 1276149341L

ROWS_PER_JOB <- as.integer(Sys.getenv("ROWS_PER_JOB", "12"))  # ceiling(7290/12) = 608 gorev
OUT_DIR      <- Sys.getenv("OUT_DIR", file.path(script_dir, "sim_results"))
FILENAME     <- "bifa"
NCORES       <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
MAX_TIME     <- Sys.getenv("MAX_TIME", "")       # orn. "21:30:00" (--time'in ~%90'i)
MAX_RAM      <- Sys.getenv("MAX_RAM", "")        # orn. "7GB"  (--mem'in ~%90'i)

## arrayID -> tasarim satirlari. Tasarim disina tasan satirlar ve zaten
## tamamlanmis kosullar atlanir; boylece --array fazla verilse de ya da eksik
## kosullar icin ayni ID'ler yeniden gonderilse de ayni satir iki kez kosmaz.
array2row <- function(arrayID) {
  rows <- ((arrayID - 1L) * ROWS_PER_JOB + 1L):(arrayID * ROWS_PER_JOB)
  rows <- rows[rows <= nrow(Design)]
  ## sprintf bos girdide bos doner (paste0 "bifa-.rds" uretir -> NA satir)
  rows[!file.exists(file.path(OUT_DIR, sprintf("%s-%d.rds", FILENAME, rows)))]
}

arrayID <- getArrayID(type = "slurm")
todo    <- array2row(arrayID)
if (!length(todo)) {
  message("arrayID ", arrayID, ": calistirilacak kosul yok (tamamlanmis ya da ",
          "tasarim disi; gereken gorev sayisi ",
          ceiling(nrow(Design) / ROWS_PER_JOB), "), cikiliyor.")
  quit(save = "no", status = 0)
}
message("arrayID ", arrayID, " -> satirlar: ", paste(todo, collapse = ", "))

control <- list(allow_na = TRUE)
if (nzchar(MAX_TIME)) control$max_time <- MAX_TIME
if (nzchar(MAX_RAM))  control$max_RAM  <- MAX_RAM

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
  filename      = FILENAME,
  dirname       = OUT_DIR,
  parallel      = NCORES > 1L,
  ncores        = NCORES,
  control       = control
)
