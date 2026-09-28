## ===========================================================================
## collect_hpc.R -- HPC ciktilarini kontrol et ve birlestir.
## ===========================================================================
library(SimDesign)

script_dir <- Sys.getenv("SIM_DIR", getwd())
source(file.path(script_dir, "functions.R"))

Design       <- make_design()
OUT_DIR      <- Sys.getenv("OUT_DIR", file.path(script_dir, "sim_results"))
FILENAME     <- "bifa"
ROWS_PER_JOB <- as.integer(Sys.getenv("ROWS_PER_JOB", "12"))

## 1) Eksik kosullar -> yeniden gonderilecek array ID'leri
files   <- file.path(OUT_DIR, paste0(FILENAME, "-", seq_len(nrow(Design)), ".rds"))
missing <- which(!file.exists(files))
if (length(missing)) {
  ids <- sort(unique(ceiling(missing / ROWS_PER_JOB)))
  message(length(missing), " kosul eksik. Yeniden gonderin:\n  sbatch --array=",
          paste(ids, collapse = ","), " submit_slurm.sh")
  message("(Uzun listede virgul siniri asilirsa parcalara bolun.)")
}
SimCheck(dir = OUT_DIR)

## 2) Birlestir (eksik varsa yalniz mevcutlar birlesir)
final <- SimCollect(dir = OUT_DIR)
final

## max_time/max_RAM nedeniyle erken kesilen kosullar: REPLICATIONS < 100
short <- final[final$REPLICATIONS < 100, c("ROW", "REPLICATIONS")]
if (nrow(short)) {
  message(nrow(short), " kosul eksik replikasyonla bitti; dosyalarini silip ",
          "ilgili array ID'lerini yeniden gonderin:")
  print(short)
}

saveRDS(final, file.path(script_dir, "bifa_sim_results.rds"))
saveRDS(to_long(final), file.path(script_dir, "bifa_sim_results_long.rds"))
