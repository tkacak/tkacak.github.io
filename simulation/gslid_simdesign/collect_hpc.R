library(SimDesign)
source("functions.R")

Design <- make_design()
hpc    <- hpc_settings()
n_rep  <- 100

files   <- file.path(hpc$out_dir, sprintf("%s-%d.rds", hpc$filename, seq_len(nrow(Design))))
missing <- which(!file.exists(files))
if (length(missing)) {
  ids <- sort(unique(ceiling(missing / hpc$rows_per_job)))
  message(length(missing), " conditions missing. Resubmit with:\n  sbatch --array=",
          paste(ids, collapse = ","), " submit_slurm.sh")
}
SimCheck(dir = hpc$out_dir)

final <- SimCollect(dir = hpc$out_dir)
final

short <- final[final$REPLICATIONS < n_rep, c("ROW", "REPLICATIONS")]
if (nrow(short)) {
  message(nrow(short), " conditions finished with fewer than ", n_rep,
          " replications; delete their files and resubmit their array IDs:")
  print(short)
}

saveRDS(final, "gslid_sim_results.rds")
saveRDS(to_long(final), "gslid_sim_results_long.rds")
