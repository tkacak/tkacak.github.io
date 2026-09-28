# Ordinal bifactor EFA simulation (SimDesign)

| File | Purpose |
|---|---|
| `functions.R` | Design, `Generate`, `Analyse`, `Summarise`, alignment/metric helpers (single source) |
| `test_small.R` | 6-condition × 10-rep smoke test; run first |
| `run_local.R` | Full design on one machine (`runSimulation`, resumable) |
| `run_hpc.R` + `submit_slurm.sh` | SLURM job array (`runArraySimulation`, 12 conditions per task → 729 tasks) |
| `collect_hpc.R` | `SimCheck` + `SimCollect`, lists array IDs to resubmit |

Metrics per method (ULS, ML, PA, REGULS, REGML, PCA), separately for general and specific factors:

- **CC**: Tucker's congruence between aligned estimated and true loading columns (mean over replications)
- **AB**: `mean_p | mean_r(λ̂_pr) − λ_p |` over salient loadings
- **RB (%)**: `mean_p (mean_r(λ̂_pr) − λ_p) / λ_p × 100` over salient loadings
- **MAE**: `mean_p mean_r | λ̂_pr − λ_p |` (what the old `metrics.R` called "AB")
- **conv_rate**: share of replications where the method returned a usable solution

Output is wide (`ULS.CC_general`, …); `to_long()` gives one row per condition × method.
