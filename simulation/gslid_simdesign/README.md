# Ordinal exploratory bifactor simulation with GSLiD (SimDesign)

Monte Carlo simulation of exploratory bifactor analysis with one or two general factors, estimated with the
GSLiD algorithm of `bifactor::bifactor()` (ULS and ML) on polychoric correlations of ordinal data. The study
is a SimDesign `Generate` → `Analyse` → `Summarise` pipeline that runs on one machine (`runSimulation`) or
as a SLURM job array (`runArraySimulation`); `analysis.R` produces the tables, ANOVAs and figures.

## Files

| File | Purpose |
|---|---|
| `functions.R` | Design, fixed objects, `Generate`, `Analyse`, `Summarise`, alignment and metric helpers, `to_long()`, HPC settings. Every other simulation script `source()`s it. |
| `test_small.R` | Smoke test: 7 conditions × 10 replications. |
| `run_local.R` | Full design on one machine. |
| `submit_slurm.sh` | SLURM job-array submission script. |
| `run_hpc.R` | What each array task runs. |
| `collect_hpc.R` | Checks and combines the array output. |
| `analysis.R` | Heywood and metric tables, factorial ANOVAs with partial ω², figures. |

All scripts expect this folder to be the working directory.

## Requirements

- Simulation: SimDesign, bifactor, fungible, latentFactoR, MASS.
  bifactor is installed from GitHub: `devtools::install_github("marcosjnez/bifactor", force = TRUE)`
  (compiled C++; the simulation was checked with version 0.1.1).
- `analysis.R`: dplyr, tidyr, purrr, stringr, broom, effectsize, knitr, kableExtra, ggplot2, writexl.

The current CRAN release of fungible (2.4.8) requires R ≥ 4.5.0, and its dependency CVXR (1.9.2) requires
Matrix ≥ 1.7. On an older R (e.g. a cluster module), install the versions you tested with locally
(`packageVersion()`, `remotes::install_version()`).

## How to run

### 1. Smoke test

```r
source("test_small.R")
```

Runs 7 conditions (one-general model with NCAT = 2 and 6 at each LSKEW level, plus the largest two-general
model with GRHO = FRHO = .5) with 10 replications and prints one row per condition × method.

### 2. Full run on one machine

```r
source("run_local.R")
```

- 19440 conditions × 100 replications; replications within a condition run in parallel on all cores but two.
- If R stops, run the script again: SimDesign resumes from the last finished condition.
- Output: `gslid_sim_results.rds` (wide) and `gslid_sim_results_long.rds` (long, used by `analysis.R`).
- `store_results = FALSE` keeps RAM low, so replication-level results are not kept. Add
  `save_results = TRUE` if summaries must be recomputable later.

### 3. Full run on a SLURM cluster

1. Copy this folder to the cluster and install the packages.
2. In `submit_slurm.sh`, set `module load r` to your cluster's R module and adjust `--cpus-per-task`,
   `--mem` and `--time` if needed.
3. `mkdir -p logs` (SLURM does not create it; jobs fail without it).
4. From this folder: `sbatch submit_slurm.sh`
5. When all tasks have finished: `Rscript collect_hpc.R`, then `Rscript analysis.R`.

- Each task runs `ROWS_PER_JOB` = 20 design rows, so the array is `1-972` (= 19440 / 20), below the common
  `MaxArraySize` of 1001. With a different `ROWS_PER_JOB` (environment variable), change the `--array`
  upper bound to ⌈19440 / ROWS_PER_JOB⌉ and use the same value when running `collect_hpc.R`.
- Each condition is saved as `gslid-<ROW>.rds` in `$SCRATCH/gslid_sim_results`, or in
  `./gslid_sim_results` when `$SCRATCH` is not set; `OUT_DIR` overrides this for both the jobs and
  `collect_hpc.R`.
- A task skips rows whose file already exists, and array IDs beyond the design exit immediately, so
  resubmitting never reruns a finished condition.
- `MAX_TIME` and `MAX_RAM` (about 90% of `--time` and `--mem`) make SimDesign stop a condition before the
  scheduler kills the job, keeping the replications finished so far.
- `collect_hpc.R` prints an `sbatch --array=…` command for missing conditions, runs `SimCheck()`, combines
  the files with `SimCollect()`, and lists conditions with fewer than 100 replications (delete those files
  and resubmit their array IDs; array ID = ⌈ROW / ROWS_PER_JOB⌉).
- Every HPC file also stores the replication-level results (`SimResults()`).

## Design (`make_design()`)

| Factor | Levels | Meaning |
|---|---|---|
| `NGEN` | 1, 2 | general factors |
| `NFAC` | 3, 4 | group factors per general factor |
| `NVAR` | 3, 5 | items per group factor |
| `GLOAD` | low, medium, high | general loadings from U(.3, .5), U(.4, .6), U(.5, .7) |
| `FLOAD` | low, medium, high | group loadings, same ranges |
| `GRHO` | 0, .2, .5 | correlation among general factors |
| `FRHO` | 0, .2, .5 | correlation among group factors |
| `NOBS` | 250, 500, 1000 | sample size |
| `LSKEW` | normal, slightly, moderately | latent distribution (see Data generation) |
| `NCAT` | 2, 3, 4, 5, 6 | response categories |

Conditions with one general factor and GRHO > 0 are removed (a single general factor has no correlation), which
leaves 19440 conditions. Two identifier columns are added:

- `ROW` — row number after that filter. It equals the index `i` of `current_condition` in the original
  scripts; NCAT varies slowest, so leaving out NCAT = 7 removed only the last rows.
- `POP_ID` (1–432) — the population model, shared by conditions with the same NGEN, NFAC, NVAR, GLOAD,
  FLOAD, GRHO and FRHO.

NCAT = 7 is not part of the design because `latentFactoR::categorize()` handles more than six categories
differently: it cuts each sample's observed range into equal-width bins instead of using fixed thresholds.

## Data generation (`Generate()`)

1. `population_model()`: `bifactor::sim_factor(..., crossloadings = 0, method = "minres")` returns the true
   loadings `lambda` (items × (NGEN + NGEN·NFAC); general factors in the first NGEN columns, each loading on
   its own half of the items when NGEN = 2) and the population correlation matrix `R`.
2. Latent data (NOBS × items) with correlation `R`: `MASS::mvrnorm` for `normal`; `fungible::monte1` with
   skewness 0.5 / kurtosis 1.5 (`slightly`) or 1.5 / 3 (`moderately`) for every item.
3. `latentFactoR::categorize(categories = NCAT, skew_value = 0)` applied column by column. With
   `skew_value = 0` the fixed thresholds give categories with zero skewness, which for 4–6 categories is
   not the same as symmetric (for normal latent data: 4 categories 12.6 / 35.8 / 28.6 / 22.9%,
   5 categories 12.4 / 24.7 / 22.7 / 20.9 / 19.2%, 6 categories 11.0 / 19.3 / 18.5 / 17.8 / 17.1 / 16.4%).

## Analysis (`Analyse()`)

1. Polychoric correlations once per replication (`bifactor::polyfast()`), shared by both estimators. If they
   fail, `Analyse()` stops and SimDesign draws a new data set (counted in `ERRORS`).
2. For each estimator (`uls`, `ml`): `bifactor::bifactor()` with `method = "GSLiD"`,
   `n_generals = NGEN`, `n_groups = NGEN·NFAC`, `projection = "oblq"`, `oblq_factors = NGEN·NFAC`,
   `maxit = 100`, `random_starts = 10`, `rot_control = list(maxit = 1e4, rotation = "oblimin")`.
   Settings other than the estimator are in `make_fixed_objects()$gslid`.
3. Factor-correlation target (`phi_target()`): the general factors may correlate freely (weight 0); every
   other correlation has target 0 and weight 1.
4. A fit that errors, or whose loading matrix has the wrong size or non-finite values, is recorded as `NA`
   for that estimator, and the data set is kept (redrawing would hide the failure and inflate its
   convergence rate). `control = list(allow_na = TRUE)` lets SimDesign accept these `NA`s. Some ML fits in
   small models stop inside bifactor with `Col::subvec(): indices out of bounds`; they count as
   non-converged.
5. The Heywood flag is `fit$efa$heywood` (any uniqueness ≤ 0 in the first-order EFA).
6. The estimate is aligned to the true loadings, and `Analyse()` returns `CC_general`, `CC_specific`,
   `heywood` and the aligned salient loadings (`g1, …`, `s1, …`), prefixed with the method name.

### Alignment (`align_loadings()`)

- General factors, one at a time: the estimated column with the highest absolute Tucker congruence with the
  true general column, reflected to point in the same direction.
- Group factors: the remaining columns are matched and reflected with
  `fungible::faAlign(..., MatchMethod = "CC")` (Hungarian algorithm when the best matches are not unique).

## Metrics (`Summarise()`)

For each method, separately for general and group loadings, over the replications in which the method
converged. With λₚ a salient true loading, λ̂ₚᵣ its aligned estimate in replication r, P the salient
loadings in the block and R the converged replications:

| Column | Definition |
|---|---|
| `conv_rate` | proportion of replications with a usable solution |
| `heywood_rate` | proportion of converged replications with a Heywood case |
| `CC_general`, `CC_specific` | Tucker's congruence φₖ = Σᵢ λᵢₖ λ̂ᵢₖ / √(Σᵢ λᵢₖ² · Σᵢ λ̂ᵢₖ²) per factor, averaged over the block's factors and then over replications |
| `AB_general`, `AB_specific` | absolute bias: (1/P) Σₚ \| mean_r(λ̂ₚᵣ) − λₚ \| |
| `RB_general`, `RB_specific` | relative bias as a proportion (−0.04 = −4%): (1/P) Σₚ (mean_r(λ̂ₚᵣ) − λₚ) / λₚ |
| `MAE_general`, `MAE_specific` | mean absolute error: (1/(P·R)) Σₚ Σᵣ \| λ̂ₚᵣ − λₚ \| |

- AB, RB and MAE use only salient true loadings (|λ| > 1e-8). With two general factors the general columns
  contain zeros (each general factor loads on half of the items); these are excluded from AB/RB/MAE but
  included in CC.
- AB is the bias of the mean estimate; MAE also contains sampling variability. MAE corresponds to the
  "ABS"/"MAB" measure of the original scripts.
- AB and RB are computed with `SimDesign::bias()`.

## Output

- Wide: one row per condition with the design columns, SimDesign's columns (`REPLICATIONS`, `SIM_TIME`,
  `ERRORS`, `WARNINGS`, …) and 2 × 10 metric columns named `<METHOD>.<metric>` (e.g. `ULS.CC_general`).
- Long (`to_long()`): one row per condition × method with a `Method` column.

## Reproducibility

- Replications: `iseed = 1276149341` (`genSeeds(Design, iseed = …)` locally, the `iseed` argument of
  `runArraySimulation()` on the cluster) gives each condition its own L'Ecuyer-CMRG stream. Keep it
  unchanged between runs and resubmissions.
- Population loadings: `population_model()` seeds `sim_factor()` with `set.seed(seed_base + POP_ID)`
  (Mersenne-Twister) inside `with_preserved_rng()`, which restores the simulation's random-number state.
  `Generate()` and `Summarise()` therefore rebuild the same matrix, and conditions sharing a population use
  the same loadings. `make_fixed_objects("row")` seeds with `seed_base + ROW` instead (1234 + i, as in the
  original data generation script).
- `fungible::monte1()` calls `set.seed()` internally; `monte1_latent()` draws that seed from the replication
  stream and restores the stream afterwards.

## `analysis.R`

Reads `gslid_sim_results_long.rds` and reshapes it to one row per condition × method × metric × order
(General/Specific). Settings at the top: `factors`, `metrics` (CC, AB, RB), `heywood_cut` (.05) and
`omega_cut` (.14).

- "Clean" conditions: condition × method combinations with `heywood_rate < heywood_cut`.
- Heywood table: mean (SD) of `heywood_rate` per factor level and method.
- Metric table: mean per factor level, method, metric and order; clean-condition means in parentheses.
- ANOVAs: for each method × metric × order, `aov(Value ~ (all ten factors)^3)` on the clean conditions,
  with partial ω² and 95% CIs from `effectsize::omega_squared()`. NGEN and GRHO are not fully crossed
  (GRHO > 0 only with two general factors), so the NGEN:GRHO interaction is not estimable and does not
  appear in the tables. The full tables go to the supplement; effects with partial ω² ≥ `omega_cut` form
  the main-text table.
- `gslid_results_noheywood.xlsx`: the clean long data.
- Figures: group-factor CC by NCAT, faceted by FLOAD (reference line at .95), and group-factor AB by NCAT,
  faceted by NFAC. `coord_cartesian()` limits the y-axis without dropping data.

## Functions in `functions.R`

| Function | Role |
|---|---|
| `make_design()` | design with `ROW` and `POP_ID` |
| `make_fixed_objects()` | estimators, GSLiD settings, salience threshold, seed settings |
| `hpc_settings()` | `ROWS_PER_JOB`, `OUT_DIR` and file prefix for the HPC scripts |
| `with_preserved_rng()` | evaluates an expression, then restores `.Random.seed` |
| `population_model()` | seeded `sim_factor()` → true loadings and population correlation matrix |
| `monte1_latent()` | non-normal latent data through `fungible::monte1()` |
| `phi_target()` | target and weight matrices for the factor correlations |
| `factor_index()` | column indices of general and group factors; checks the matrix size |
| `salient_mask()` | which true loadings are salient |
| `tucker_cc()`, `cc_against()` | Tucker's congruence |
| `extract_lambda()`, `extract_heywood()` | loadings and Heywood flag from a `bifactor()` fit |
| `align_loadings()` | factor matching and reflection |
| `rep_output()` | one method's output for one replication |
| `Generate()`, `Analyse()`, `Summarise()` | SimDesign steps |
| `summarise_block()` | AB, RB and MAE for one block of loadings |
| `to_long()` | wide results → one row per condition × method |
