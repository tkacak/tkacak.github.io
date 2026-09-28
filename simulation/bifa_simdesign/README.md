# Ordinal bifactor EFA simulation (SimDesign)

Monte Carlo simulation comparing six factor-extraction methods for direct Schmid–Leiman bifactor analysis
(`fungible::BiFAD`) on polychoric correlations of ordinal data. The study is written as a SimDesign
`Generate` → `Analyse` → `Summarise` pipeline and runs either on one machine (`runSimulation`) or as a
SLURM job array (`runArraySimulation`).

## Files

| File | Purpose |
|---|---|
| `functions.R` | Everything the simulation needs: design, fixed objects, `Generate`, `Analyse`, `Summarise`, alignment and metric helpers, `to_long()`, HPC settings. Every other script `source()`s it. |
| `test_small.R` | Smoke test: 7 conditions × 10 replications. |
| `run_local.R` | Full design on one machine. |
| `submit_slurm.sh` | SLURM job-array submission script. |
| `run_hpc.R` | What each array task runs. |
| `collect_hpc.R` | Checks and combines the array output. |

All scripts expect this folder to be the working directory.

## Requirements

SimDesign, bifactor (<https://github.com/Marcosjnez/bifactor>, compiled C++), fungible, latentFactoR and MASS.
On a cluster, install them into your user library before submitting.

The current CRAN release of fungible (2.4.8) requires R ≥ 4.5.0, and its dependency CVXR (1.9.2) requires
Matrix ≥ 1.7. On an older R, `install.packages()` will not install them; install the versions you tested
with locally instead (check with `packageVersion()`, install with `remotes::install_version()`).
SimDesign 2.27 needs qs2 → stringfish → RcppParallel ≥ 6.1.1.

## How to run

### 1. Smoke test

```r
source("test_small.R")
```

Runs 7 conditions (NCAT = 2 and 6 at each LSKEW level, plus the largest model with FRHO = .5) with
10 replications and prints one row per condition × method. `conv_rate` should be near 1 for every method.
A method with `conv_rate = 0` means BiFAD's output does not contain the loading matrix that
`extract_lambda()` looks for (`loading_slot`).

### 2. Full run on one machine

```r
source("run_local.R")
```

- Runs all 7290 conditions × 100 replications; within each condition the replications run in parallel on
  all cores but two.
- SimDesign keeps a temporary file after every condition. If R stops, run the script again and it
  resumes from the last finished condition.
- Output: `bifa_sim_results.rds` (SimDesign object, wide) and `bifa_sim_results_long.rds` (long).
- `store_results = FALSE`: replication-level results are not kept, because 7290 × 100 rows do not fit
  in RAM. Summaries therefore cannot be recomputed later without rerunning. Add `save_results = TRUE` to
  write them to disk (roughly 1–2 GB) if you need that.

### 3. Full run on a SLURM cluster

1. Copy this folder to the cluster and install the packages.
2. In `submit_slurm.sh`, adjust `module load r` to your cluster's R module and, if needed,
   `--cpus-per-task`, `--mem` and `--time`.
3. `mkdir -p logs` — SLURM does not create the log folder, and jobs fail without it.
4. From this folder: `sbatch submit_slurm.sh`
5. When all tasks have finished: `Rscript collect_hpc.R`

How the array works:

- Each task runs `ROWS_PER_JOB` = 12 design rows, so the array is `1-608` (= ⌈7290 / 12⌉). If your
  cluster's `MaxArraySize` is smaller, raise `ROWS_PER_JOB` (environment variable read by
  `hpc_settings()`), lower the `--array` upper bound to match, and use the same value when running
  `collect_hpc.R`.
- Each condition is saved as `bifa-<ROW>.rds` in `$SCRATCH/bifa_sim_results`, or in
  `./bifa_sim_results` when `$SCRATCH` is not set. The `OUT_DIR` environment variable overrides this;
  set it for both the jobs and `collect_hpc.R`.
- A task skips rows whose file already exists, and array IDs beyond the design exit immediately, so
  resubmitting an array ID never reruns a finished condition.
- Within a task, replications run in parallel on `SLURM_CPUS_PER_TASK` cores.
- `MAX_TIME` and `MAX_RAM` (set in `submit_slurm.sh` to about 90% of `--time` and `--mem`) make SimDesign
  stop a condition before the scheduler kills the job, keeping the replications finished so far.
- `collect_hpc.R` prints an `sbatch --array=…` command for missing conditions (split it if it is too long
  for your scheduler), runs `SimCheck()`, combines the files with `SimCollect()`, and lists conditions that
  finished with fewer than 100 replications. Delete those files and resubmit their array IDs
  (array ID = ⌈ROW / ROWS_PER_JOB⌉). Output: `bifa_sim_results.rds` and `bifa_sim_results_long.rds`.
- Unlike the local run, every HPC file also stores the replication-level results (`SimResults()`), so new
  summaries can be computed later without rerunning the simulation.

## Design (`make_design()`)

| Factor | Levels | Meaning |
|---|---|---|
| `NGEN` | 1 | general factors |
| `NFAC` | 2, 3, 4 | group factors per general factor |
| `NVAR` | 4, 8 | items per group factor |
| `GLOAD` | low, medium, high | general-factor loadings drawn from U(.3, .5), U(.4, .6), U(.5, .7) |
| `FLOAD` | low, medium, high | group-factor loadings, same ranges |
| `FRHO` | 0, .3, .5 | correlation among group factors |
| `GRHO` | 0 | correlation among general factors |
| `NOBS` | 250, 500, 1000 | sample size |
| `CROSS` | 0 | cross-loadings |
| `LSKEW` | normal, slightly, moderately | latent distribution (see Data generation) |
| `NCAT` | 2, 3, 4, 5, 6 | response categories |

Fully crossed: 3 × 2 × 3 × 3 × 3 × 3 × 3 × 5 = 7290 conditions. Two identifier columns are added:

- `ROW` — row number. The design is built with `expand.grid` in the same factor order as the original
  scripts, so `ROW` equals the original condition number. NCAT varies slowest, so leaving out NCAT = 7
  removed only the last 1458 rows (7291–8748).
- `POP_ID` (1–162) — the population model. Conditions with the same NGEN, NFAC, NVAR, GLOAD, FLOAD, FRHO,
  GRHO and CROSS share it; NOBS, LSKEW and NCAT do not change the population.

NCAT = 7 is not part of the design because `latentFactoR::categorize()` treats more than six categories
differently: instead of fixed thresholds it cuts each sample's observed range into equal-width bins, so
category proportions would depend on N and on the sample extremes.

## Fixed objects (`make_fixed_objects()`)

- `methods` — output prefix → `facMethod` of `fungible::BiFAD`: ULS = `"fals"`, ML = `"faml"`,
  PA = `"fapa"`, REGULS = `"faregLS"`, REGML = `"faregML"`, PCA = `"pca"`.
- `loading_slot = "BstarSL"` — the BiFAD element used as the estimate: direct Schmid–Leiman loadings,
  items × (1 + NFAC), general factor in the first column.
- `min_true = 1e-8` — true loadings with |λ| above this are salient.
- `seed_base = 1234` and `pop_seed` — see Reproducibility.

## Data generation (`Generate()`)

1. `population_model()` calls `bifactor::sim_factor(..., method = "minres")` and returns the true loading
   matrix `lambda` (items × (1 + NFAC), general factor first) and the population correlation matrix `R`.
2. Latent continuous data (NOBS × items) with correlation `R`:
   - `normal`: `MASS::mvrnorm`
   - `slightly`: `fungible::monte1` with skewness 0.5 and kurtosis 1.5 for every item (`monte1_latent()`)
   - `moderately`: `fungible::monte1` with skewness 1.5 and kurtosis 3
3. Categorization: `latentFactoR::categorize(categories = NCAT, skew_value = 0)` applied column by column.
   The function is written for a single vector, and latentFactoR's own `simulate_factors()` calls it the
   same way. With `skew_value = 0`, the thresholds come from latentFactoR's table and are fixed; they give
   categories with zero skewness, which for 4–6 categories is not the same as symmetric. Proportions for
   normal latent data:

   | NCAT | % per category |
   |---|---|
   | 2 | 50 / 50 |
   | 3 | 30 / 40 / 30 |
   | 4 | 12.6 / 35.8 / 28.6 / 22.9 |
   | 5 | 12.4 / 24.7 / 22.7 / 20.9 / 19.2 |
   | 6 | 11.0 / 19.3 / 18.5 / 17.8 / 17.1 / 16.4 |

   With skewed latent data the same fixed thresholds produce skewed categories.

`Generate()` returns `list(dat = <ordinal data>, lambda = <true loadings>)`.

## Analysis (`Analyse()`)

1. Polychoric correlations are computed once per replication with `bifactor::polyfast()` and shared by all
   six methods. If they fail or contain non-finite values, `Analyse()` stops; SimDesign then draws a new
   data set and counts the error in the `ERRORS` column.
2. For each method: `fungible::BiFAD(R, facMethod = m, numFactors = NFAC)`. `numFactors` is the number of
   group factors; BiFAD adds the general factor.
3. A method that fails, or returns a loading matrix of the wrong size or with non-finite values, is
   recorded as `NA` for that replication, and the data set is kept. Redrawing it would hide that method's
   failures and inflate its convergence rate. `control = list(allow_na = TRUE)` lets SimDesign accept these
   `NA`s.
4. The estimate is aligned to the true loadings (`align_loadings()`).
5. For each method, `Analyse()` returns `CC_general`, `CC_specific` and the aligned salient loadings
   (`g1, g2, …` general, `s1, s2, …` group), prefixed with the method name: `ULS.CC_general`, `ULS.g1`, ….

### Alignment (`align_loadings()`)

Factor order and sign are arbitrary in EFA, so every estimate is matched to the true matrix first:

- General factor: the estimated column with the highest absolute Tucker congruence with the true general
  column, reflected to point in the same direction.
- Group factors: the remaining columns are matched and reflected with
  `fungible::faAlign(..., MatchMethod = "CC")`, which matches on congruence and uses the Hungarian
  algorithm when the best matches are not unique.

## Metrics (`Summarise()`)

Computed for each method, separately for the general and the group (specific) loadings, over the
replications in which the method converged. With λₚ a salient true loading, λ̂ₚᵣ its aligned estimate in
replication r, P the number of salient loadings in the block and R the number of converged replications:

| Column | Definition |
|---|---|
| `conv_rate` | proportion of replications with a usable solution |
| `CC_general`, `CC_specific` | Tucker's congruence between each aligned estimated column and the true column, averaged over the block's factors and then over replications |
| `AB_general`, `AB_specific` | absolute bias: (1/P) Σₚ \| mean_r(λ̂ₚᵣ) − λₚ \| |
| `RB_general`, `RB_specific` | relative bias, as a proportion (−0.04 = −4%): (1/P) Σₚ (mean_r(λ̂ₚᵣ) − λₚ) / λₚ |
| `MAE_general`, `MAE_specific` | mean absolute error: (1/(P·R)) Σₚ Σᵣ \| λ̂ₚᵣ − λₚ \| |

Tucker's congruence for factor k: φₖ = Σᵢ λᵢₖ λ̂ᵢₖ / √(Σᵢ λᵢₖ² · Σᵢ λ̂ᵢₖ²).

- AB, RB and MAE use only salient true loadings (|λ| > 1e-8), since RB is undefined for zero loadings.
  CC uses whole columns, zeros included.
- AB is the bias of the mean estimate; MAE also contains sampling variability.
- AB and RB are computed with `SimDesign::bias()`.

## Output

- Wide: one row per condition with the design columns, SimDesign's columns (`REPLICATIONS`, `SIM_TIME`,
  `ERRORS`, `WARNINGS`, …) and 6 × 9 metric columns named `<METHOD>.<metric>`, e.g. `ULS.CC_general`,
  `PCA.RB_specific`.
- Long (`to_long()`): one row per condition × method, with a `Method` column and unprefixed metric
  columns — the format for ANOVA and plots.
- Error and warning messages: `SimExtract(res, what = "errors")`, `SimExtract(res, what = "warnings")`.

## Reproducibility

- Replications: both runs use `iseed = 1276149341` (`genSeeds(Design, iseed = …)` locally, the `iseed`
  argument of `runArraySimulation()` on the cluster), which gives each condition its own L'Ecuyer-CMRG
  stream. Keep `iseed` unchanged between runs and resubmissions.
- Population loadings: `sim_factor()` draws them at random, so `population_model()` seeds it with
  `set.seed(seed_base + POP_ID)` (Mersenne-Twister) inside `with_preserved_rng()`, which restores the
  simulation's random-number state afterwards. As a result `Generate()` and `Summarise()` rebuild the same
  matrix, conditions that share a population use the same loadings, and the replication streams are
  untouched.
- `make_fixed_objects("row")` seeds with `seed_base + ROW` instead (1234 + i, as in the original data
  generation script) and reproduces the original population matrices. Each condition then has its own
  loading draw, so NOBS, LSKEW and NCAT effects are mixed with differences between population draws.
- `fungible::monte1()` calls `set.seed()` internally. `monte1_latent()` draws that seed from the current
  replication stream and restores the stream afterwards, so every replication and condition gets
  different latent data.

## Functions in `functions.R`

| Function | Role |
|---|---|
| `make_design()` | design data frame with `ROW` and `POP_ID` |
| `make_fixed_objects()` | methods, loading slot, salience threshold, seed settings |
| `hpc_settings()` | `ROWS_PER_JOB`, `OUT_DIR` and file prefix shared by `run_hpc.R` and `collect_hpc.R` |
| `with_preserved_rng()` | evaluates an expression, then restores `.Random.seed` |
| `population_model()` | seeded `sim_factor()` → true loadings and population correlation matrix |
| `monte1_latent()` | non-normal latent data through `fungible::monte1()` |
| `factor_index()` | column indices of general and group factors; checks the matrix size |
| `salient_mask()` | which true loadings are salient |
| `tucker_cc()`, `cc_against()` | Tucker's congruence |
| `extract_lambda()` | loading matrix from a BiFAD fit |
| `align_loadings()` | factor matching and reflection |
| `rep_output()` | one method's output for one replication |
| `Generate()`, `Analyse()`, `Summarise()` | SimDesign steps |
| `summarise_block()` | AB, RB and MAE for one block of loadings |
| `to_long()` | wide results → one row per condition × method |
