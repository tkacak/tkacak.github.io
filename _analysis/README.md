# Analysis scripts

Working R scripts kept alongside the site source. The folder name starts with `_`, so
`rmarkdown::render_site()` ignores it and nothing here is published to the website.

## `pisa2022_creative_thinking_eirm.R`

Explanatory IRT / GLMM analysis of the PISA 2022 Creative Thinking assessment
(English-language testing, `LANGTEST_COG == 313`), 18 dichotomised DT items.

### Model ladder

All four models share the random-effects structure `(1 | ID) + (1 | Items)`, so they are
properly nested and the likelihood ratio tests are valid.

| Model | Fixed effects |
|:--|:--|
| Model 0 | intercept only |
| Model 1 | + `Complexity_c` (ATOS item complexity, grand-mean centred) |
| Model 2 | + `PVREAD`, `Language`, `Homework`, `ESCS`, `Gender` |
| Model 3 | + `Complexity_c ×` each student-level predictor |

A separate `-1 + Items` Rasch model is fitted for descriptive item difficulties only; it is
not part of the ladder.

### Plausible values

Models 0 and 1 contain no plausible value and are fitted once. Models 2 and 3 are fitted
ten times, once per READ plausible value, and combined with Rubin's rules. The script
checks that all models are fitted on an identical N, which is what makes the AIC/BIC and
LRT comparisons across the ladder legitimate.

- **Fixed effects** — pooled with Rubin's rules (`rubin_pool_fixed()`), cross-checked
  against `mitml::testEstimates()`.
- **Likelihood ratio tests** — Model 0 → 1 is an ordinary LRT; Model 1 → 2 and Model 2 → 3
  use the D2 statistic (Li, Meng, Raghunathan & Rubin, 1991) via `mitml::testModels()`,
  falling back to `mice::D2()`.
- **AIC / BIC / logLik / R²** — reported as the arithmetic mean over the ten plausible
  value models. These are *not* Rubin-pooled; pooling likelihood-based fit indices is not
  defined. Tables label this explicitly and report the SD across plausible values.

### Weights

Models are fitted **unweighted**. Passing non-integer sampling weights to `glmer()` yields
a pseudo-likelihood, which would invalidate the AIC, BIC, logLik and LRT comparisons the
model ladder is built on. Normalised weights (`WEIGHT`, `WEIGHT_GLOBAL`, `WEIGHT_SENATE`)
are still computed and are used for the descriptive weighted correlations.

### Requirements

Required: `haven`, `tidyr`, `dplyr`, `lme4`, `purrr`, `tibble`, `broom.mixed`,
`performance`, `mice`, `mitml`, `ggplot2`.

Optional: `eirm` (used for estimation and its `print`/`plot` helpers when installed;
otherwise the script falls back to `lme4::glmer` with identical results), `sjPlot`
(HTML table), `expss` (label stripping; a base-R fallback is built in).

### Running

Put `STU_QQQ.SAV` and `CRT_COG.SAV` in the working directory and run:

```r
source("pisa2022_creative_thinking_eirm.R")
```

Output goes to `model_outputs_revised/` (full sample) and
`model_outputs_revised/EO/` (English-official subgroup analyses).
