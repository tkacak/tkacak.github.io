## ===========================================================================
## functions.R
## SimDesign bilesenleri: Design, Generate, Analyse, Summarise + yardimcilar.
## run_local.R ve run_hpc.R bu dosyayi source() eder; tek kaynak burasidir.
##
## Metrikler (her yontem icin, genel / ozel faktor ayri):
##   CC  : Tucker uyum katsayisi (hizalanmis kestirim vs gercek yuk sutunu)
##   AB  : Mutlak yanlilik = ortalama_p | mean_r(lambda_hat_pr) - lambda_p |
##   RB  : Goreli yanlilik (%) = ortalama_p [ (mean_r(lambda_hat_pr) - lambda_p) / lambda_p ] * 100
##   MAE : ortalama_p mean_r | lambda_hat_pr - lambda_p |  (eski metrics.R'deki "AB")
##   conv_rate : yontemin gecerli cozum dondurdugu replikasyon orani
## AB/RB/MAE yalniz salient (sifir olmayan) gercek yukler uzerinden hesaplanir.
## ===========================================================================

## --- Tasarim ----------------------------------------------------------------
## expand.grid kullaniliyor: satir numaralari eski calismadaki kosul
## numaralariyla (sim_data_current_condition_<i>.RDS) birebir ayni kalir.
make_design <- function() {
  d <- expand.grid(
    NGEN  = 1,
    NFAC  = c(2, 3, 4),
    NVAR  = c(4, 8),
    GLOAD = c("low", "medium", "high"),
    FLOAD = c("low", "medium", "high"),
    FRHO  = c(0.00, 0.30, 0.50),
    GRHO  = 0.00,
    NOBS  = c(250, 500, 1000),
    CROSS = 0.00,
    LSKEW = c("normal", "slightly", "moderately"), # https://doi.org/10.3758/s13428-015-0619-7
    NCAT  = c(2, 3, 4, 5, 6, 7),
    stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  d$ROW <- seq_len(nrow(d))
  ## Ayni evren parametrelerini paylasan kosullar ayni POP_ID'yi alir
  ## (NOBS, LSKEW, NCAT evreni degistirmez).
  pop_key  <- do.call(paste, d[c("NGEN", "NFAC", "NVAR", "GLOAD", "FLOAD",
                                 "FRHO", "GRHO", "CROSS")])
  d$POP_ID <- match(pop_key, unique(pop_key))
  d
}

## --- Sabit nesneler (fixed_objects) -----------------------------------------
make_fixed_objects <- function(pop_seed = c("population", "row")) {
  list(
    ## isim = cikti oneki, deger = fungible::BiFAD facMethod
    methods   = c(ULS = "fals", ML = "faml", PA = "fapa",
                  REGULS = "faregLS", REGML = "faregML", PCA = "pca"),
    loading_slot = "BstarSL",   # BiFAD ciktisindan alinacak yuk matrisi
    min_true  = 1e-8,           # |lambda| bunun ustundeyse salient
    seed_base = 1234L,
    ## "population": evren seed'i = seed_base + POP_ID (onerilen)
    ## "row"       : evren seed'i = seed_base + ROW (eski kodu birebir yeniden uretir)
    pop_seed  = match.arg(pop_seed)
  )
}

## --- RNG yardimcilari -------------------------------------------------------
## expr'i calistirir, sonra .Random.seed'i (ve dolayisiyla RNGkind'i) eski
## haline getirir. SimDesign'in replikasyon akisi bozulmaz.
with_preserved_rng <- function(expr) {
  had  <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old  <- if (had) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(if (had) assign(".Random.seed", old, envir = .GlobalEnv)
          else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  expr
}

## --- Evren modeli -----------------------------------------------------------
## sim_factor yukleri rastgele cektigi icin sabit bir Mersenne-Twister seed'i ile
## deterministik uretilir; Generate ve Summarise ayni matrisi gorur.
population_model <- function(condition, fixed_objects) {
  id <- if (fixed_objects$pop_seed == "row") condition$ROW else condition$POP_ID
  with_preserved_rng({
    set.seed(fixed_objects$seed_base + id, kind = "Mersenne-Twister",
             normal.kind = "Inversion", sample.kind = "Rejection")
    model <- bifactor::sim_factor(
      n_generals         = condition$NGEN,
      groups_per_general = condition$NFAC,
      items_per_group    = condition$NVAR,
      loadings_g         = condition$GLOAD,
      loadings_s         = condition$FLOAD,
      crossloadings      = condition$CROSS,
      generals_rho       = condition$GRHO,
      groups_rho         = condition$FRHO,
      method             = "minres")
  })
  lambda <- as.matrix(model$lambda); storage.mode(lambda) <- "double"
  list(lambda = lambda, R = model$R)
}

## fungible::monte1 kendi icinde set.seed() cagirir. Seed SimDesign akisindan
## cekilir (her replikasyon/kosul farkli), sonra akis geri yuklenir.
## Eski kodda seed = 1234 + j idi: ayni j ayni ham sayilari tum kosullarda
## tekrar kullaniyordu.
monte1_latent <- function(nsub, R, skew, kurt) {
  s <- sample.int(.Machine$integer.max, 1L)
  p <- nrow(R)
  with_preserved_rng(
    fungible::monte1(seed = s, nvar = p, nsub = nsub, cormat = R,
                     skewvec = rep(skew, p), kurtvec = rep(kurt, p))$data)
}

## --- Faktor indeksleri ------------------------------------------------------
factor_index <- function(condition, TL) {
  n_gen  <- as.integer(condition$NGEN)
  n_spec <- n_gen * as.integer(condition$NFAC)
  n_item <- n_spec * as.integer(condition$NVAR)
  if (ncol(TL) != n_gen + n_spec || nrow(TL) != n_item)
    stop(sprintf("true lambda %dx%d, beklenen %dx%d", nrow(TL), ncol(TL),
                 n_item, n_gen + n_spec))
  list(g = seq_len(n_gen), s = seq.int(n_gen + 1L, n_gen + n_spec))
}

salient_mask <- function(TL, cols, min_true) abs(TL[, cols, drop = FALSE]) > min_true

## --- Hizalama ve Tucker CC --------------------------------------------------
tucker_cc <- function(A, B) {
  A <- as.matrix(A); B <- as.matrix(B)
  colSums(A * B) / sqrt(colSums(A^2) * colSums(B^2))
}

cc_against <- function(target, M) abs(tucker_cc(matrix(target, nrow(M), ncol(M)), M))

extract_lambda <- function(fit, slot) {
  lam <- if (is.list(fit)) fit[[slot]] else NULL
  if (is.list(lam) && !is.null(lam$lambda)) lam <- lam$lambda
  if (!is.matrix(lam)) return(NULL)
  unname(lam)
}

## Genel faktor(ler): congruence ile greedy eslestirme + isaret duzeltme.
## Ozel faktorler: fungible::faAlign (varsa MatchMethod = "CC").
align_loadings <- function(lam, TL, g_cols, s_cols) {
  if (is.null(lam) || !identical(dim(lam), dim(TL)) || !all(is.finite(lam)))
    return(NULL)
  avail  <- seq_len(ncol(lam))
  g_pick <- integer(length(g_cols))
  for (j in seq_along(g_cols)) {
    k <- which.max(cc_against(TL[, g_cols[j]], lam[, avail, drop = FALSE]))
    g_pick[j] <- avail[k]
    avail <- avail[-k]
  }
  g_est <- lam[, g_pick, drop = FALSE]
  sgn   <- sign(colSums(TL[, g_cols, drop = FALSE] * g_est))
  g_est <- sweep(g_est, 2, ifelse(sgn == 0, 1, sgn), `*`)

  args <- list(F1 = TL[, s_cols, drop = FALSE], F2 = lam[, avail, drop = FALSE])
  if ("MatchMethod" %in% names(formals(fungible::faAlign))) args$MatchMethod <- "CC"
  out_s <- tryCatch(do.call(fungible::faAlign, args), error = function(e) NULL)
  if (is.null(out_s) || is.null(out_s$F2)) return(NULL)

  lam_a <- cbind(g_est, out_s$F2)
  if (!identical(dim(lam_a), dim(TL))) return(NULL)
  unname(lam_a)
}

## --- Tek yontem, tek replikasyon --------------------------------------------
## Cikti uzunlugu kosula bagli ama sabit: 2 CC + salient genel + salient ozel.
rep_output <- function(lam_a, TL, idx, min_true) {
  mg <- salient_mask(TL, idx$g, min_true)
  ms <- salient_mask(TL, idx$s, min_true)
  if (is.null(lam_a)) {
    cc_g <- cc_s <- NA_real_
    eg <- rep(NA_real_, sum(mg)); es <- rep(NA_real_, sum(ms))
  } else {
    cc   <- tucker_cc(TL, lam_a)
    cc_g <- mean(cc[idx$g]); cc_s <- mean(cc[idx$s])
    eg   <- lam_a[, idx$g, drop = FALSE][mg]
    es   <- lam_a[, idx$s, drop = FALSE][ms]
  }
  c(CC_general = cc_g, CC_specific = cc_s,
    stats::setNames(eg, paste0("g", seq_along(eg))),
    stats::setNames(es, paste0("s", seq_along(es))))
}

## ===========================================================================
## SimDesign fonksiyonlari
## ===========================================================================
Generate <- function(condition, fixed_objects) {
  pop <- population_model(condition, fixed_objects)
  latent <- switch(condition$LSKEW,
    normal     = MASS::mvrnorm(condition$NOBS, rep(0, nrow(pop$R)), pop$R),
    slightly   = monte1_latent(condition$NOBS, pop$R, skew = 0.5, kurt = 1.5),
    moderately = monte1_latent(condition$NOBS, pop$R, skew = 1.5, kurt = 3),
    stop("Bilinmeyen LSKEW: ", condition$LSKEW))
  dat <- latentFactoR::categorize(data = latent, categories = condition$NCAT,
                                  skew_value = 0)
  list(dat = dat, lambda = pop$lambda)
}

## Polikorik matris replikasyon basina BIR KEZ hesaplanir, 6 yontem ayni
## matrisi kullanir. Polikorik basarisizsa stop() -> SimDesign veriyi yeniden
## ceker (ERRORS sutununda sayilir). Tek bir yontemin basarisizligi NA olarak
## kaydedilir; veri yeniden cekilmez (aksi halde conv_rate sisirilir).
Analyse <- function(condition, dat, fixed_objects) {
  R <- bifactor::polyfast(dat$dat)$correlation
  if (!is.matrix(R) || !all(is.finite(R))) stop("polyfast gecersiz korelasyon matrisi")
  TL  <- dat$lambda
  idx <- factor_index(condition, TL)
  n_group <- as.integer(condition$NGEN * condition$NFAC)

  out <- lapply(fixed_objects$methods, function(m) {
    fit <- tryCatch(
      fungible::BiFAD(R = R, B = NULL, facMethod = m, numFactors = n_group),
      error = function(e) NULL)
    lam_a <- align_loadings(extract_lambda(fit, fixed_objects$loading_slot),
                            TL, idx$g, idx$s)
    rep_output(lam_a, TL, idx, fixed_objects$min_true)
  })
  unlist(out)   # isimler: "ULS.CC_general", "ULS.g1", ..., "PCA.s32"
}

summarise_block <- function(est, truth) {
  if (!nrow(est)) return(c(AB = NA_real_, RB = NA_real_, MAE = NA_real_))
  c(AB  = mean(SimDesign::bias(est, parameter = truth, abs = TRUE)),
    RB  = mean(SimDesign::bias(est, parameter = truth, type = "relative",
                               percent = TRUE)),
    MAE = mean(abs(sweep(est, 2, truth))))
}

Summarise <- function(condition, results, fixed_objects) {
  TL  <- population_model(condition, fixed_objects)$lambda
  idx <- factor_index(condition, TL)
  tg  <- TL[, idx$g, drop = FALSE][salient_mask(TL, idx$g, fixed_objects$min_true)]
  ts  <- TL[, idx$s, drop = FALSE][salient_mask(TL, idx$s, fixed_objects$min_true)]
  res <- as.matrix(results)

  out <- lapply(names(fixed_objects$methods), function(m) {
    cc_g <- res[, paste0(m, ".CC_general")]
    cc_s <- res[, paste0(m, ".CC_specific")]
    ok   <- is.finite(cc_g) & is.finite(cc_s)
    eg   <- res[ok, paste0(m, ".g", seq_along(tg)), drop = FALSE]
    es   <- res[ok, paste0(m, ".s", seq_along(ts)), drop = FALSE]
    bg   <- summarise_block(eg, tg)
    bs   <- summarise_block(es, ts)
    c(conv_rate   = mean(ok),
      CC_general  = if (any(ok)) mean(cc_g[ok]) else NA_real_,
      CC_specific = if (any(ok)) mean(cc_s[ok]) else NA_real_,
      AB_general  = bg[["AB"]],  AB_specific  = bs[["AB"]],
      RB_general  = bg[["RB"]],  RB_specific  = bs[["RB"]],
      MAE_general = bg[["MAE"]], MAE_specific = bs[["MAE"]])
  })
  names(out) <- names(fixed_objects$methods)
  unlist(out)   # isimler: "ULS.CC_general", ..., "PCA.MAE_specific"
}

## --- Genis sonucu uzun formata cevir (ANOVA / grafik icin) ------------------
to_long <- function(res, methods = names(make_fixed_objects()$methods)) {
  res <- as.data.frame(res)
  pat <- paste0("^(", paste(methods, collapse = "|"), ")\\.")
  mcols <- grep(pat, names(res), value = TRUE)
  design_cols <- setdiff(names(res), mcols)
  do.call(rbind, lapply(methods, function(m) {
    cols <- grep(paste0("^", m, "\\."), mcols, value = TRUE)
    blk  <- res[cols]; names(blk) <- sub(paste0("^", m, "\\."), "", cols)
    cbind(res[design_cols], Method = m, blk, row.names = NULL)
  }))
}
