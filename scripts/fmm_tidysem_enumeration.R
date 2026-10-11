# =============================================================================
# Faktör Karışım Modeli — 1–5 sınıflı çözümlerin karşılaştırılması
# OpenMx (kestirim motoru) + tidySEM (model kurma, başlangıç değerleri,
# uyum tablosu, sınıf olasılıkları, VLMR/LMR ve BLRT)
#
# Gereksinimler: OpenMx >= 2.21.13, tidySEM >= 0.2.12, MASS
#   install.packages(c("OpenMx", "tidySEM", "MASS"))
#
# Veri üretimi fmm_openmx_enumeration.R ile birebir aynıdır (aynı seed ->
# aynı veri). O script'in bulduğu -2LL değerleri aşağıda referans olarak
# kullanılır; iki script aynı modeli kestirdiği için değerler eşleşmelidir.
# =============================================================================

library(OpenMx)
library(tidySEM)
library(MASS)

# -----------------------------------------------------------------------------
# 0. Ayarlar
# -----------------------------------------------------------------------------
seed       <- 20261011
classes    <- 1:5     # karşılaştırılacak sınıf sayıları
nStarts    <- 20      # k > 1 için ek başlangıç sayısı (k-means/hclust/rastgele)
runLMR     <- TRUE    # VLMR/LMR (tidySEM::lr_lmr; nonnest2 paketi gerekir)
runBLRT    <- TRUE    # BLRT uzun sürer; hızlı deneme için FALSE yap
blrtMaxK   <- 3       # BLRT yalnızca 1..blrtMaxK sınıflı modeller arasında
                      # (1 vs 2, 2 vs 3); 5'e kadar çıkmak saatler sürebilir
nBoot      <- 100     # BLRT bootstrap replikasyonu
nWorkers   <- max(1, parallel::detectCores() - 1)
outDir     <- "fmm_tidysem_output"

# fmm_openmx_enumeration.R'nin aynı veride bulduğu en iyi -2LL değerleri
refM2LL <- c(`1` = 25276.774, `2` = 25272.653, `3` = 25265.943,
             `4` = 25263.941, `5` = 25258.667)

# -----------------------------------------------------------------------------
# 1. Veri üretimi (popülasyon: 2 sınıf, 2 faktör, faktör başına 5 madde)
# -----------------------------------------------------------------------------
set.seed(seed)

N        <- 1000
pi2      <- 0.20
nFac     <- 2
nPerFac  <- 5
nItem    <- nFac * nPerFac
lambda   <- 0.70
facCor   <- 0.40

items    <- paste0("y", 1:nItem)
facNames <- paste0("F", 1:nFac)

Lambda <- matrix(0, nItem, nFac, dimnames = list(items, facNames))
for (f in 1:nFac) Lambda[((f - 1) * nPerFac + 1):(f * nPerFac), f] <- lambda

Phi   <- matrix(c(1, facCor, facCor, 1), 2, 2, dimnames = list(facNames, facNames))
Theta <- diag(1 - lambda^2, nItem)
tau   <- rep(0, nItem)
kappa <- list(class1 = c(0, 0), class2 = c(1, 1))

trueClass <- rbinom(N, 1, pi2) + 1L
eta <- matrix(NA_real_, N, nFac)
for (k in 1:2) {
  idx <- which(trueClass == k)
  eta[idx, ] <- mvrnorm(length(idx), mu = kappa[[k]], Sigma = Phi)
}
eps <- mvrnorm(N, mu = rep(0, nItem), Sigma = Theta)
Y   <- sweep(eta %*% t(Lambda), 2, tau, "+") + eps
colnames(Y) <- items
simData <- as.data.frame(Y)

cat("Gerçek sınıf dağılımı:\n"); print(table(trueClass))

# -----------------------------------------------------------------------------
# 2. Model sözdizimi (lavaan; {C} sınıf numarasıyla değiştirilir)
#
# Tanımlama — tüm sınıflar simetrik olsun diye:
#   - faktör varyansları 1 (her sınıfta), faktör korelasyonu serbest ve eşit
#   - yükler, hata varyansları ve kesişimler sınıflar arası eşit (ortak etiket)
#   - işaretçi maddelerin (y1, y6) kesişimi 0; diğer kesişimler serbest
#   - faktör ortalamaları HER sınıfta serbest (k{C}F1, k{C}F2)
# Tek bir {C} sözdizimi her sınıfa aynı kalıbı uyguladığı için "yalnızca
# sınıf 1'in ortalaması 0" doğrudan yazılamaz. Bu tanımlama, "sınıf 1
# ortalamaları 0, kesişimler serbest" modeliyle AYNI olabilirliği verir
# (tau' = tau + Lambda*kappa_1, kappa'_c = kappa_c - kappa_1 dönüşümü).
# Raporlamada sınıflar büyüklüğe göre sıralanır ve ortalamalar en büyük
# sınıfa göre verilir; yani en büyük sınıf "sınıf 1 = 0" referansıdır.
# -----------------------------------------------------------------------------
markers <- c("y1", "y6")
loadLine <- function(f) {
  its <- items[Lambda[, f] != 0]
  paste0(facNames[f], " =~ NA*", its[1], " + ",
         paste0("lam_", its, "*", its, collapse = " + "))
}
fmmSyntax <- paste(c(
  sapply(1:nFac, loadLine),
  paste0(facNames, " ~~ 1*", facNames),
  "F1 ~~ phi21*F2",
  paste0(items, " ~~ th_", items, "*", items),
  paste0(markers, " ~ 0*1"),
  paste0(setdiff(items, markers), " ~ tau_", setdiff(items, markers), "*1"),
  paste0(facNames, " ~ k{C}", facNames, "*1")
), collapse = "\n")
cat("\nModel sözdizimi (sınıf {C}):\n", fmmSyntax, "\n")

# Beklenen serbest parametre sayısı:
# yük 10 + phi 1 + hata 10 + kesişim 8 + faktör ortalaması 2k + ağırlık (k-1)
expectedNpar <- function(k) nItem + 1 + nItem + (nItem - length(markers)) +
  nFac * k + (k - 1)

checkModel <- function(m, k) {
  np <- length(omxGetParameters(m))
  if (np != expectedNpar(k))
    stop(sprintf(paste("k = %d: %d serbest parametre bulundu, %d bekleniyordu.",
                       "tidySEM sözdizimini beklenenden farklı yorumlamış olabilir;",
                       "omxGetParameters() çıktısını kontrol et."),
                 k, np, expectedNpar(k)))
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# 3. Kestirim
#   3a. tidySEM varsayılanı: k-means ile başlangıç + run_mx
#   3b. Ek başlangıçlar: farklı bölünmelerle mixture_starts + run_mx;
#       daha iyi -2LL bulunursa 3a sonucunun yerine konur.
# -----------------------------------------------------------------------------
m2ll <- function(fit) fit$output$minimum
step <- function(txt) cat(sprintf("\n[%s] ADIM: %s\n", format(Sys.time(), "%H:%M:%S"), txt))

cat(sprintf("\nR %s | OpenMx %s | tidySEM %s\n", getRversion(),
            packageVersion("OpenMx"), packageVersion("tidySEM")))

step("3a mx_mixture (tidySEM varsayılanları)")
set.seed(seed + 1)
res <- mx_mixture(model = fmmSyntax, classes = classes, data = simData, run = TRUE)
cat("mx_mixture sınıfı:", paste(class(res), collapse = ", "), "\n")
for (i in seq_along(classes)) checkModel(res[[i]], classes[i])

makeSplits <- function(k, n) {
  X <- scale(simData)
  sp <- list(hclust = cutree(hclust(dist(X), method = "ward.D2"), k = k))
  for (s in seq_len(n - 1)) {
    sp[[paste0("kmeans_", s)]] <- if (s %% 2 == 1) {
      kmeans(X, centers = k)$cluster
    } else {
      sample.int(k, nrow(X), replace = TRUE)
    }
  }
  sp
}

startLog <- list()
oracle   <- NULL
for (i in which(classes > 1)) {
  k    <- classes[i]
  step(sprintf("3b ek başlangıçlar, k = %d (%d başlangıç)", k, nStarts))
  base <- mx_mixture(model = fmmSyntax, classes = k, data = simData, run = FALSE)
  checkModel(base, k)

  splits <- makeSplits(k, nStarts)
  fits <- lapply(splits, function(sp)
    tryCatch(suppressWarnings(run_mx(mixture_starts(base, splits = sp))),
             error = function(e) NULL))
  ll <- sapply(fits, function(f) if (is.null(f)) NA_real_ else m2ll(f))

  startLog[[paste0("k", k)]] <- data.frame(
    k = k, start = names(splits), m2LL = round(ll, 3), row.names = NULL)

  best <- which.min(ll)
  if (length(best) && ll[best] < m2ll(res[[i]]) - 1e-3) {
    cat(sprintf("k = %d: ek başlangıç daha iyi çözüm buldu (%.3f < %.3f)\n",
                k, ll[best], m2ll(res[[i]])))
    res[[i]] <- fits[[best]]
  }
  cat(sprintf("k = %d: en iyi -2LL'yi tekrarlayan başlangıç sayısı: %d / %d\n",
              k, sum(abs(ll - min(ll, na.rm = TRUE)) < 1e-3, na.rm = TRUE),
              sum(!is.na(ll))))

  # Yalnızca simülasyon tanısı: gerçek sınıflardan başlatılan çözüm.
  # Gerçek veride bu bilgi yoktur; seçim için KULLANILMAZ.
  if (k == 2) {
    oracle <- tryCatch(suppressWarnings(run_mx(mixture_starts(base, splits = trueClass))),
                       error = function(e) NULL)
    if (!is.null(oracle))
      cat(sprintf("k = 2: gerçek sınıflardan başlatılan çözüm -2LL = %.3f (seçilen: %.3f)\n",
                  m2ll(oracle), m2ll(res[[i]])))
  }
}

# Ara kayıt: sonraki adımlar hata verirse kestirimler kaybolmasın
dir.create(outDir, showWarnings = FALSE)
saveRDS(list(data = simData, trueClass = trueClass, fits = res, oracle = oracle),
        file.path(outDir, "fits.rds"))

# OpenMx script'iyle tutarlılık kontrolü
refTab <- data.frame(classes = classes, tidySEM = round(sapply(res, m2ll), 3),
                     OpenMx = unname(refM2LL[as.character(classes)]))
refTab$fark <- round(refTab$tidySEM - refTab$OpenMx, 3)
refTab$yorum <- ifelse(is.na(refTab$fark), "",
                ifelse(refTab$fark >  0.01, "tidySEM yerel maksimumda: nStarts'ı artır",
                ifelse(refTab$fark < -0.01, "OpenMx script'i daha iyi çözümü kaçırmış",
                       "eşleşiyor")))
cat("\n-2LL tutarlılık kontrolü (fmm_openmx_enumeration.R ile):\n")
print(refTab, row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Uyum ve sınıf sayısı kararı
# -----------------------------------------------------------------------------
step("4 table_fit")
fitTab <- table_fit(res)
cat("\nUyum tablosu (tidySEM::table_fit):\n"); print(fitTab)

# VLMR/LMR'nin k-1 vs k karşılaştırmasındaki kuramsal temeli tartışmalıdır
# (Jeffries, 2003); karar için BIC ve BLRT'yi esas al.
lmrTab <- NULL
if (runLMR) {
  step("4 lr_lmr (VLMR/LMR)")
  lmrTab <- tryCatch(lr_lmr(res), error = function(e) {
    message("lr_lmr hata verdi: ", conditionMessage(e)); NULL })
  print(lmrTab)
}

blrtTab <- NULL
if (runBLRT) {
  keep <- which(classes <= blrtMaxK)
  resB <- res[keep]
  class(resB) <- class(res)   # `[` sınıf özniteliğini düşürür
  step(sprintf("4 BLRT (k = %s; %d replikasyon, %d çekirdek) — uzun sürebilir",
               paste(classes[keep], collapse = ","), nBoot, nWorkers))
  if (requireNamespace("future", quietly = TRUE) && nWorkers > 1)
    future::plan(future::multisession, workers = nWorkers)
  set.seed(seed + 2)
  blrtTab <- tryCatch(BLRT(resB, replications = nBoot), error = function(e) {
    message("BLRT hata verdi: ", conditionMessage(e)); NULL })
  if (requireNamespace("future", quietly = TRUE)) future::plan(future::sequential)
  print(blrtTab)
}

# -----------------------------------------------------------------------------
# 5. Sınıf olasılıkları, gerçek sınıflarla karşılaştırma, faktör ortalamaları
# -----------------------------------------------------------------------------
posteriorMatrix <- function(fit, k) {
  cp <- class_prob(fit, type = "individual")
  if (is.list(cp) && !is.data.frame(cp)) cp <- cp[["individual"]]
  cp <- as.matrix(cp)
  probCols <- grep("^CPROB", colnames(cp))
  if (length(probCols) != k) probCols <- seq_len(k)
  cp[, probCols, drop = FALSE]
}

step("5 sınıf olasılıkları ve sınıf ortalamaları")
classSummary <- list()
for (i in which(classes > 1)) {
  k   <- classes[i]
  fit <- res[[i]]
  cat(sprintf("\n================ %d sınıflı çözüm ================\n", k))
  print(class_prob(fit, type = c("sum.posterior", "avg.mostlikely")))

  post <- posteriorMatrix(fit, k)
  size <- colMeans(post)
  ord  <- order(size, decreasing = TRUE)   # yeni sınıf 1 = en büyük sınıf
  post <- post[, ord, drop = FALSE]
  size <- size[ord]

  modal <- max.col(post)
  cat("\nGerçek sınıf x atanan sınıf (sınıflar büyüklüğe göre sıralı):\n")
  print(table(true = trueClass, assigned = factor(modal, levels = 1:k)))

  p   <- omxGetParameters(fit)
  kap <- t(sapply(1:k, function(c) p[paste0("k", c, facNames)]))[ord, , drop = FALSE]
  # En büyük sınıfa göre fark = "sınıf 1 ortalamaları 0" parametrizasyonu;
  # faktör varyansları 1 olduğundan standart fark olarak okunur.
  kapRel <- sweep(kap, 2, kap[1, ], "-")

  cs <- data.frame(k = k, class = 1:k, size = round(size, 3),
                   nModal = tabulate(modal, nbins = k),
                   round(kapRel, 3), row.names = NULL)
  names(cs)[5:(4 + nFac)] <- paste0("mean_", facNames)
  cs$flag <- ifelse(cs$size < 0.05 | cs$nModal < 20, "küçük/boş sınıf", "")
  classSummary[[paste0("k", k)]] <- cs
  cat("\nSınıf büyüklükleri ve faktör ortalamaları (sınıf 1 = 0 referansı):\n")
  print(cs, row.names = FALSE)
}
cat("\nPopülasyon: sınıf oranları .80/.20, sınıf 2 - sınıf 1 ortalama farkı 1 / 1\n")

# Ölçme parametreleri (2 sınıflı çözüm) — popülasyonla karşılaştırma
i2 <- which(classes == 2)
if (length(i2)) {
  p2 <- omxGetParameters(res[[i2]])
  meas <- p2[!grepl("^k[0-9]", names(p2)) & names(p2) %in%
               c(paste0("lam_", items), "phi21", paste0("th_", items),
                 paste0("tau_", items))]
  truth <- c(setNames(rep(lambda, nItem), paste0("lam_", items)), phi21 = facCor,
             setNames(rep(1 - lambda^2, nItem), paste0("th_", items)),
             setNames(tau, paste0("tau_", items)))[names(meas)]
  cat("\n2 sınıflı çözüm — ölçme parametreleri:\n")
  print(data.frame(param = names(meas), true = truth, est = round(meas, 3),
                   bias = round(meas - truth, 3), row.names = NULL))
}

# -----------------------------------------------------------------------------
# 6. Özet ve çıktıların kaydı
# -----------------------------------------------------------------------------
summ <- data.frame(classes = classes,
                   m2LL = sapply(res, m2ll),
                   npar = sapply(res, function(f) length(omxGetParameters(f))))
summ$BIC  <- summ$m2LL + log(N) * summ$npar
summ$aBIC <- summ$m2LL + log((N + 2) / 24) * summ$npar
cat("\nÖzet (kendi hesabımla, table_fit ile tutarlı olmalı):\n")
print(summ, row.names = FALSE)
cat(sprintf("BIC'e göre seçilen sınıf sayısı: %d (popülasyon: 2)\n",
            summ$classes[which.min(summ$BIC)]))

write.csv(refTab, file.path(outDir, "m2ll_check.csv"), row.names = FALSE)
write.csv(as.data.frame(fitTab), file.path(outDir, "fit_table.csv"), row.names = FALSE)
write.csv(do.call(rbind, startLog), file.path(outDir, "starts.csv"), row.names = FALSE)
write.csv(do.call(rbind, classSummary), file.path(outDir, "class_summary.csv"), row.names = FALSE)
if (!is.null(lmrTab))  write.csv(as.data.frame(lmrTab),  file.path(outDir, "lmr.csv"))
if (!is.null(blrtTab)) write.csv(as.data.frame(blrtTab), file.path(outDir, "blrt.csv"))
saveRDS(list(data = simData, trueClass = trueClass, fits = res, oracle = oracle),
        file.path(outDir, "fits.rds"))
cat("\nÇıktılar:", normalizePath(outDir), "\n")
