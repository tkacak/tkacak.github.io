# =============================================================================
# Faktör Karışım Modeli (FMM) — aynı veride 1–5 sınıflı çözümlerin
# OpenMx ile kestirimi ve karşılaştırılması
#
# Popülasyon modeli
#   - 2 sınıf, 2 faktör, faktör başına 5 madde (toplam 10 madde)
#   - Faktör yükleri = .70 (basit yapı), kesişimler = 0 (her iki sınıfta)
#   - Hata varyansları = 1 - .70^2 = .51
#   - Faktörler arası korelasyon = .40 (her iki sınıfta)
#   - Sınıf 1 (referans): faktör ortalamaları 0, varyanslar 1
#   - Sınıf 2           : faktör ortalamaları 1, varyanslar 1
#   - P(sınıf 2) = .20, N = 1000
#
# Analiz modeli (K = 1, ..., 5)
#   - Yükler, hata varyansları, kesişimler ve faktör korelasyonu sınıflar
#     arası eşit (ortak etiket)
#   - Faktör varyansları her sınıfta 1'e sabit
#   - Sınıf 1'in faktör ortalamaları 0; diğer sınıfların ortalamaları serbest
#   - K = 1, tek sınıflı 2 faktörlü CFA'dır
#
# NOT: Clark ve ark. (2013) numaralandırmasında bu yapı FMM-2'ye karşılık
# gelir (FMM-1 orada sınıf içi faktör varyansı 0 olan modeldir).
# =============================================================================

library(OpenMx)
library(MASS)

# Başlangıçlar paralel koşulduğu için her kestirim tek iş parçacığı kullanır
mxOption(key = "Number of Threads", value = 1)

set.seed(20261011)

# -----------------------------------------------------------------------------
# 0. Ayarlar
# -----------------------------------------------------------------------------
Kmax    <- 5     # en fazla sınıf sayısı
nStarts <- 30    # K > 1 için başlangıç sayısı
nCores  <- if (.Platform$OS.type == "windows") 1L else max(1L, parallel::detectCores() - 1L)
outDir  <- "fmm_enumeration_output"

# -----------------------------------------------------------------------------
# 1. Popülasyon parametreleri
# -----------------------------------------------------------------------------
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

kappa <- list(class1 = c(0, 0),   # referans
              class2 = c(1, 1))

Sigma_pop <- Lambda %*% Phi %*% t(Lambda) + Theta
dMu       <- Lambda %*% (kappa$class2 - kappa$class1)
D_obs     <- sqrt(drop(t(dMu) %*% solve(Sigma_pop) %*% dMu))
cat(sprintf("Gözlenen düzeyde Mahalanobis D (sınıf ayrışması): %.3f\n", D_obs))

# -----------------------------------------------------------------------------
# 2. Veri üretimi
#    (mxGenerateData mxExpectationMixture için tanımlı değil; elle üretiyoruz)
# -----------------------------------------------------------------------------
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

cat("\nGerçek sınıf dağılımı:\n"); print(table(trueClass))

# -----------------------------------------------------------------------------
# 3. Model kurucular
# -----------------------------------------------------------------------------
lamLabels   <- ifelse(Lambda != 0, paste0("lam_", rep(items, nFac), "_",
                                          rep(facNames, each = nItem)), NA)
kappaLabels <- function(k) paste0("kappa", k, "_", facNames)
logitLabels <- function(K) paste0("logit_", 2:K)

makeClass <- function(k) {
  meanFree <- k != 1
  mxModel(paste0("class", k),
    mxMatrix("Full", nItem, nFac, free = Lambda != 0, values = Lambda,
             labels = lamLabels, dimnames = list(items, facNames), name = "L"),
    mxMatrix("Symm", nFac, nFac, free = c(FALSE, TRUE, FALSE),
             values = c(1, facCor, 1), labels = c(NA, "phi21", NA),
             lbound = c(NA, -0.99, NA), ubound = c(NA, 0.99, NA), name = "P"),
    mxMatrix("Diag", nItem, nItem, free = TRUE, values = 1 - lambda^2,
             labels = paste0("th_", items), lbound = 1e-4, name = "Th"),
    mxMatrix("Full", 1, nItem, free = TRUE, values = 0,
             labels = paste0("tau_", items), dimnames = list(NULL, items),
             name = "Tau"),
    mxMatrix("Full", 1, nFac, free = meanFree, values = 0,
             labels = if (meanFree) kappaLabels(k) else NA, name = "K"),
    mxAlgebra(L %*% P %*% t(L) + Th, dimnames = list(items, items), name = "S"),
    mxAlgebra(Tau + K %*% t(L), dimnames = list(NULL, items), name = "M"),
    mxExpectationNormal(covariance = "S", means = "M"),
    mxFitFunctionML(vector = TRUE)
  )
}

buildFMM <- function(K) {
  if (K == 1) {
    m <- mxModel(makeClass(1), name = "FMM_1class", mxFitFunctionML(),
                 mxData(simData, type = "raw"))
  } else {
    m <- mxModel(paste0("FMM_", K, "class"),
      lapply(1:K, makeClass),
      mxData(simData, type = "raw"),
      # softmax: sınıf 1 (referans) ham ağırlığı 0'a sabit
      mxMatrix("Full", K, 1, free = c(FALSE, rep(TRUE, K - 1)), values = 0,
               labels = c(NA, logitLabels(K)), name = "weights"),
      mxExpectationMixture(paste0("class", 1:K), scale = "softmax"),
      mxFitFunctionML())
  }
  # Çoklu başlangıçta Hessian hesaplanmaz; yalnızca en iyi çözüm için hesaplanır
  m <- mxOption(m, "Calculate Hessian", "No")
  mxOption(m, "Standard Errors", "No")
}

setStart <- function(model, values) {
  values <- values[names(values) %in% names(omxGetParameters(model))]
  omxSetParameters(model, labels = names(values), values = values)
}

# -----------------------------------------------------------------------------
# 4. Başlangıç değeri üreticileri
#   (a) k-means: madde ortalaması tabanlı faktör puanlarında kümeleme;
#       en büyük küme sınıf 1 olur, diğerlerinin ortalaması ona göre
#   (b) bir önceki K'nın en iyi çözümü + rastgele yeni sınıf
#   (c) 1 sınıflı çözüm + rastgele sınıf ortalamaları ve oranları
# -----------------------------------------------------------------------------
kmeansStart <- function(K, measBase) {
  lamBar <- mean(measBase[grep("^lam_", names(measBase))])
  sc <- sapply(1:nFac, function(f) rowMeans(simData[, Lambda[, f] != 0]) / lamBar)
  km  <- kmeans(sc, centers = K, nstart = 20)
  ord <- order(km$size, decreasing = TRUE)
  ctr <- km$centers[ord, , drop = FALSE]
  kap <- sweep(ctr, 2, ctr[1, ])
  tauStart <- colMeans(simData[km$cluster == ord[1], ])
  c(setNames(as.vector(t(kap[-1, , drop = FALSE])),
             unlist(lapply(2:K, kappaLabels))),
    setNames(log(km$size[ord][-1] / km$size[ord][1]), logitLabels(K)),
    setNames(tauStart, paste0("tau_", items)))
}

randomStart <- function(K) {
  c(setNames(runif(nFac * (K - 1), -2, 2), unlist(lapply(2:K, kappaLabels))),
    setNames(runif(K - 1, -2, 1), logitLabels(K)))
}

# -----------------------------------------------------------------------------
# 5. Kestirim
# -----------------------------------------------------------------------------
runOne <- function(model) {
  tryCatch(mxRun(model, silent = TRUE, suppressWarnings = TRUE),
           error = function(e) NULL)
}
m2llOf <- function(f) if (!inherits(f, "MxModel")) NA_real_ else f$output$minimum

fits     <- vector("list", Kmax)
startLog <- list()
t0 <- Sys.time()

for (K in 1:Kmax) {
  tk <- Sys.time()
  if (K == 1) {
    best <- runOne(buildFMM(1))
    measBase <- omxGetParameters(best)
    ll <- m2llOf(best)
  } else {
    base <- setStart(buildFMM(K), measBase)
    prev <- omxGetParameters(fits[[K - 1]])
    starts <- c(
      list(kmeans = kmeansStart(K, measBase)),
      list(previous = c(prev, setNames(runif(nFac, -2, 2), kappaLabels(K)),
                        setNames(-2, paste0("logit_", K)))),
      lapply(seq_len(nStarts - 2), function(s) randomStart(K))
    )
    names(starts)[-(1:2)] <- paste0("random_", seq_len(nStarts - 2))

    res <- parallel::mclapply(starts, function(st) runOne(setStart(base, st)),
                              mc.cores = nCores, mc.set.seed = TRUE)
    ll  <- sapply(res, m2llOf)
    startLog[[K]] <- data.frame(K = K, start = names(starts),
                                m2LL = round(ll, 3), row.names = NULL)
    best <- res[[which.min(ll)]]
  }

  # En iyi çözüm için Hessian ve standart hatalar
  best <- mxOption(best, "Calculate Hessian", "Yes")
  best <- mxOption(best, "Standard Errors", "Yes")
  best <- mxRun(best, silent = TRUE, suppressWarnings = TRUE)
  fits[[K]] <- best

  nRep <- sum(abs(ll - min(ll, na.rm = TRUE)) < 0.01, na.rm = TRUE)
  cat(sprintf("K = %d | -2LL = %.3f | en iyiyi tekrarlayan: %d/%d | süre: %.0f sn\n",
              K, best$output$minimum, nRep, sum(!is.na(ll)),
              as.numeric(difftime(Sys.time(), tk, units = "secs"))))
}
cat(sprintf("Toplam süre: %.1f dk\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# -----------------------------------------------------------------------------
# 6. Her çözüm için özet bilgiler
# -----------------------------------------------------------------------------
classInfo <- function(fit, K) {
  if (K == 1) {
    return(list(p = 1, post = matrix(1, N, 1),
                kap = matrix(0, 1, nFac, dimnames = list("class1", facNames))))
  }
  w <- fit$weights$values[, 1]
  p <- exp(w) / sum(exp(w))
  rowLik <- sapply(1:K, function(k) fit[[paste0("class", k)]]$fitfunction$result[, 1])
  post <- sweep(rowLik, 2, p, "*")
  post <- post / rowSums(post)
  est <- omxGetParameters(fit)
  kap <- rbind(0, t(sapply(2:K, function(k) est[kappaLabels(k)])))
  dimnames(kap) <- list(paste0("class", 1:K), facNames)
  list(p = p, post = post, kap = kap)
}

info <- lapply(1:Kmax, function(K) classInfo(fits[[K]], K))

fitTab <- do.call(rbind, lapply(1:Kmax, function(K) {
  f  <- fits[[K]]
  m2 <- f$output$minimum
  np <- length(omxGetParameters(f))
  ci <- info[[K]]
  ent <- if (K == 1) NA else
    1 - sum(-ci$post * log(pmax(ci$post, 1e-300))) / (N * log(K))
  modalN <- tabulate(max.col(ci$post), nbins = K)
  se <- tryCatch(summary(f)$parameters$Std.Error, error = function(e) NA)
  # Sınıflar arası en küçük ortalama farkı (aynılaşmış sınıfları yakalamak için)
  minDist <- if (K == 1) NA else min(dist(ci$kap))

  data.frame(K = K, m2LL = round(m2, 2), npar = np,
             AIC  = round(m2 + 2 * np, 2),
             BIC  = round(m2 + log(N) * np, 2),
             aBIC = round(m2 + log((N + 2) / 24) * np, 2),
             entropy = round(ent, 3),
             minProp = round(min(ci$p), 3),
             minModalN = min(modalN),
             minMeanDist = round(minDist, 3),
             status = f$output$status$code,
             SEok = all(is.finite(se)))
}))

fitTab$flag <- with(fitTab, ifelse(
  K > 1 & (minProp < 0.05 | minModalN < 20), "kucuk_sinif",
  ifelse(K > 1 & minMeanDist < 0.20, "ayni_siniflar",
         ifelse(!SEok, "SE_yok", ""))))

cat("\n================ Uyum karşılaştırması ================\n")
print(fitTab, row.names = FALSE)
cat(sprintf("\nBIC'e göre: K = %d | aBIC'e göre: K = %d | AIC'e göre: K = %d (popülasyon: K = 2)\n",
            fitTab$K[which.min(fitTab$BIC)], fitTab$K[which.min(fitTab$aBIC)],
            fitTab$K[which.min(fitTab$AIC)]))

# Ardışık modeller arasındaki olabilirlik oranı farkı.
# p-değeri bilerek verilmiyor: K-1 ile K sınıflı karışım karşılaştırmasında
# LR istatistiği ki-kare dağılmaz (parametre sınırda; düzenlilik koşulları
# sağlanmaz). Test gerekiyorsa BLRT kullan.
lrTab <- data.frame(comparison = paste0(1:(Kmax - 1), " vs ", 2:Kmax),
                    LR  = round(-diff(fitTab$m2LL), 2),
                    dfDiff = diff(fitTab$npar))
cat("\nArdışık -2LL farkları (p-değeri yok; yukarıdaki nota bak):\n")
print(lrTab, row.names = FALSE)

# -----------------------------------------------------------------------------
# 7. Sınıf büyüklükleri, faktör ortalamaları ve gerçek sınıflarla karşılaştırma
# -----------------------------------------------------------------------------
for (K in 2:Kmax) {
  ci <- info[[K]]
  cat(sprintf("\n================ K = %d ================\n", K))
  print(data.frame(class = 1:K, prop = round(ci$p, 3), round(ci$kap, 3)),
        row.names = FALSE)
  cat("Gerçek sınıf x atanan sınıf:\n")
  print(table(true = trueClass, assigned = factor(max.col(ci$post), levels = 1:K)))
}
cat("\nPopülasyon: sınıf oranları .80/.20, sınıf 2 - sınıf 1 ortalama farkı 1/1\n")

# Ölçme parametreleri (K = 2) — popülasyonla karşılaştırma
est2 <- omxGetParameters(fits[[2]])
se2  <- summary(fits[[2]])$parameters
meas <- est2[!grepl("^kappa|^logit", names(est2))]
truth <- c(setNames(rep(lambda, nItem), lamLabels[!is.na(lamLabels)]),
           phi21 = facCor,
           setNames(rep(1 - lambda^2, nItem), paste0("th_", items)),
           setNames(tau, paste0("tau_", items)))[names(meas)]
cat("\nK = 2 — ölçme parametreleri:\n")
print(data.frame(param = names(meas), true = truth, est = round(meas, 3),
                 bias = round(meas - truth, 3),
                 SE = round(se2$Std.Error[match(names(meas), se2$name)], 3),
                 row.names = NULL))

# -----------------------------------------------------------------------------
# 8. Kayıt
# -----------------------------------------------------------------------------
dir.create(outDir, showWarnings = FALSE)
write.csv(fitTab, file.path(outDir, "fit_table.csv"), row.names = FALSE)
write.csv(do.call(rbind, startLog), file.path(outDir, "starts.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(2:Kmax, function(K)
  data.frame(K = K, class = 1:K, prop = info[[K]]$p, info[[K]]$kap))),
  file.path(outDir, "class_means.csv"), row.names = FALSE)
saveRDS(list(data = simData, trueClass = trueClass, fits = fits, info = info),
        file.path(outDir, "fits.rds"))
cat("\nÇıktılar:", normalizePath(outDir), "\n")
