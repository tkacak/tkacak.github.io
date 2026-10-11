# =============================================================================
# Faktör Karışım Modeli (FMM) — OpenMx ile veri üretimi ve analiz
#
# Popülasyon modeli
#   - 2 sınıf, 2 faktör, faktör başına 5 madde (toplam 10 madde)
#   - Faktör yükleri = .70 (basit yapı), kesişimler = 0 (her iki sınıfta)
#   - Hata varyansları = 1 - .70^2 = .51
#   - Faktörler arası korelasyon = .40 (her iki sınıfta)
#   - Sınıf 1 (referans): faktör ortalamaları 0, varyanslar 1
#   - Sınıf 2           : faktör ortalamaları 1, varyanslar 1
#   - P(sınıf 2) = .30, N = 1000
#
# Ölçme değişmezliği tam (yük, kesişim, hata varyansı sınıflar arası eşit);
# sınıflar yalnızca faktör ortalamalarında ayrışır.
#
# NOT: Bu yapı Lubke & Neale numaralandırmasındaki "ortalama farkı + değişmez
# ölçme" modeline karşılık gelir. Clark ve ark. (2013) numaralandırmasında
# FMM-1, sınıf içi faktör varyansı 0'a sabit modeldir (LCFA); oradaki karşılığı
# FMM-2'dir. Raporlarken hangi taksonomiyi kullandığını belirt.
# =============================================================================

library(OpenMx)
library(MASS)

set.seed(20261011)

# -----------------------------------------------------------------------------
# 1. Popülasyon parametreleri
# -----------------------------------------------------------------------------
N        <- 1000
pi2      <- 0.30
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

# Kuramsal sınıf ayrışması (gözlenen düzeyde Mahalanobis D)
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
cat("\nSınıf içi madde ortalamaları:\n")
print(round(aggregate(simData, list(class = trueClass), mean), 2))

# -----------------------------------------------------------------------------
# 3. Analiz modelleri
#    Tüm ölçme parametreleri sınıflar arası aynı etiketle eşitlenir.
#    Tanımlama: sınıf 1 faktör ortalamaları 0; faktör varyansları 1 (her iki
#    sınıfta), faktör korelasyonu serbest ve eşit; sınıf 2 ortalamaları serbest.
#    fixLoadings = TRUE yaparsan yükler .70'e sabitlenir (popülasyon modeli).
# -----------------------------------------------------------------------------
fixLoadings <- FALSE

lamFree   <- Lambda != 0 & !fixLoadings
lamLabels <- ifelse(Lambda != 0, paste0("lam_", rep(items, nFac), "_",
                                        rep(facNames, each = nItem)), NA)

makeClass <- function(k, kappaStart = c(0, 0)) {
  meanFree <- k != 1
  mxModel(paste0("class", k),
    mxMatrix("Full", nItem, nFac, free = lamFree, values = Lambda,
             labels = lamLabels, dimnames = list(items, facNames), name = "L"),
    mxMatrix("Symm", nFac, nFac, free = c(FALSE, TRUE, FALSE),
             values = c(1, facCor, 1), labels = c(NA, "phi21", NA),
             lbound = c(NA, -0.99, NA), ubound = c(NA, 0.99, NA), name = "P"),
    mxMatrix("Diag", nItem, nItem, free = TRUE, values = 1 - lambda^2,
             labels = paste0("th_", items), lbound = 1e-4, name = "Th"),
    mxMatrix("Full", 1, nItem, free = TRUE, values = 0,
             labels = paste0("tau_", items), dimnames = list(NULL, items),
             name = "Tau"),
    mxMatrix("Full", 1, nFac, free = meanFree, values = kappaStart,
             labels = if (meanFree) paste0("kappa", k, "_", facNames) else NA,
             name = "K"),
    mxAlgebra(L %*% P %*% t(L) + Th, dimnames = list(items, items), name = "S"),
    mxAlgebra(Tau + K %*% t(L), dimnames = list(NULL, items), name = "M"),
    mxExpectationNormal(covariance = "S", means = "M"),
    mxFitFunctionML(vector = TRUE)
  )
}

# 3a. Tek sınıflı, 2 faktörlü CFA (karşılaştırma modeli)
cfa1 <- makeClass(1)
cfa1 <- mxModel(cfa1, name = "CFA_1class", mxFitFunctionML(),
                mxData(simData, type = "raw"))
cfa1Fit <- mxRun(cfa1, silent = TRUE)

# 3b. İki sınıflı FMM
buildFMM <- function(kappaStart, logitStart) {
  mxModel("FMM_2class",
    makeClass(1),
    makeClass(2, kappaStart),
    mxData(simData, type = "raw"),
    # softmax: sınıf 1 (referans) ham ağırlığı 0'a sabit
    mxMatrix("Full", 2, 1, free = c(FALSE, TRUE), values = c(0, logitStart),
             labels = c(NA, "logit_pi2"), name = "weights"),
    mxExpectationMixture(paste0("class", 1:2), scale = "softmax"),
    mxFitFunctionML()
  )
}

# Karışım modellerinde olabilirlik yüzeyi çok tepelidir: çoklu rastgele başlangıç
nStarts <- 30
fits <- vector("list", nStarts)
for (s in 1:nStarts) {
  ks <- if (s == 1) c(1, 1) else runif(2, -2, 2)
  ls <- if (s == 1) qlogis(pi2) else runif(1, -2, 1)
  fits[[s]] <- tryCatch(mxRun(buildFMM(ks, ls), silent = TRUE, suppressWarnings = TRUE),
                        error = function(e) NULL)
}
m2ll   <- sapply(fits, function(f) if (is.null(f)) NA else f$output$minimum)
codes  <- sapply(fits, function(f) if (is.null(f)) NA else f$output$status$code)
best   <- which.min(m2ll)
fmmFit <- fits[[best]]

cat(sprintf("\nBaşlangıçlar: %d | en iyi -2LL = %.3f | en iyi değeri tekrarlayan: %d\n",
            nStarts, m2ll[best], sum(abs(m2ll - m2ll[best]) < 1e-3, na.rm = TRUE)))
cat("Optimizer durum kodları:\n"); print(table(codes, useNA = "ifany"))

# Farklı yerel çözümleri göster (başlangıç 1 = popülasyon değerleri)
sol <- t(sapply(fits, function(f) {
  if (is.null(f)) return(rep(NA, 4))
  p <- omxGetParameters(f)
  c(m2LL = f$output$minimum, pi2 = plogis(p[["logit_pi2"]]),
    k1 = p[["kappa2_F1"]], k2 = p[["kappa2_F2"]])
}))
sol <- as.data.frame(round(sol, 3)); sol$start <- seq_len(nStarts)
sol <- sol[order(sol$m2LL), ]
cat("\nYerel çözümler (-2LL'ye göre, tekrarlar atıldı):\n")
print(sol[!duplicated(round(sol$m2LL, 2)), ], row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Sonuçlar
# -----------------------------------------------------------------------------
w    <- fmmFit$weights$values[, 1]
pHat <- exp(w) / sum(exp(w))   # softmax
cat("\nTahmin edilen sınıf oranları:", round(pHat, 3), "\n")

est <- omxGetParameters(fmmFit)
se  <- tryCatch(summary(fmmFit)$parameters[, c("name", "Std.Error")],
                error = function(e) NULL)

truth <- c(setNames(rep(lambda, nItem), lamLabels[!is.na(lamLabels)]),
           phi21 = facCor,
           setNames(rep(1 - lambda^2, nItem), paste0("th_", items)),
           setNames(tau, paste0("tau_", items)),
           kappa2_F1 = 1, kappa2_F2 = 1,
           logit_pi2 = qlogis(pi2))
truth <- truth[names(est)]
recovery <- data.frame(param = names(est), true = truth, est = round(est, 3),
                       bias = round(est - truth, 3), row.names = NULL)
if (!is.null(se)) recovery$SE <- round(se$Std.Error[match(names(est), se$name)], 3)
cat("\nParametre geri kazanımı:\n"); print(recovery)

# Etiket değişimi uyarısı: sınıf 2 ortalaması negatif ve pi2 > .5 ise
# sınıflar yer değiştirmiştir; yorum o yönde yapılmalı.
if (all(est[c("kappa2_F1", "kappa2_F2")] < 0))
  cat("UYARI: Etiket değişimi olası (sınıf 2 ortalamaları negatif).\n")

# Bireysel posterior sınıf olasılıkları
rowLik <- sapply(1:2, function(k) fmmFit[[paste0("class", k)]]$fitfunction$result[, 1])
post   <- sweep(rowLik, 2, pHat, "*")
post   <- post / rowSums(post)
modal  <- max.col(post)

entropy <- 1 - sum(-post * log(pmax(post, 1e-300))) / (N * log(2))
cat(sprintf("\nGöreli entropi: %.3f\n", entropy))
cat("Gerçek x atanan sınıf:\n"); print(table(true = trueClass, assigned = modal))
acc <- max(mean(modal == trueClass), mean((3 - modal) == trueClass))
cat(sprintf("Sınıflama doğruluğu (etiket değişimine dayanıklı): %.3f\n", acc))
cat("Ortalama posterior olasılık (atanan sınıfa göre):\n")
print(round(sapply(1:2, function(k) colMeans(post[modal == k, , drop = FALSE])), 3))

# Model karşılaştırma: 1 sınıflı CFA vs 2 sınıflı FMM
ic <- function(fit, label) {
  m2 <- fit$output$minimum
  k  <- length(omxGetParameters(fit))
  data.frame(model = label, m2LL = round(m2, 2), npar = k,
             AIC  = round(m2 + 2 * k, 2),
             BIC  = round(m2 + log(N) * k, 2),
             aBIC = round(m2 + log((N + 2) / 24) * k, 2))
}
cat("\nModel karşılaştırması:\n")
print(rbind(ic(cfa1Fit, "CFA_1class"), ic(fmmFit, "FMM_2class")))
