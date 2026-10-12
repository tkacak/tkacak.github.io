# k3.inp karşılığı — OpenMx
library(OpenMx)
source("fmm_k3_generate.R")          # sim, items, L, nFac, nPer, J
mxOption(key = "Number of Threads", value = 1)

K       <- 3
nStarts <- 50
nCores  <- if (.Platform$OS.type == "windows") 1L else max(1L, parallel::detectCores() - 1L)
dat     <- sim$data
facs    <- paste0("f", 1:nFac)

lamLab <- matrix(NA, J, nFac); lamLab[L != 0] <- paste0("l", 1:J)
covLab <- matrix(NA, nFac, nFac)
for (a in 1:nFac) for (b in 1:nFac) if (a != b) covLab[a, b] <- paste0("c", min(a, b), max(a, b))

makeClass <- function(k) mxModel(paste0("class", k),
  mxMatrix("Full", J, nFac, free = L != 0, values = (L != 0) * 0.4, labels = lamLab, name = "L"),
  mxMatrix("Symm", nFac, nFac, free = !diag(nFac), values = ifelse(diag(nFac) == 1, 1, 0.3),
           labels = covLab, lbound = ifelse(diag(nFac) == 1, NA, -0.99),
           ubound = ifelse(diag(nFac) == 1, NA, 0.99), name = "P"),                      # f@1, f WITH f*
  mxMatrix("Diag", J, J, free = TRUE, values = 0.8, labels = paste0("e", 1:J), lbound = 1e-4, name = "Th"),
  mxMatrix("Full", 1, J, free = TRUE, values = 0, labels = paste0("i", 1:J), name = "Tau"),  # (i1-i30)
  mxMatrix("Full", 1, nFac, free = k != 1, values = 0,
           labels = if (k != 1) paste0("m", k, "_", facs) else NA, name = "K"),            # c#1 = 0
  mxAlgebra(L %*% P %*% t(L) + Th, dimnames = list(items, items), name = "S"),
  mxAlgebra(Tau + K %*% t(L), dimnames = list(NULL, items), name = "M"),
  mxExpectationNormal("S", "M"), mxFitFunctionML(vector = TRUE))

buildFMM <- function(K) {
  m <- mxModel(paste0("FMM_K", K), lapply(1:K, makeClass), mxData(dat, "raw"),
    mxMatrix("Full", 1, K, free = c(FALSE, rep(TRUE, K - 1)), values = 0,
             labels = c(NA, paste0("w", 2:K)), name = "weights"),
    mxExpectationMixture(paste0("class", 1:K), scale = "softmax"), mxFitFunctionML())
  mxOption(mxOption(m, "Calculate Hessian", "No"), "Standard Errors", "No")
}

fitFMM <- function(K, nStarts) {
  base <- buildFMM(K)
  fits <- parallel::mclapply(seq_len(nStarts), function(s) {
    st <- c(setNames(runif(nFac * (K - 1), -1, 3), unlist(lapply(2:K, function(k) paste0("m", k, "_", facs)))),
            setNames(runif(K - 1, -1, 1), paste0("w", 2:K)))
    tryCatch(mxRun(omxSetParameters(base, labels = names(st), values = st),
                   silent = TRUE, suppressWarnings = TRUE), error = function(e) NULL)
  }, mc.cores = nCores, mc.set.seed = TRUE)
  ll   <- sapply(fits, function(f) if (inherits(f, "MxModel")) f$output$minimum else NA)
  best <- fits[[which.min(ll)]]
  best <- mxRun(mxOption(mxOption(best, "Calculate Hessian", "Yes"), "Standard Errors", "Yes"),
                silent = TRUE)
  attr(best, "replicated") <- sum(abs(ll - min(ll, na.rm = TRUE)) < 0.01, na.rm = TRUE)
  best
}

set.seed(1)
fit <- fitFMM(K, nStarts)

w    <- fit$weights$values[1, ]; p <- exp(w) / sum(exp(w))
rl   <- sapply(1:K, function(k) fit[[paste0("class", k)]]$fitfunction$result[, 1])
post <- sweep(rl, 2, p, "*"); post <- post / rowSums(post)
m2   <- fit$output$minimum; np <- length(omxGetParameters(fit)); N <- nrow(dat)

cat(sprintf("LL = %.3f | Npar = %d | BIC = %.2f | aBIC = %.2f | entropi = %.3f | best LL tekrar: %d/%d\n",
            -m2 / 2, np, m2 + log(N) * np, m2 + log((N + 2) / 24) * np,
            1 - sum(-post * log(pmax(post, 1e-300))) / (N * log(K)),
            attr(fit, "replicated"), nStarts))
est <- omxGetParameters(fit)
print(data.frame(class = 1:K, prop = round(p, 3),
                 round(rbind(0, t(sapply(2:K, function(k) est[paste0("m", k, "_", facs)]))), 3)))
print(table(true = sim$class, assigned = max.col(post)))
summary(fit)$parameters[, c("name", "Estimate", "Std.Error")]
