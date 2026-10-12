# Mplus MONTECARLO (generate q3 J10 K3 N562) karşılığı
library(MASS)

lam <- c(0.3254, 0.3004, 0.3074, 0.3518, 0.4173, 0.4898, 0.4750, 0.3811, 0.4460, 0.4445,
         0.3767, 0.3755, 0.4860, 0.3294, 0.4492, 0.3163, 0.4481, 0.3375, 0.4267, 0.4227,
         0.3488, 0.4875, 0.4232, 0.4923, 0.4276, 0.3103, 0.4407, 0.3728, 0.3053, 0.4056)
th  <- c(0.8941, 0.9098, 0.9055, 0.8763, 0.8258, 0.7601, 0.7744, 0.8547, 0.8011, 0.8025,
         0.8581, 0.8590, 0.7638, 0.8915, 0.7982, 0.8999, 0.7992, 0.8861, 0.8179, 0.8214,
         0.8783, 0.7624, 0.8209, 0.7577, 0.8172, 0.9037, 0.8058, 0.8610, 0.9068, 0.8355)
nFac <- 3; nPer <- 10; J <- nFac * nPer
items <- paste0("y", 1:J)
L <- matrix(0, J, nFac); for (f in 1:nFac) L[(f - 1) * nPer + 1:nPer, f] <- lam[(f - 1) * nPer + 1:nPer]
Phi <- matrix(0.4987, nFac, nFac); diag(Phi) <- 1
kap <- rbind(0, 1.1568, 2.3136)[, rep(1, nFac)]   # sınıf 1..3 faktör ortalamaları
logit <- c(0, 0, 0)                                # [c#1*0 c#2*0], c#3 = 0

genK3 <- function(N = 562, seed = 47631781) {
  set.seed(seed)
  cls <- sample(1:3, N, replace = TRUE, prob = exp(logit) / sum(exp(logit)))
  eta <- kap[cls, ] + mvrnorm(N, rep(0, nFac), Phi)
  Y   <- eta %*% t(L) + mvrnorm(N, rep(0, J), diag(th))   # kesişimler 0
  colnames(Y) <- items
  list(data = as.data.frame(Y), class = cls)
}

sim <- genK3()
write.table(sim$data, "data.dat", row.names = FALSE, col.names = FALSE)  # Mplus için
table(sim$class)
