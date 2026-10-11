# =============================================================================
# Latent GOLD (6.x Syntax modülü) dosyaları — faktör karışım modeli
#
# fmm_openmx_enumeration.R'deki modelin LG karşılığını yazar:
#   - 2 sürekli faktör (F1: y1-y5, F2: y6-y10), basit yapı, bütün yükler serbest
#   - yükler, kesişimler, hata varyansları ve F1-F2 kovaryansı sınıflar arası eşit
#   - faktör varyansları her sınıfta 1
#   - faktör ortalamaları sınıf 1'de 0, diğer sınıflarda serbest
#   - Bayes sabitleri 0 (saf ML; OpenMx -2LL değerleriyle karşılaştırma için)
#
# Üretilen dosyalar (outDir içinde, Windows satır sonlarıyla):
#   fmm_data_lg.txt        veri (ilk satır değişken adları; LG bunu ister)
#   fmm_K1_check_n10.lgs   tanı: K = 1, quadrature nodes = 10
#   fmm_K1_check_n50.lgs   tanı: K = 1, quadrature nodes = 50
#   fmm_K2.lgs             K = 2, sınıflama dosyası yazar
#   fmm_K1to5_range.lgs    'Class nominal 1:5' — BLRT (LLdiff) ve VLMR
#   fmm_K1to5_blocks.lgs   yedek: K = 1..5 ayrı model blokları
#
# Kullanım: önce fmm_openmx_enumeration.R'yi çalıştır (simData oluşur ya da
# fmm_enumeration_output/fits.rds kaydedilir), sonra bu dosyayı source() et.
# Kendi verin için: lgData ve factorItems'ı değiştir.
# =============================================================================

outDir    <- "latentgold_files"
lgVersion <- "6.0"   # kurulu sürümün etiketi, ör. "6.1"
nodes     <- 50L     # 10 YETERSİZ olabilir; aşağıdaki tanı adımına bak

# -----------------------------------------------------------------------------
# 1. Veri
# -----------------------------------------------------------------------------
if (!exists("simData")) {
  rds <- file.path("fmm_enumeration_output", "fits.rds")
  if (!file.exists(rds))
    stop("simData yok: önce fmm_openmx_enumeration.R'yi çalıştır.")
  saved     <- readRDS(rds)
  simData   <- saved$data
  trueClass <- saved$trueClass
}

# Faktör -> madde eşlemesi (kendi verin için burayı değiştir)
factorItems <- list(F1 = paste0("y", 1:5), F2 = paste0("y", 6:10))
items <- unlist(factorItems, use.names = FALSE)

lgData <- data.frame(id = seq_len(nrow(simData)), simData[, items])
if (exists("trueClass")) lgData$trueclass <- trueClass
keepVars <- intersect(c("id", "trueclass"), names(lgData))

dir.create(outDir, showWarnings = FALSE)
eol <- "\r\n"   # LG Windows programı; CRLF en güvenli seçenek
write.table(format(lgData, digits = 15, scientific = FALSE, trim = TRUE),
            file.path(outDir, "fmm_data_lg.txt"), sep = " ", na = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE, eol = eol)

chk <- read.table(file.path(outDir, "fmm_data_lg.txt"), header = TRUE)
stopifnot(nrow(chk) == nrow(simData),
          max(abs(as.matrix(chk[, items]) - as.matrix(simData[, items]))) < 1e-10)

# -----------------------------------------------------------------------------
# 2. Sözdizimi
# -----------------------------------------------------------------------------
header <- c(sprintf("//LG%s//", lgVersion), sprintf("version = %s", lgVersion),
            "infile 'fmm_data_lg.txt'", "")

# BLRT her bootstrap örnekleminde K sınıflı modeli yeniden kestirir; az
# başlangıç seti yerel maksimuma takılıp p-değerini büyütür.
blrtSets <- function(K) if (K <= 3) 32L else 64L

# 'dependent' satırını şablondaki gibi satır başına 5 değişkene böl
depLines <- function(v) {
  ch  <- split(paste(v, "continuous"), ceiling(seq_along(v) / 5))
  txt <- vapply(ch, paste, "", collapse = ", ")
  end <- c(rep(",", length(txt) - 1), ";")
  paste0(c("   dependent ", rep("      ", length(txt) - 1)), txt, end)
}

block <- function(title, classes, lldiff = NULL, outfile = NULL, append = NULL,
                  sets = 200L, iterations = 100L, mcSets = 32L, nodesK = nodes) {
  facs <- names(factorItems)
  c("model",
    sprintf("title '%s';", title),
    "// FMM: class-invariant loadings, intercepts, residual variances, factor covariance;",
    "// factor variances fixed at 1; factor means 0 in class 1, free in other classes;",
    "// Bayes constants 0 = plain ML.",
    "options",
    "   maxthreads=all;",
    "   algorithm",
    "      tolerance=1e-008 emtolerance=0.01 emiterations=250 nriterations=100;",
    "   startvalues",
    sprintf("      seed=0 sets=%d tolerance=1e-005 iterations=%d;", sets, iterations),
    "   bayes",
    "      categorical=0 variances=0 latent=0 poisson=0;",
    "   montecarlo",
    if (is.null(lldiff)) sprintf("      seed=0 sets=%d replicates=500 tolerance=1e-008;", mcSets)
    else c(sprintf("      seed=0 sets=%d replicates=500 tolerance=1e-008", mcSets),
           sprintf("      %s;", lldiff)),
    sprintf("   quadrature  nodes=%d;", nodesK),
    "   missing  includeall;",
    "   output",
    "      parameters=first standarderrors profile probmeans=posterior",
    if (is.null(append)) "      estimatedvalues=model;"
    else c("      estimatedvalues=model", sprintf("      append = '%s';", append)),
    if (!is.null(outfile))
      c(sprintf("   outfile '%s'", outfile),
        sprintf("      classification=posterior keep %s;", paste(keepVars, collapse = ", "))),
    "variables",
    depLines(items),
    "   latent",
    sprintf("      Class nominal %s,", classes),
    sprintf("      %s continuous%s", facs, c(rep(",", length(facs) - 1), ";")),
    "equations",
    "   Class <- 1;",
    sprintf("   (1) %s;", facs),
    "   F1 <-> F2;",
    sprintf("   %s <- (%s) 1 | Class;", facs, letters[seq_along(facs)]),
    unlist(lapply(facs, function(f) sprintf("   %s <- 1 + %s;", factorItems[[f]], f))),
    sprintf("   %s;", items),
    sprintf("   %s[1] = 0;", letters[seq_along(facs)]),
    "end model", "")
}

files <- list(
  fmm_K1_check_n10.lgs = c(header, block("FMM K1 nodes10", "1", nodesK = 10L)),
  fmm_K1_check_n50.lgs = c(header, block("FMM K1 nodes50", "1", nodesK = 50L)),
  fmm_K2.lgs = c(header, block("FMM K2", "2", outfile = "fmm_K2_classification.txt")),
  fmm_K1to5_range.lgs = c(header, block("FMM K1-K5", "1:5", lldiff = "LLdiff",
                                        append = "fmm_lg_summary.csv",
                                        mcSets = 64L)),
  fmm_K1to5_blocks.lgs = c(header, unlist(lapply(1:5, function(K)
    block(sprintf("FMM K%d", K), as.character(K),
          lldiff = if (K > 1) sprintf("LLdiff='FMM K%d'", K - 1),
          outfile = sprintf("fmm_K%d_classification.txt", K),
          append = "fmm_lg_summary.csv",
          sets = if (K <= 2) 200L else 500L,
          iterations = if (K <= 2) 100L else 250L,
          mcSets = blrtSets(K)))))
)

for (f in names(files)) {
  txt <- files[[f]]
  stopifnot(all(utf8ToInt(paste(txt, collapse = "")) < 128))   # yalnızca ASCII
  con <- file(file.path(outDir, f), "wb")
  writeLines(txt, con, sep = eol)
  close(con)
}
cat("Yazıldı:", normalizePath(outDir), "\n ", paste(c("fmm_data_lg.txt", names(files)),
                                                   collapse = "\n  "), "\n")
