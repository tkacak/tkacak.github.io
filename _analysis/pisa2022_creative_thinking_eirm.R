####################################################################
## EXPLANATORY IRT / GLMM
## PISA 2022 CREATIVE THINKING
##
## MODEL MERDIVENI
##   Model 0 : Null model
##   Model 1 : Model 0 + madde duzeyi yordayici (Complexity_c)
##   Model 2 : Model 1 + ogrenci duzeyi yordayicilar
##   Model 3 : Model 2 + moderasyon (Complexity_c x ogrenci duzeyi)
##
## Dort model de AYNI rastgele etki yapisini kullanir:
##   (1 | ID) + (1 | Items)
## Bu nedenle ic ice gecmislerdir (nested) ve LRT gecerlidir.
##
## 10 READ plausible value Rubin's rules ile pool edilir.
##
## ICINDEKILER
##   00. Paketler, ayarlar, cikti klasorleri
##   01. Verilerin okunmasi
##   02. Ogrenci verisi
##   03. Creative Thinking bilissel verisi
##   04. Madde yanit sayilari
##   05. Ikili madde kodlama
##   06. ATOS complexity degiskeni
##   07. Verilerin birlestirilmesi ve EO tanimi
##   08. Ornek agirliklarinin normalizasyonu
##   09. Long veri formati
##   10. Faktor kodlamalari
##   11. Listwise deletion ve veri kontrolleri
##   12. Person duzeyinde standardizasyon (PV ve Homework)
##   13. Ulke ve EO bazinda ornek buyuklukleri
##   14. Model altyapisi (yardimci fonksiyonlar)
##   15. Rubin's rules altyapisi
##   16. Pooled LRT (D2) altyapisi
##   17. Model merdiveni fonksiyonu (Model 0-3)
##   18. Betimsel Rasch modeli (madde guclukleri)
##   19. Tum ornek icin Model 0-3
##   20. Pooled sabit etkilerin raporlanmasi
##   21. Varyans bilesenleri ve aciklanan varyans
##   22. Complexity etkisinin modellere gore ozeti
##   23. Pooled marjinal etki grafikleri
##   24. Model 0-1 icin HTML tablo
##   25. Betimsel agirlikli korelasyonlar
##   26. Ulke bazinda pooled korelasyon
##   27. EO grup analizleri
##   28. Forest plot
##   29. Yakinsama ve singularity kontrolleri
##   30. Calisma ciktilarinin kaydedilmesi
##   31. Session info
##
## ------------------------------------------------------------------
## YONTEMSEL NOTLAR
##
## 1) AGIRLIKLAR
##    Modeller AGIRLIKSIZ tahmin edilir.
##    Nedeni: glmer'e tam sayi olmayan ornekleme agirligi verildiginde
##    elde edilen sey bir pseudo-likelihood'dur. Bu durumda AIC, BIC,
##    logLik ve LRT degerleri gecerli likelihood tabanli olcutler
##    olmaktan cikar. Bu betikte model karsilastirmasi AIC/BIC/LL/LRT
##    uzerine kurulu oldugu icin modeller agirliksiz tahmin edilmistir.
##    Normalize agirliklar yine de hesaplanir ve BETIMSEL istatistikler
##    (agirlikli korelasyonlar) icin kullanilir.
##
## 2) PLAUSIBLE VALUE'LAR
##    Model 0 ve Model 1 PV icermez; tek kez tahmin edilir.
##    Model 2 ve Model 3 PVREAD icerir; 10 kez tahmin edilip
##    Rubin's rules ile pool edilir.
##    Listwise deletion PV'lerden ONCE yapildigi ve PV'lerde eksik
##    veri olmadigi icin dort modelin N'i ayni olur; betik bunu
##    ayrica kontrol eder.
##
## 3) AIC / BIC / logLik
##    Model 2 ve Model 3 icin raporlanan degerler 10 modelin
##    ARITMETIK ORTALAMASIDIR. Bunlar Rubin-pooled degildir;
##    likelihood tabanli olcutlerin Rubin kurallari ile birlestirilmesi
##    tanimli degildir. Tablolarda bu durum acikca belirtilir.
##
## 4) LRT
##    Model 0 - Model 1 karsilastirmasi tek modelden yapildigi icin
##    klasik LRT'dir (anova).
##    Model 1 - Model 2 ve Model 2 - Model 3 karsilastirmalari
##    PV'ler uzerinden D2 yontemi ile pool edilir
##    (Li, Meng, Raghunathan & Rubin, 1991).
##
## 5) R-KARE
##    Nakagawa & Schielzeth (2013) marjinal ve kosullu R2 degerleri
##    performance::r2_nakagawa ile hesaplanir. Model 2 ve Model 3 icin
##    10 PV modelinin ortalamasi raporlanir.
####################################################################


####################################################################
## 00. PAKETLER, AYARLAR VE CIKTI KLASORLERI
####################################################################

library(haven)
library(tidyr)
library(dplyr)
library(lme4)
library(purrr)
library(tibble)
library(broom.mixed)   ## mice::pool'un merMod'u tanimasi icin gerekli
library(performance)
library(mice)
library(mitml)
library(ggplot2)

## Istege bagli paketler --------------------------------------------
HAS_EIRM  <- requireNamespace("eirm",  quietly = TRUE)
HAS_SJP   <- requireNamespace("sjPlot", quietly = TRUE)
HAS_EXPSS <- requireNamespace("expss", quietly = TRUE)

if (!HAS_EIRM) {
  message(
    "eirm paketi bulunamadi. Modeller dogrudan lme4::glmer ile ",
    "tahmin edilecek. Sonuclar aynidir; yalnizca eirm'e ozgu ",
    "print/plot yardimcilari devre disi kalir."
  )
}

## Cikti klasorleri -------------------------------------------------
OUT_MAIN <- "model_outputs_revised"
OUT_EO   <- file.path(OUT_MAIN, "EO")

dir.create(OUT_MAIN, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_EO,   showWarnings = FALSE, recursive = TRUE)

## Kisa yol yardimcisi
out_main <- function(...) file.path(OUT_MAIN, ...)
out_eo   <- function(...) file.path(OUT_EO, ...)

## Etiketleri guvenli sekilde dusuren yardimci ----------------------
drop_labs <- function(x) {
  if (HAS_EXPSS) {
    expss::drop_var_labs(x)
  } else {
    as.data.frame(lapply(x, function(col) {
      attr(col, "label")  <- NULL
      attr(col, "labels") <- NULL
      if (inherits(col, "haven_labelled")) col <- haven::zap_labels(col)
      col
    }), stringsAsFactors = FALSE)
  }
}


####################################################################
## 01. VERILERIN OKUNMASI
####################################################################

stu <- read_sav("STU_QQQ.SAV")
crt <- read_sav("CRT_COG.SAV")


####################################################################
## 02. OGRENCI VERISI
##     10 READ plausible value korunur
####################################################################

pv_names <- paste0("PV", 1:10, "READ")
M_PV     <- length(pv_names)

student <- stu %>%
  dplyr::select(
    CNTSTUID,
    W_FSTUWT,
    LANGTEST_COG,
    ST004D01T,
    ST022Q01TA,
    ST296Q02JA,
    PV1READ:PV10READ,
    ESCS
  ) %>%
  filter(LANGTEST_COG == 313) %>%
  dplyr::select(-LANGTEST_COG) %>%
  drop_labs() %>%
  mutate(
    CNTSTUID   = as.character(CNTSTUID),
    W_FSTUWT   = as.numeric(W_FSTUWT),
    ST004D01T  = as.numeric(ST004D01T),
    ST022Q01TA = as.numeric(ST022Q01TA),
    ST296Q02JA = as.numeric(ST296Q02JA),
    ESCS       = as.numeric(ESCS),
    across(all_of(pv_names), as.numeric)
  )

## PISA'da gecersiz kodlar haven tarafindan NA'ya cevrilir.
## Yine de beklenmeyen kodlara karsi acik bir guvenlik kontrolu:
student <- student %>%
  mutate(
    ST004D01T  = ifelse(ST004D01T  %in% c(1, 2), ST004D01T,  NA_real_),
    ST022Q01TA = ifelse(ST022Q01TA %in% c(1, 2), ST022Q01TA, NA_real_)
  )


####################################################################
## 03. CREATIVE THINKING BILISSEL VERISI
####################################################################

creative <- crt %>%
  dplyr::select(
    CNT,
    CNTSTUID,
    LANGTEST_COG,
    starts_with("DT")
  ) %>%
  filter(LANGTEST_COG == 313) %>%
  dplyr::select(-LANGTEST_COG)

data_crt <- creative %>%
  dplyr::select(
    CNT,
    CNTSTUID,

    DT200Q01C2,  # ScienceFairPoster - Q01
    DT200Q02C2,  # ScienceFairPoster - Q02

    DT240Q01C2,  # SpaceComic - Q01
    DT240Q02C,   # SpaceComic - Q02

    DT690Q01C,   # SaveTheRiver - Q01
    DT690Q02C2,  # SaveTheRiver - Q02

    DT300Q01C2,  # IllustrationTitles - Q01
    DT300Q02C,   # IllustrationTitles - Q02

    DT400Q01C,   # SaveTheBees - Q01
    DT400Q02C2,  # SaveTheBees - Q02
    DT400Q03C2,  # SaveTheBees - Q03

    DT500Q01C,   # LibraryAccessibility - Q01
    DT500Q02C2,  # LibraryAccessibility - Q02

    DT570Q01C,   # RobotStory - Q01
    DT570Q02C2,  # RobotStory - Q02
    DT570Q03C2,  # RobotStory - Q03

    DT370Q01C2,  # 2983 - Q01
    DT630Q01C2   # Carpooling - Q01
  ) %>%
  drop_labs() %>%
  mutate(
    CNT      = as.character(CNT),
    CNTSTUID = as.character(CNTSTUID),
    across(-c(CNT, CNTSTUID), as.numeric)
  )

item_cols <- setdiff(colnames(data_crt), c("CNT", "CNTSTUID"))

## CNTSTUID birlestirme anahtari olarak tekil olmali
if (anyDuplicated(data_crt$CNTSTUID) > 0) {
  stop(
    "CNTSTUID bilissel veride tekil degil. ",
    "Birlestirme CNT + CNTSTUID uzerinden yapilmalidir."
  )
}


####################################################################
## 04. HER MADDEYE YANIT VEREN BIREY SAYISI
##
## Ikili kodlamadan ONCE hesaplanir.
## Uygulanmayan maddeler NA olarak korunur.
####################################################################

item_response_counts <- data_crt %>%
  summarise(
    across(
      all_of(item_cols),
      list(
        n_responded = ~ sum(!is.na(.x)),
        n_missing   = ~ sum(is.na(.x)),
        n_correct   = ~ sum(.x >= 1, na.rm = TRUE)
      ),
      .names = "{.col}__{.fn}"
    )
  ) %>%
  pivot_longer(
    cols      = everything(),
    names_to  = c("Item", ".value"),
    names_sep = "__"
  ) %>%
  mutate(
    n_total   = nrow(data_crt),
    p_admin   = ifelse(n_total > 0, n_responded / n_total, NA_real_),
    p_correct = ifelse(n_responded > 0, n_correct / n_responded, NA_real_),
    p_admin   = round(p_admin, 3),
    p_correct = round(p_correct, 3)
  ) %>%
  arrange(desc(n_responded))

print(as.data.frame(item_response_counts))

write.csv(
  item_response_counts,
  out_main("item_response_counts.csv"),
  row.names = FALSE
)


####################################################################
## 04.1. ULKE X MADDE YANIT SAYISI
####################################################################

item_response_counts_by_country <- data_crt %>%
  group_by(CNT) %>%
  summarise(across(all_of(item_cols), ~ sum(!is.na(.x))), .groups = "drop")

print(as.data.frame(item_response_counts_by_country))

write.csv(
  item_response_counts_by_country,
  out_main("item_response_counts_by_country.csv"),
  row.names = FALSE
)


####################################################################
## 05. IKILI MADDE KODLAMA
##
## 0 puan       -> 0
## 1 veya uzeri -> 1  (kismi puan dogru kabul edilir)
## NA           -> NA  (uygulanmayan madde yanlis olarak kodlanmaz)
####################################################################

data_crt <- data_crt %>%
  mutate(
    across(
      all_of(item_cols),
      ~ case_when(
        is.na(.x) ~ NA_real_,
        .x >= 1   ~ 1,
        TRUE      ~ 0
      )
    )
  )


####################################################################
## 06. ATOS COMPLEXITY
####################################################################

atos_complexity <- c(
  "DT200Q01C2" = 8.07,
  "DT200Q02C2" = 8.27,
  "DT240Q01C2" = 7.45,
  "DT240Q02C"  = 5.37,
  "DT690Q01C"  = 7.97,
  "DT690Q02C2" = 7.88,
  "DT300Q01C2" = 7.21,
  "DT300Q02C"  = 5.52,
  "DT400Q01C"  = 8.00,
  "DT400Q02C2" = 7.80,
  "DT400Q03C2" = 6.83,
  "DT500Q01C"  = 7.88,
  "DT500Q02C2" = 8.27,
  "DT570Q01C"  = 7.34,
  "DT570Q02C2" = 7.20,
  "DT570Q03C2" = 8.16,
  "DT370Q01C2" = 5.71,
  "DT630Q01C2" = 9.14
)

complexity_df <- tibble(
  Item       = item_cols,
  Complexity = as.numeric(atos_complexity[item_cols])
)

if (anyNA(complexity_df$Complexity)) {
  stop(
    "ATOS complexity degeri bulunamayan maddeler: ",
    paste(complexity_df$Item[is.na(complexity_df$Complexity)], collapse = ", ")
  )
}

complexity_mean <- mean(complexity_df$Complexity, na.rm = TRUE)

complexity_df <- complexity_df %>%
  mutate(Complexity_c = Complexity - complexity_mean)

cat("Complexity ortalamasi (madde havuzu):",
    round(complexity_mean, 3), "\n")


####################################################################
## 07. VERILERIN BIRLESTIRILMESI VE EO TANIMI
####################################################################

## Ingilizce'nin resmi dil oldugu sistemler
english_official <- c(
  "AUS",
  "CAN",
  "HKG",
  "JAM",
  "MLT",
  "PHL",
  "SGP"
)

merged_data <- data_crt %>%
  left_join(student, by = "CNTSTUID") %>%
  mutate(
    CNT = factor(CNT),
    EO  = ifelse(as.character(CNT) %in% english_official, 1, 0)
  )

## Birlestirme satir sayisini degistirmemeli
if (nrow(merged_data) != nrow(data_crt)) {
  stop(
    "left_join satir sayisini degistirdi: ogrenci verisinde ",
    "tekrarli CNTSTUID var."
  )
}

## EO siniflamasinin gozle kontrolu:
## listede olmayan ulkeler burada gorunur.
cat("\n--- EO siniflamasi ---\n")
cat("Formal (EO = 1):",
    paste(sort(unique(as.character(merged_data$CNT[merged_data$EO == 1]))),
          collapse = ", "), "\n")
cat("Non-Formal (EO = 0):",
    paste(sort(unique(as.character(merged_data$CNT[merged_data$EO == 0]))),
          collapse = ", "), "\n")
cat("Listede olup veride bulunmayan kodlar:",
    paste(setdiff(english_official, unique(as.character(merged_data$CNT))),
          collapse = ", "), "\n\n")


####################################################################
## 08. ORNEK AGIRLIKLARININ NORMALIZASYONU
##
## Normalizasyon kisi duzeyinde yapilir (merged_data satirlari
## zaten kisi duzeyindedir).
##
##   WEIGHT        : ulke icinde ortalama = 1
##   WEIGHT_GLOBAL : tum ornekte ortalama = 1
##   WEIGHT_SENATE : her ulkenin toplam katkisi = 1000
##
## NOT: Bu agirliklar modellerde KULLANILMAZ (bkz. bas taraftaki
## yontemsel not 1). Yalnizca betimsel agirlikli korelasyonlarda
## kullanilirlar.
####################################################################

merged_data <- merged_data %>%
  group_by(CNT) %>%
  mutate(
    WEIGHT        = W_FSTUWT / mean(W_FSTUWT, na.rm = TRUE),
    WEIGHT_SENATE = 1000 * W_FSTUWT / sum(W_FSTUWT, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    WEIGHT_GLOBAL = W_FSTUWT / mean(W_FSTUWT, na.rm = TRUE)
  )


####################################################################
## 09. LONG VERI FORMATI
##
## PV'ler ve Homework burada standardize EDILMEZ.
## Standardizasyon listwise deletion sonrasinda, kisi duzeyinde
## bir kez yapilir (bkz. bolum 12).
####################################################################

id_cols <- c(
  "CNTSTUID",
  "CNT",
  "W_FSTUWT",
  "WEIGHT",
  "WEIGHT_GLOBAL",
  "WEIGHT_SENATE",
  "ST004D01T",
  "ST022Q01TA",
  "ST296Q02JA",
  "EO",
  "ESCS",
  pv_names
)

data_long <- merged_data %>%
  pivot_longer(
    cols      = all_of(item_cols),
    names_to  = "Item",
    values_to = "Response"
  ) %>%
  dplyr::select(all_of(id_cols), Item, Response) %>%
  left_join(complexity_df, by = "Item") %>%
  rename(
    ID           = CNTSTUID,
    Gender       = ST004D01T,
    Language     = ST022Q01TA,
    Homework_raw = ST296Q02JA,
    Items        = Item
  ) %>%
  mutate(
    ID           = factor(ID),
    CNT          = factor(CNT),
    Items        = factor(Items),
    Response     = as.numeric(Response),
    Homework_raw = as.numeric(Homework_raw),
    ESCS         = as.numeric(ESCS),
    Complexity   = as.numeric(Complexity),
    Complexity_c = as.numeric(Complexity_c),
    across(all_of(pv_names), as.numeric)
  )


####################################################################
## 10. FAKTOR KODLAMALARI
####################################################################

data_long <- data_long %>%
  mutate(
    Gender = factor(
      as.numeric(Gender),
      levels = c(1, 2),
      labels = c("Female", "Male")
    ),
    Language = factor(
      as.numeric(Language),
      levels = c(1, 2),
      labels = c("Test", "Other")
    ),
    EO = factor(
      as.numeric(EO),
      levels = c(0, 1),
      labels = c("Non-Formal", "Formal")
    )
  ) %>%
  mutate(
    Gender   = relevel(Gender,   ref = "Female"),
    Language = relevel(Language, ref = "Test"),
    EO       = relevel(EO,       ref = "Non-Formal")
  )


####################################################################
## 11. LISTWISE DELETION VE VERI KONTROLLERI
##
## PV'ler icin listwise deletion yapilmaz: PISA'da PV'ler eksiksizdir
## ve zaten tum modellerde ayni kisiler yer almalidir.
####################################################################

n_before <- n_distinct(data_long$ID)

data_long <- data_long %>%
  filter(
    !is.na(ID),
    !is.na(CNT),
    !is.na(Items),
    !is.na(Response),
    !is.na(Homework_raw),
    !is.na(Language),
    !is.na(ESCS),
    !is.na(Gender),
    !is.na(EO),
    !is.na(WEIGHT),
    !is.na(Complexity_c)
  ) %>%
  droplevels()

## PV'lerde eksik veri kalmamali: aksi halde Model 2-3'un N'i
## Model 0-1'den farkli olur ve LRT gecersizlesir.
pv_missing <- vapply(
  data_long[pv_names],
  function(x) sum(is.na(x)),
  numeric(1)
)

if (any(pv_missing > 0)) {
  stop(
    "PV kolonlarinda eksik veri var; Model 0-1 ile Model 2-3'un N'i ",
    "farklilasir. Eksik sayilari: ",
    paste(names(pv_missing), pv_missing, sep = "=", collapse = ", ")
  )
}

cat("Listwise deletion oncesi birey:", n_before, "\n")
cat("Listwise deletion sonrasi birey:", n_distinct(data_long$ID), "\n")
cat("Toplam ulke:", n_distinct(data_long$CNT), "\n")
cat("Toplam yanit (satir):", nrow(data_long), "\n")

## Matrix sinifinda kolon kalmamali (scale() tuzagi)
matrix_columns <- names(data_long)[vapply(data_long, is.matrix, logical(1))]

cat("Matrix sinifindaki sutunlar:",
    if (length(matrix_columns) == 0) "YOK" else paste(matrix_columns, collapse = ", "),
    "\n")

if (length(matrix_columns) > 0) {
  stop("Matrix sinifinda sutunlar bulundu: ",
       paste(matrix_columns, collapse = ", "))
}


####################################################################
## 12. PERSON DUZEYINDE STANDARDIZASYON
##
## PV'ler ve Homework kisi duzeyinde standardize edilir.
## Satir (yanit) duzeyinde standardizasyon, cok madde yanitlayan
## ogrencilere fazladan agirlik verirdi.
##
## Standardizasyon TUM analiz orneklemi uzerinden yapilir. Boylece
## EO gruplarinin katsayilari ortak bir olcekte kalir ve
## karsilastirilabilir olur.
####################################################################

person_level <- data_long %>%
  distinct(ID, .keep_all = TRUE)

standardize_with <- function(x, mu, s) {
  if (is.na(s) || s == 0) {
    stop("Standardizasyon icin standart sapma sifir veya NA.")
  }
  as.numeric((x - mu) / s)
}

## Homework -------------------------------------------------------
hw_mu <- mean(person_level$Homework_raw, na.rm = TRUE)
hw_sd <- sd(person_level$Homework_raw,   na.rm = TRUE)

data_long$Homework <- standardize_with(data_long$Homework_raw, hw_mu, hw_sd)

cat("Homework standardizasyonu: mean =", round(hw_mu, 3),
    "sd =", round(hw_sd, 3), "\n")

## Plausible values ------------------------------------------------
pv_z_names <- paste0(pv_names, "_z")

pv_moments <- tibble(
  PV   = pv_names,
  mean = vapply(person_level[pv_names], mean, numeric(1), na.rm = TRUE),
  sd   = vapply(person_level[pv_names], sd,   numeric(1), na.rm = TRUE)
)

print(as.data.frame(pv_moments))

for (k in seq_along(pv_names)) {
  data_long[[pv_z_names[k]]] <- standardize_with(
    data_long[[pv_names[k]]],
    pv_moments$mean[k],
    pv_moments$sd[k]
  )
}

write.csv(pv_moments, out_main("pv_standardization_moments.csv"),
          row.names = FALSE)

## Son kontrol: hicbir kolon matrix olmamali
stopifnot(!any(vapply(data_long, is.matrix, logical(1))))

str(data_long[c("Homework", "ESCS", pv_z_names[1:2])])

saveRDS(data_long, out_main("data_long.RDS"))


####################################################################
## 13. ULKE VE EO BAZINDA ORNEK BUYUKLUKLERI
####################################################################

country_counts <- data_long %>%
  group_by(CNT) %>%
  summarise(
    n_students   = n_distinct(ID),
    n_responses  = n(),
    n_items_mean = round(n() / n_distinct(ID), 2),
    .groups = "drop"
  ) %>%
  mutate(pct_students = round(100 * n_students / sum(n_students), 2)) %>%
  arrange(desc(n_students))

print(as.data.frame(country_counts))

cat("TOPLAM birey:", sum(country_counts$n_students),
    "| TOPLAM ulke:", nrow(country_counts), "\n")

write.csv(country_counts, out_main("country_counts.csv"), row.names = FALSE)

country_counts_by_EO <- data_long %>%
  group_by(EO, CNT) %>%
  summarise(n_students = n_distinct(ID), .groups = "drop") %>%
  arrange(EO, desc(n_students))

print(as.data.frame(country_counts_by_EO))

write.csv(country_counts_by_EO, out_main("country_counts_by_EO.csv"),
          row.names = FALSE)

raw_country_counts <- merged_data %>%
  group_by(CNT) %>%
  summarise(n_students_raw = n_distinct(CNTSTUID), .groups = "drop") %>%
  arrange(desc(n_students_raw))

print(as.data.frame(raw_country_counts))

write.csv(raw_country_counts, out_main("raw_country_counts.csv"),
          row.names = FALSE)


####################################################################
## 14. MODEL ALTYAPISI
####################################################################

ctrl_glmer <- glmerControl(
  optimizer   = "bobyqa",
  calc.derivs = FALSE,
  optCtrl     = list(maxfun = 200000)
)

## Model formulleri ------------------------------------------------
##
## Dort modelin rastgele etki yapisi AYNIDIR; boylece
## Model 0 < Model 1 < Model 2 < Model 3 ic ice gecmistir.

re_terms      <- "(1 | ID) + (1 | Items)"
item_terms    <- "Complexity_c"
person_terms  <- c("PVREAD", "Language", "Homework", "ESCS", "Gender")

f_m0 <- paste("Response ~ 1 +", re_terms)

f_m1 <- paste("Response ~ 1 +", item_terms, "+", re_terms)

f_m2 <- paste(
  "Response ~ 1 +", item_terms, "+",
  paste(person_terms, collapse = " + "), "+", re_terms
)

f_m3 <- paste0(
  "Response ~ 1 + ", item_terms,
  " * (", paste(person_terms, collapse = " + "), ") + ", re_terms
)

## Betimsel Rasch modeli: madde guclukleri icin
## (model merdivenine DAHIL DEGILDIR)
f_rasch <- "Response ~ -1 + Items + (1 | ID)"

model_labels <- c(
  m0 = "Model 0: Null",
  m1 = "Model 1: + item-level (Complexity)",
  m2 = "Model 2: + student-level",
  m3 = "Model 3: + moderation"
)

## eirm varsa eirm, yoksa dogrudan glmer --------------------------
fit_model <- function(formula_text, data, control = ctrl_glmer) {
  if (HAS_EIRM) {
    eirm::eirm(
      formula = formula_text,
      data    = as.data.frame(data),
      control = control
    )
  } else {
    lme4::glmer(
      stats::as.formula(formula_text),
      data    = as.data.frame(data),
      family  = binomial(link = "logit"),
      control = control
    )
  }
}

## eirm nesnesinden glmer nesnesini cikarir
as_glmer <- function(x) {
  if (inherits(x, "eirm")) x$model else x
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

## Yakinsama / singularity kontrolu --------------------------------
check_glmer_model <- function(model, model_name) {

  g   <- as_glmer(model)
  msg <- g@optinfo$conv$lme4$messages

  tibble(
    Model               = model_name,
    Singular            = lme4::isSingular(g, tol = 1e-4),
    Convergence_message = if (is.null(msg)) "No convergence message"
                          else paste(msg, collapse = " | "),
    N_observations      = stats::nobs(g),
    LogLik              = as.numeric(stats::logLik(g)),
    AIC                 = stats::AIC(g),
    BIC                 = stats::BIC(g)
  )
}

## Tek modelin uyum istatistikleri ---------------------------------
r2_safe <- function(model) {
  out <- tryCatch(
    suppressWarnings(performance::r2_nakagawa(model)),
    error = function(e) NULL
  )

  if (is.null(out)) {
    return(c(R2_marginal = NA_real_, R2_conditional = NA_real_))
  }

  c(
    R2_marginal    = unname(as.numeric(out$R2_marginal    %||% NA_real_)),
    R2_conditional = unname(as.numeric(out$R2_conditional %||% NA_real_))
  )
}

fit_stats <- function(model) {
  g  <- as_glmer(model)
  r2 <- r2_safe(g)

  tibble(
    n_obs          = stats::nobs(g),
    n_fixed        = length(lme4::fixef(g)),
    npar           = attr(stats::logLik(g), "df"),
    AIC            = stats::AIC(g),
    BIC            = stats::BIC(g),
    logLik         = as.numeric(stats::logLik(g)),
    R2_marginal    = unname(r2["R2_marginal"]),
    R2_conditional = unname(r2["R2_conditional"])
  )
}

## PV'ler uzerinde model tahmini -----------------------------------
##
## Her PV icin PVREAD kolonu degistirilerek model yeniden tahmin
## edilir. Bellek kullanimini dusurmek icin PV kolonlari dislanir.
fit_across_pvs <- function(formula_text,
                           dat,
                           pv_cols = pv_z_names,
                           control = ctrl_glmer,
                           label   = "") {

  keep <- setdiff(names(dat), c(pv_names, pv_z_names))

  lapply(seq_along(pv_cols), function(k) {
    d_k         <- dat[keep]
    d_k$PVREAD  <- dat[[pv_cols[k]]]

    message(sprintf("  [%s] PV %d / %d tahmin ediliyor...",
                    label, k, length(pv_cols)))

    fit <- as_glmer(fit_model(formula_text, d_k, control))
    rm(d_k)
    fit
  })
}

## PV modellerinin ortalama uyum istatistikleri --------------------
##
## DIKKAT: Bu degerler Rubin-pooled DEGILDIR; 10 modelin aritmetik
## ortalamasidir.
fit_stats_across <- function(fit_list) {
  by_pv <- purrr::map_dfr(seq_along(fit_list), function(i) {
    fit_stats(fit_list[[i]]) %>% mutate(PV = i, .before = 1)
  })

  ## DIKKAT: summarise() ifadeleri sirayla degerlendirir. Standart
  ## sapmalar ortalamalardan ONCE hesaplanmalidir; aksi halde sd()
  ## tek bir ortalama degerin uzerinde calisir ve NA doner.
  avg <- by_pv %>%
    summarise(
      n_obs          = dplyr::first(n_obs),
      n_fixed        = dplyr::first(n_fixed),
      npar           = dplyr::first(npar),
      AIC_sd         = sd(AIC,    na.rm = TRUE),
      BIC_sd         = sd(BIC,    na.rm = TRUE),
      logLik_sd      = sd(logLik, na.rm = TRUE),
      AIC            = mean(AIC,    na.rm = TRUE),
      BIC            = mean(BIC,    na.rm = TRUE),
      logLik         = mean(logLik, na.rm = TRUE),
      R2_marginal    = mean(R2_marginal,    na.rm = TRUE),
      R2_conditional = mean(R2_conditional, na.rm = TRUE)
    ) %>%
    dplyr::select(
      n_obs, n_fixed, npar,
      AIC, BIC, logLik,
      R2_marginal, R2_conditional,
      AIC_sd, BIC_sd, logLik_sd
    )

  list(by_pv = by_pv, avg = avg)
}


####################################################################
## 15. RUBIN'S RULES ALTYAPISI
##
## Sabit etkiler Rubin (1987) kurallari ile birlestirilir:
##   qbar   = ortalama tahmin
##   ubar   = ortalama within varyans
##   b      = between varyans
##   t      = ubar + (1 + 1/m) * b
##   riv    = (1 + 1/m) * b / ubar
##   lambda = (1 + 1/m) * b / t
##   df     = (m - 1) / lambda^2            (dfcom = Inf)
##   fmi    = (riv + 2 / (df + 3)) / (riv + 1)
##
## Bu uygulama mitml::testEstimates() ile makine hassasiyetinde
## ayni sonucu verir; mice::pool()'a bagimli olmadigi icin paket
## surumlerinden etkilenmez.
####################################################################

rubin_pool_fixed <- function(fit_list, conf = 0.95, dfcom = Inf) {

  m   <- length(fit_list)
  if (m < 2) stop("Rubin birlestirmesi icin en az 2 model gerekir.")

  nms <- names(lme4::fixef(fit_list[[1]]))

  ## Tum modellerde ayni parametreler olmali
  for (i in seq_along(fit_list)) {
    if (!identical(names(lme4::fixef(fit_list[[i]])), nms)) {
      stop("PV modelleri farkli sabit etki parametrelerine sahip (model ", i, ").")
    }
  }

  Q <- t(vapply(fit_list, lme4::fixef, numeric(length(nms))))
  U <- t(vapply(
    fit_list,
    function(f) diag(as.matrix(stats::vcov(f))),
    numeric(length(nms))
  ))

  qbar   <- colMeans(Q)
  ubar   <- colMeans(U)
  b      <- apply(Q, 2, stats::var)
  tvar   <- ubar + (1 + 1 / m) * b
  riv    <- (1 + 1 / m) * b / ubar
  lambda <- (1 + 1 / m) * b / tvar

  df_old <- (m - 1) / lambda^2

  if (is.finite(dfcom)) {
    df_obs <- ((dfcom + 1) / (dfcom + 3)) * dfcom * (1 - lambda)
    dfv    <- df_old * df_obs / (df_old + df_obs)
  } else {
    dfv <- df_old
  }

  fmi  <- (riv + 2 / (dfv + 3)) / (riv + 1)
  se   <- sqrt(tvar)
  stat <- qbar / se
  crit <- stats::qt(1 - (1 - conf) / 2, df = dfv)

  tibble(
    term      = nms,
    estimate  = qbar,
    std.error = se,
    statistic = stat,
    df        = dfv,
    p.value   = 2 * stats::pt(-abs(stat), df = dfv),
    conf.low  = qbar - crit * se,
    conf.high = qbar + crit * se,
    ubar      = ubar,
    b         = b,
    t         = tvar,
    riv       = riv,
    lambda    = lambda,
    fmi       = fmi,
    M         = m
  )
}

## Sabit etkilerin TAM pooled kovaryans matrisi
## (marjinal etki grafiklerinde lineer bilesimlerin SE'si icin)
rubin_total_vcov <- function(fit_list) {
  m      <- length(fit_list)
  B_list <- lapply(fit_list, lme4::fixef)
  V_list <- lapply(fit_list, function(f) as.matrix(stats::vcov(f)))

  b_bar <- Reduce(`+`, B_list) / m
  W     <- Reduce(`+`, V_list) / m
  Bvar  <- Reduce(`+`, lapply(B_list, function(b) tcrossprod(b - b_bar))) / (m - 1)

  list(b = b_bar, T = W + (1 + 1 / m) * Bvar)
}

## Rastgele etki varyanslarinin PV'ler uzerindeki ortalamasi
ran_pars_across <- function(fit_list) {
  purrr::map_dfr(seq_along(fit_list), function(i) {
    as.data.frame(lme4::VarCorr(fit_list[[i]])) %>%
      dplyr::select(grp, var1, var2, vcov, sdcor) %>%
      mutate(PV = i)
  }) %>%
    group_by(grp, var1, var2) %>%
    summarise(
      vcov  = mean(vcov,  na.rm = TRUE),
      sdcor = mean(sdcor, na.rm = TRUE),
      .groups = "drop"
    )
}


####################################################################
## 16. POOLED LRT (D2) ALTYAPISI
##
## Li, Meng, Raghunathan & Rubin (1991) D2 istatistigi kullanilir.
## Hesaplama mitml::testModels() ile yapilir; basarisiz olursa
## mice::D2() denenir. Her iki paket de ayni istatistigi uygular.
##
## Hicbir yol calismazsa PV bazindaki LRT'ler raporlanir ve pooled
## degerler NA birakilir. Formul elle yeniden yazilmaz.
####################################################################

## mitml / mice cikti nesnelerinden test degerlerini guvenli cikarma
extract_test_row <- function(x) {
  if (is.matrix(x)) x <- x[1, ]
  nm <- names(x)
  v  <- as.numeric(x)

  pick <- function(pattern, default_idx) {
    idx <- grep(pattern, nm)
    if (length(idx) >= 1) v[idx[1]] else
      if (!is.na(default_idx) && default_idx <= length(v)) v[default_idx] else NA_real_
  }

  c(
    statistic = pick("^F|statistic", 1L),
    df1       = pick("^df1$",        2L),
    df2       = pick("^df2$",        3L),
    p.value   = pick("^P\\(|p\\.value", 4L),
    riv       = pick("^RIV$|^riv$",  5L)
  )
}

## full  : PV basina tahmin edilmis TAM model listesi
## reduce: PV basina tahmin edilmis KISITLI model listesi
##         (PV icermeyen bir model tek nesne olarak da verilebilir)
pooled_lrt <- function(full, reduced, comparison = "") {

  if (!is.list(full)) full <- list(full)

  ## Kisitli model PV'den bagimsizsa ayni nesne tekrarlanir
  if (!is.list(reduced)) reduced <- rep(list(reduced), length(full))
  if (length(reduced) == 1 && length(full) > 1) {
    reduced <- rep(reduced, length(full))
  }

  stopifnot(length(full) == length(reduced))

  m <- length(full)

  ## PV bazinda klasik LRT ----------------------------------------
  by_pv <- purrr::map_dfr(seq_len(m), function(i) {
    ll1 <- as.numeric(stats::logLik(full[[i]]))
    ll0 <- as.numeric(stats::logLik(reduced[[i]]))
    df1 <- attr(stats::logLik(full[[i]]),    "df")
    df0 <- attr(stats::logLik(reduced[[i]]), "df")

    tibble(
      PV      = i,
      chisq   = 2 * (ll1 - ll0),
      df      = df1 - df0,
      p.value = stats::pchisq(2 * (ll1 - ll0), df1 - df0, lower.tail = FALSE)
    )
  })

  ## Tum modellerde N ayni olmali
  n_obs <- vapply(c(full, reduced), stats::nobs, numeric(1))
  if (length(unique(n_obs)) != 1L) {
    stop("LRT gecersiz: karsilastirilan modellerin N'i farkli (",
         paste(unique(n_obs), collapse = ", "), ").")
  }

  ## Tek model varsa klasik LRT yeterlidir ------------------------
  if (m == 1) {
    return(list(
      pooled = tibble(
        Comparison = comparison,
        Method     = "LRT (tek model, PV yok)",
        statistic  = by_pv$chisq[1],
        df1        = by_pv$df[1],
        df2        = NA_real_,
        p.value    = by_pv$p.value[1],
        riv        = NA_real_
      ),
      by_pv = by_pv
    ))
  }

  ## D2: once mitml, sonra mice ------------------------------------
  res    <- NULL
  method <- NA_character_

  res <- tryCatch({
    tst <- mitml::testModels(
      full, reduced,
      method = "D2",
      use    = "likelihood"
    )
    method <- "D2 (likelihood, mitml::testModels)"
    extract_test_row(tst$test)
  }, error = function(e) {
    message("mitml::testModels basarisiz: ", conditionMessage(e))
    NULL
  })

  if (is.null(res)) {
    res <- tryCatch({
      d2 <- mice::D2(
        mice::as.mira(full),
        mice::as.mira(reduced),
        use = "likelihood"
      )
      method <- "D2 (likelihood, mice::D2)"
      extract_test_row(d2$result)
    }, error = function(e) {
      message("mice::D2 basarisiz: ", conditionMessage(e))
      NULL
    })
  }

  if (is.null(res)) {
    message(
      "Pooled LRT hesaplanamadi (", comparison, "). ",
      "Yalnizca PV bazindaki LRT'ler raporlanacak."
    )
    res    <- c(statistic = NA_real_, df1 = by_pv$df[1],
                df2 = NA_real_, p.value = NA_real_, riv = NA_real_)
    method <- "Hesaplanamadi"
  }

  list(
    pooled = tibble(
      Comparison = comparison,
      Method     = method,
      statistic  = unname(res["statistic"]),
      df1        = unname(res["df1"]),
      df2        = unname(res["df2"]),
      p.value    = unname(res["p.value"]),
      riv        = unname(res["riv"])
    ),
    by_pv = by_pv
  )
}


####################################################################
## 17. MODEL MERDIVENI FONKSIYONU
##
## Model 0 - Model 3 dizisini tahmin eder, Rubin ile pool eder,
## uyum istatistiklerini ve LRT'leri hesaplar.
## Hem tum ornek hem de EO gruplari icin kullanilir.
####################################################################

run_model_ladder <- function(dat, label, out_dir) {

  cat("\n==================================================\n")
  cat("MODEL MERDIVENI:", label, "\n")
  cat("Birey:", n_distinct(dat$ID),
      "| Yanit:", nrow(dat),
      "| Ulke:", n_distinct(dat$CNT), "\n")
  cat("==================================================\n")

  dat <- droplevels(as.data.frame(dat))

  ## --- Model 0 ve Model 1: PV icermez, tek kez tahmin edilir ----
  ## Ham nesneler de saklanir: eirm kullaniliyorsa eirm'e ozgu
  ## print/plot yardimcilari icin gerekli olurlar.
  message("[", label, "] Model 0 tahmin ediliyor...")
  m0_raw <- fit_model(f_m0, dat)
  m0     <- as_glmer(m0_raw)

  message("[", label, "] Model 1 tahmin ediliyor...")
  m1_raw <- fit_model(f_m1, dat)
  m1     <- as_glmer(m1_raw)

  ## --- Model 2 ve Model 3: 10 PV ---------------------------------
  message("[", label, "] Model 2 tahmin ediliyor...")
  m2_fits <- fit_across_pvs(f_m2, dat, label = paste(label, "M2"))

  message("[", label, "] Model 3 tahmin ediliyor...")
  m3_fits <- fit_across_pvs(f_m3, dat, label = paste(label, "M3"))

  ## --- N kontrolu: LRT ve AIC/BIC icin sart ----------------------
  n_all <- c(
    stats::nobs(m0),
    stats::nobs(m1),
    vapply(m2_fits, stats::nobs, numeric(1)),
    vapply(m3_fits, stats::nobs, numeric(1))
  )

  if (length(unique(n_all)) != 1L) {
    stop(
      "Modellerin N'i farkli: ", paste(unique(n_all), collapse = ", "),
      ". AIC/BIC/LRT karsilastirmalari gecersiz olur."
    )
  }

  cat("Tum modellerde N (yanit) =", unique(n_all), "\n")

  ## --- Rubin pooled sabit etkiler --------------------------------
  pooled_m2 <- rubin_pool_fixed(m2_fits)
  pooled_m3 <- rubin_pool_fixed(m3_fits)

  ## --- Uyum istatistikleri ---------------------------------------
  s0 <- fit_stats(m0)
  s1 <- fit_stats(m1)
  s2 <- fit_stats_across(m2_fits)
  s3 <- fit_stats_across(m3_fits)

  comparison <- bind_rows(
    s0 %>% mutate(Model = unname(model_labels["m0"]),
                  PV_based = FALSE, .before = 1),
    s1 %>% mutate(Model = unname(model_labels["m1"]),
                  PV_based = FALSE, .before = 1),
    s2$avg %>% mutate(Model = unname(model_labels["m2"]),
                      PV_based = TRUE, .before = 1),
    s3$avg %>% mutate(Model = unname(model_labels["m3"]),
                      PV_based = TRUE, .before = 1)
  ) %>%
    mutate(
      Group   = label,
      dAIC    = AIC - min(AIC, na.rm = TRUE),
      dBIC    = BIC - min(BIC, na.rm = TRUE),
      .before = 1
    ) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))

  ## --- LRT'ler ----------------------------------------------------
  ## M0 - M1 : PV icermez, klasik LRT
  ## M1 - M2 : PV'ler uzerinden D2
  ## M2 - M3 : PV'ler uzerinden D2
  lrt_01 <- pooled_lrt(m1, m0, "Model 0 -> Model 1")
  lrt_12 <- pooled_lrt(m2_fits, m1, "Model 1 -> Model 2")
  lrt_23 <- pooled_lrt(m3_fits, m2_fits, "Model 2 -> Model 3")

  lrt_table <- bind_rows(
    lrt_01$pooled,
    lrt_12$pooled,
    lrt_23$pooled
  ) %>%
    mutate(
      Group = label,
      Sig = case_when(
        is.na(p.value) ~ NA_character_,
        p.value < .001 ~ "***",
        p.value < .01  ~ "**",
        p.value < .05  ~ "*",
        TRUE           ~ ""
      ),
      .before = 1
    ) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))

  lrt_by_pv <- bind_rows(
    lrt_01$by_pv %>% mutate(Comparison = "Model 0 -> Model 1"),
    lrt_12$by_pv %>% mutate(Comparison = "Model 1 -> Model 2"),
    lrt_23$by_pv %>% mutate(Comparison = "Model 2 -> Model 3")
  ) %>%
    mutate(Group = label, .before = 1)

  ## --- Varyans bilesenleri ---------------------------------------
  vc <- bind_rows(
    as.data.frame(lme4::VarCorr(m0)) %>%
      dplyr::select(grp, var1, var2, vcov, sdcor) %>%
      mutate(Model = unname(model_labels["m0"])),
    as.data.frame(lme4::VarCorr(m1)) %>%
      dplyr::select(grp, var1, var2, vcov, sdcor) %>%
      mutate(Model = unname(model_labels["m1"])),
    ran_pars_across(m2_fits) %>%
      mutate(Model = unname(model_labels["m2"])),
    ran_pars_across(m3_fits) %>%
      mutate(Model = unname(model_labels["m3"]))
  ) %>%
    mutate(Group = label, .before = 1)

  ## --- Aciklanan varyans (proportional reduction) ----------------
  get_var <- function(model_label, group_name) {
    v <- vc$vcov[vc$Model == model_label & vc$grp == group_name]
    if (length(v) == 0) NA_real_ else v[1]
  }

  v_item_m0 <- get_var(model_labels["m0"], "Items")
  v_item_m1 <- get_var(model_labels["m1"], "Items")
  v_id_m1   <- get_var(model_labels["m1"], "ID")
  v_id_m2   <- get_var(model_labels["m2"], "ID")

  variance_explained <- tibble(
    Group = label,
    Level = c("Item (Items)", "Person (ID)"),
    Source = c(
      "Complexity_c (Model 0 -> Model 1)",
      "Student-level predictors (Model 1 -> Model 2)"
    ),
    Var_before = c(v_item_m0, v_id_m1),
    Var_after  = c(v_item_m1, v_id_m2),
    R2_level   = c(
      (v_item_m0 - v_item_m1) / v_item_m0,
      (v_id_m1   - v_id_m2)   / v_id_m1
    )
  ) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))

  ## --- Yakinsama kontrolleri -------------------------------------
  checks <- bind_rows(
    check_glmer_model(m0, unname(model_labels["m0"])),
    check_glmer_model(m1, unname(model_labels["m1"])),
    purrr::map_dfr(seq_along(m2_fits), function(i) {
      check_glmer_model(m2_fits[[i]], unname(model_labels["m2"])) %>%
        mutate(PV = i)
    }),
    purrr::map_dfr(seq_along(m3_fits), function(i) {
      check_glmer_model(m3_fits[[i]], unname(model_labels["m3"])) %>%
        mutate(PV = i)
    })
  ) %>%
    mutate(Group = label, .before = 1)

  ## --- Pooled sabit etki tablosu ---------------------------------
  pooled_fixed <- bind_rows(
    pooled_m2 %>% mutate(Model = unname(model_labels["m2"])),
    pooled_m3 %>% mutate(Model = unname(model_labels["m3"]))
  ) %>%
    mutate(
      Group = label,
      OR    = exp(estimate),
      OR_lo = exp(conf.low),
      OR_hi = exp(conf.high),
      Sig = case_when(
        p.value < .001 ~ "***",
        p.value < .01  ~ "**",
        p.value < .05  ~ "*",
        TRUE           ~ ""
      )
    ) %>%
    dplyr::select(
      Group, Model, term, estimate, std.error, df, statistic, p.value, Sig,
      conf.low, conf.high, OR, OR_lo, OR_hi, riv, lambda, fmi, M
    ) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))

  ## --- Model 1'in sabit etkileri (tek model) ---------------------
  m1_fixed <- broom.mixed::tidy(m1, effects = "fixed", conf.int = TRUE) %>%
    mutate(
      Group = label,
      Model = unname(model_labels["m1"]),
      OR    = exp(estimate),
      OR_lo = exp(conf.low),
      OR_hi = exp(conf.high)
    )

  ## --- Kaydet -----------------------------------------------------
  tag <- gsub("[^A-Za-z0-9]+", "_", label)

  write.csv(comparison,
            file.path(out_dir, paste0(tag, "_model_comparison.csv")),
            row.names = FALSE)
  write.csv(lrt_table,
            file.path(out_dir, paste0(tag, "_lrt.csv")),
            row.names = FALSE)
  write.csv(lrt_by_pv,
            file.path(out_dir, paste0(tag, "_lrt_by_pv.csv")),
            row.names = FALSE)
  write.csv(pooled_fixed,
            file.path(out_dir, paste0(tag, "_pooled_fixed_effects.csv")),
            row.names = FALSE)
  write.csv(m1_fixed,
            file.path(out_dir, paste0(tag, "_model1_fixed_effects.csv")),
            row.names = FALSE)
  write.csv(vc,
            file.path(out_dir, paste0(tag, "_variance_components.csv")),
            row.names = FALSE)
  write.csv(variance_explained,
            file.path(out_dir, paste0(tag, "_variance_explained.csv")),
            row.names = FALSE)
  write.csv(checks,
            file.path(out_dir, paste0(tag, "_model_checks.csv")),
            row.names = FALSE)
  write.csv(bind_rows(s2$by_pv %>% mutate(Model = unname(model_labels["m2"])),
                      s3$by_pv %>% mutate(Model = unname(model_labels["m3"]))),
            file.path(out_dir, paste0(tag, "_fit_by_pv.csv")),
            row.names = FALSE)

  capture.output(summary(m0),
                 file = file.path(out_dir, paste0(tag, "_model0_summary.txt")))
  capture.output(summary(m1),
                 file = file.path(out_dir, paste0(tag, "_model1_summary.txt")))
  capture.output(summary(m2_fits[[1]]),
                 file = file.path(out_dir, paste0(tag, "_model2_PV1_summary.txt")))
  capture.output(summary(m3_fits[[1]]),
                 file = file.path(out_dir, paste0(tag, "_model3_PV1_summary.txt")))

  ## --- Ekrana ozet -----------------------------------------------
  cat("\n--- Model karsilastirmasi:", label, "---\n")
  print(as.data.frame(comparison))

  cat("\n--- Likelihood ratio testleri:", label, "---\n")
  print(as.data.frame(lrt_table))

  cat("\n--- Aciklanan varyans:", label, "---\n")
  print(as.data.frame(variance_explained))

  list(
    label              = label,
    n_students         = n_distinct(dat$ID),
    n_responses        = nrow(dat),
    n_countries        = n_distinct(dat$CNT),
    m0                 = m0,
    m1                 = m1,
    m0_raw             = m0_raw,
    m1_raw             = m1_raw,
    m2_fits            = m2_fits,
    m3_fits            = m3_fits,
    pooled_m2          = pooled_m2,
    pooled_m3          = pooled_m3,
    pooled_fixed       = pooled_fixed,
    m1_fixed           = m1_fixed,
    comparison         = comparison,
    fit_by_pv          = bind_rows(s2$by_pv, s3$by_pv),
    lrt_table          = lrt_table,
    lrt_by_pv          = lrt_by_pv,
    variance_components = vc,
    variance_explained = variance_explained,
    checks             = checks
  )
}


####################################################################
## 18. BETIMSEL RASCH MODELI
##
## Madde guclukleri icin; model merdivenine dahil DEGILDIR.
## -1 + Items parametrelemesi madde kolayliklarini verir;
## eirm print(difficulty = TRUE) ile gucluge cevrilir.
####################################################################

message("Betimsel Rasch modeli tahmin ediliyor...")

rasch_model <- fit_model(f_rasch, data_long)
rasch_glmer <- as_glmer(rasch_model)

if (HAS_EIRM) {
  print(rasch_model, difficulty = TRUE)
  try(plot(rasch_model), silent = TRUE)
}

item_difficulty <- tibble(
  term       = names(lme4::fixef(rasch_glmer)),
  easiness   = unname(lme4::fixef(rasch_glmer)),
  std.error  = unname(sqrt(diag(as.matrix(stats::vcov(rasch_glmer)))))
) %>%
  mutate(
    Item       = sub("^Items", "", term),
    difficulty = -easiness
  ) %>%
  left_join(complexity_df, by = "Item") %>%
  dplyr::select(Item, easiness, difficulty, std.error, Complexity, Complexity_c) %>%
  arrange(difficulty) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

print(as.data.frame(item_difficulty))

write.csv(item_difficulty, out_main("item_difficulty_rasch.csv"),
          row.names = FALSE)

saveRDS(rasch_model, out_main("rasch_model.rds"))

capture.output(summary(rasch_glmer),
               file = out_main("rasch_model_summary.txt"))

## Madde guclugu ile ATOS complexity arasindaki iliski
complexity_difficulty_cor <- cor(
  item_difficulty$difficulty,
  item_difficulty$Complexity,
  use = "complete.obs"
)

cat("Madde guclugu - ATOS complexity korelasyonu:",
    round(complexity_difficulty_cor, 3), "\n")


####################################################################
## 19. TUM ORNEK ICIN MODEL 0 - MODEL 3
####################################################################

main_res <- run_model_ladder(data_long, "Full sample", OUT_MAIN)

model_comparison <- main_res$comparison
lrt_table        <- main_res$lrt_table
pooled_report    <- main_res$pooled_fixed

saveRDS(main_res$m0,      out_main("model0.rds"))
saveRDS(main_res$m1,      out_main("model1.rds"))
saveRDS(main_res$m2_fits, out_main("model2_fits.rds"))
saveRDS(main_res$m3_fits, out_main("model3_fits.rds"))


####################################################################
## 20. POOLED SABIT ETKILERIN RAPORLANMASI
####################################################################

cat("\n--- Rubin-pooled sabit etkiler (Model 2 ve Model 3) ---\n")
print(as.data.frame(pooled_report))

write.csv(pooled_report,
          out_main("pooled_fixed_effects_rubin.csv"),
          row.names = FALSE)

## Capraz kontrol: mitml::testEstimates ayni sonucu vermeli
mitml_check <- tryCatch({
  capture.output(
    {
      cat("### Model 2 ###\n")
      print(mitml::testEstimates(main_res$m2_fits, extra.pars = TRUE))
      cat("\n\n### Model 3 ###\n")
      print(mitml::testEstimates(main_res$m3_fits, extra.pars = TRUE))
    },
    file = out_main("pooled_fixed_effects_mitml_check.txt")
  )
  TRUE
}, error = function(e) {
  message("mitml::testEstimates capraz kontrolu basarisiz: ",
          conditionMessage(e))
  FALSE
})


####################################################################
## 21. VARYANS BILESENLERI VE ACIKLANAN VARYANS
####################################################################

var_components     <- main_res$variance_components
item_R2            <- main_res$variance_explained

cat("\n--- Varyans bilesenleri ---\n")
print(as.data.frame(var_components))

cat("\n--- Aciklanan varyans ---\n")
print(as.data.frame(item_R2))


####################################################################
## 22. COMPLEXITY ETKISININ MODELLERE GORE OZETI
####################################################################

complexity_effect_table <- bind_rows(
  ## Model 1 tek modeldir: p degeri Wald testinden gelir, df yoktur.
  main_res$m1_fixed %>%
    filter(term == "Complexity_c") %>%
    dplyr::select(Model, term, estimate, std.error, p.value,
                  conf.low, conf.high, OR, OR_lo, OR_hi) %>%
    mutate(df = NA_real_),

  main_res$pooled_fixed %>%
    filter(term == "Complexity_c") %>%
    dplyr::select(Model, term, estimate, std.error, df, p.value,
                  conf.low, conf.high, OR, OR_lo, OR_hi)
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

cat("\n--- Complexity_c etkisinin modellere gore ozeti ---\n")
print(as.data.frame(complexity_effect_table))

write.csv(complexity_effect_table,
          out_main("complexity_effect_table.csv"),
          row.names = FALSE)


####################################################################
## 23. POOLED MARJINAL ETKI GRAFIKLERI
##
## Grafikler POOLED katsayilar ve POOLED kovaryans matrisi
## kullanilarak cizilir (tek bir PV modeli degil).
####################################################################

## Pooled tahmin ve SE hesaplayan yardimci ------------------------
pooled_predict <- function(fit_list, grid) {

  pooled <- rubin_total_vcov(fit_list)
  b_bar  <- pooled$b
  Tvar   <- pooled$T

  tt <- stats::delete.response(stats::terms(fit_list[[1]], fixed.only = TRUE))
  X  <- stats::model.matrix(tt, grid)

  missing_cols <- setdiff(names(b_bar), colnames(X))
  if (length(missing_cols) > 0) {
    stop("Tahmin matrisinde eksik kolonlar: ",
         paste(missing_cols, collapse = ", "))
  }

  X <- X[, names(b_bar), drop = FALSE]

  eta <- drop(X %*% b_bar)
  se  <- sqrt(rowSums((X %*% Tvar) * X))

  grid$eta   <- eta
  grid$se    <- se
  grid$prob  <- stats::plogis(eta)
  grid$lower <- stats::plogis(eta - 1.96 * se)
  grid$upper <- stats::plogis(eta + 1.96 * se)
  grid
}

## Ortak kovaryans degerleri --------------------------------------
ref_values <- list(
  ESCS     = mean(person_level$ESCS, na.rm = TRUE),
  Homework = 0,   ## standardize edildigi icin ortalama = 0
  Language = factor("Test",   levels = levels(data_long$Language)),
  Gender   = factor("Female", levels = levels(data_long$Gender))
)

complexity_seq <- seq(
  min(data_long$Complexity_c, na.rm = TRUE),
  max(data_long$Complexity_c, na.rm = TRUE),
  length.out = 100
)

## --- Model 2: Complexity ana etkisi -----------------------------
grid_m2 <- data.frame(
  Complexity_c = complexity_seq,
  PVREAD       = 0,
  ESCS         = ref_values$ESCS,
  Homework     = ref_values$Homework,
  Language     = ref_values$Language,
  Gender       = ref_values$Gender
)

grid_m2 <- pooled_predict(main_res$m2_fits, grid_m2)

plot_m2 <- ggplot(grid_m2, aes(x = Complexity_c, y = prob)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title    = "Model 2: item complexity effect",
    subtitle = "Rubin-pooled across 10 READ plausible values",
    x        = "Item complexity (centered)",
    y        = "Predicted probability of a correct response"
  ) +
  theme_bw()

ggsave(out_main("model2_complexity_effect.png"), plot_m2,
       width = 7, height = 5, dpi = 300)

## --- Model 3: Complexity x PVREAD etkilesimi --------------------
grid_m3 <- expand.grid(
  Complexity_c = complexity_seq,
  PVREAD       = c(-1, 0, 1)
)

grid_m3$ESCS     <- ref_values$ESCS
grid_m3$Homework <- ref_values$Homework
grid_m3$Language <- ref_values$Language
grid_m3$Gender   <- ref_values$Gender

grid_m3 <- pooled_predict(main_res$m3_fits, grid_m3)

grid_m3$PVREAD_f <- factor(
  grid_m3$PVREAD,
  levels = c(-1, 0, 1),
  labels = c("-1 SD", "Mean", "+1 SD")
)

plot_m3 <- ggplot(
  grid_m3,
  aes(x = Complexity_c, y = prob,
      group = PVREAD_f, linetype = PVREAD_f, fill = PVREAD_f)
) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title    = "Model 3: item complexity x reading achievement",
    subtitle = "Rubin-pooled across 10 READ plausible values",
    x        = "Item complexity (centered)",
    y        = "Predicted probability of a correct response",
    linetype = "Reading (PVREAD)",
    fill     = "Reading (PVREAD)"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(out_main("model3_complexity_by_reading.png"), plot_m3,
       width = 7, height = 5, dpi = 300)

write.csv(grid_m2, out_main("model2_predicted_probabilities.csv"),
          row.names = FALSE)
write.csv(grid_m3, out_main("model3_predicted_probabilities.csv"),
          row.names = FALSE)

## eirm'e ozgu marjinal etki grafigi (varsa).
## Model yeniden tahmin edilmez; merdivende tahmin edilen
## Model 1 nesnesi kullanilir.
if (HAS_EIRM) {
  try(
    print(eirm::marginalplot(main_res$m1_raw, predictors = "Complexity_c")),
    silent = TRUE
  )
}


####################################################################
## 24. MODEL 0 - MODEL 1 ICIN HTML TABLO
##
## Model 2 ve Model 3 pooled sonuclari CSV dosyalarindadir;
## tab_model tek bir modeli gosterebildigi icin buraya alinmaz.
####################################################################

if (HAS_SJP) {
  try({
    sjPlot::tab_model(
      main_res$m0,
      main_res$m1,
      dv.labels      = c(unname(model_labels["m0"]), unname(model_labels["m1"])),
      show.aic       = TRUE,
      show.icc       = TRUE,
      show.obs       = TRUE,
      show.se        = TRUE,
      show.p         = FALSE,
      show.intercept = TRUE,
      p.style        = "stars",
      show.r2        = TRUE,
      show.loglik    = TRUE,
      digits         = 3,
      file           = out_main("model0_model1.html")
    )
  }, silent = TRUE)
}


####################################################################
## 25. BETIMSEL AGIRLIKLI KORELASYONLAR
##
## Modellerden farkli olarak burada ornekleme agirliklari
## KULLANILIR (betimsel istatistik oldugu icin).
####################################################################

person_dat <- data_long %>%
  distinct(ID, .keep_all = TRUE) %>%
  dplyr::select(ID, CNT, WEIGHT, ESCS, Homework, all_of(pv_z_names)) %>%
  mutate(
    ESCS     = as.numeric(ESCS),
    Homework = as.numeric(Homework),
    WEIGHT   = as.numeric(WEIGHT),
    across(all_of(pv_z_names), as.numeric)
  )

cat("\nBetimsel korelasyonlar icin kisi sayisi:", nrow(person_dat), "\n")

## Agirlikli korelasyon -------------------------------------------
wtd_cor <- function(x, y, w) {

  x <- as.numeric(x); y <- as.numeric(y); w <- as.numeric(w)

  ok <- stats::complete.cases(x, y, w) & w > 0
  x  <- x[ok]; y <- y[ok]; w <- w[ok]

  n_obs <- length(x)
  if (n_obs < 4L) return(list(r = NA_real_, n = n_obs))

  mx <- sum(w * x) / sum(w)
  my <- sum(w * y) / sum(w)

  cov_xy <- sum(w * (x - mx) * (y - my)) / sum(w)
  sx     <- sqrt(sum(w * (x - mx)^2) / sum(w))
  sy     <- sqrt(sum(w * (y - my)^2) / sum(w))

  if (sx == 0 || sy == 0) return(list(r = NA_real_, n = n_obs))

  r_value <- cov_xy / (sx * sy)
  r_value <- max(min(r_value, 0.999999), -0.999999)

  list(r = r_value, n = n_obs)
}

## PV korelasyonlarinin Fisher-z olceginde Rubin ile pool edilmesi
pool_correlation_pv <- function(r_vec, n_vec, conf = 0.95) {

  valid <- !is.na(r_vec) & !is.na(n_vec) & n_vec > 3 & abs(r_vec) < 1
  r_vec <- r_vec[valid]; n_vec <- n_vec[valid]

  if (length(r_vec) == 0L) {
    return(tibble(r = NA_real_, conf.low = NA_real_, conf.high = NA_real_,
                  p.value = NA_real_, df = NA_real_, fmi = NA_real_, M = 0L))
  }

  z_vec <- atanh(r_vec)
  u_vec <- 1 / (n_vec - 3)

  if (length(z_vec) == 1L) {
    qbar      <- z_vec[1]
    total_var <- u_vec[1]
    df_value  <- n_vec[1] - 3
    fmi_value <- 0
  } else {
    ps <- mice::pool.scalar(Q = z_vec, U = u_vec, n = mean(n_vec), k = 1)
    qbar      <- ps$qbar
    total_var <- ps$t
    df_value  <- ps$df
    fmi_value <- ps$f
  }

  se   <- sqrt(total_var)
  crit <- stats::qt(1 - (1 - conf) / 2, df = df_value)
  stat <- qbar / se

  tibble(
    r         = tanh(qbar),
    conf.low  = tanh(qbar - crit * se),
    conf.high = tanh(qbar + crit * se),
    p.value   = 2 * stats::pt(-abs(stat), df = df_value),
    df        = df_value,
    fmi       = fmi_value,
    M         = length(z_vec)
  )
}

cont_all <- c("PVREAD", "ESCS", "Homework")
n_cont   <- length(cont_all)

R_pool <- matrix(1, n_cont, n_cont, dimnames = list(cont_all, cont_all))
P_pool <- matrix(NA_real_, n_cont, n_cont, dimnames = list(cont_all, cont_all))
diag(P_pool) <- 0

cor_long <- list()

for (i in seq_len(n_cont)) {
  for (j in seq_len(n_cont)) {

    if (j <= i) next

    v1 <- cont_all[i]
    v2 <- cont_all[j]

    if (v1 == "PVREAD" || v2 == "PVREAD") {

      other <- if (v1 == "PVREAD") v2 else v1

      cc <- purrr::map(pv_z_names, ~ wtd_cor(person_dat[[.x]],
                                             person_dat[[other]],
                                             person_dat$WEIGHT))

      out <- pool_correlation_pv(
        r_vec = purrr::map_dbl(cc, "r"),
        n_vec = purrr::map_dbl(cc, "n")
      )

    } else {

      cc  <- wtd_cor(person_dat[[v1]], person_dat[[v2]], person_dat$WEIGHT)
      out <- pool_correlation_pv(r_vec = cc$r, n_vec = cc$n)
    }

    R_pool[i, j] <- R_pool[j, i] <- out$r
    P_pool[i, j] <- P_pool[j, i] <- out$p.value

    cor_long[[length(cor_long) + 1]] <- out %>%
      mutate(Var1 = v1, Var2 = v2, .before = 1)
  }
}

cor_table <- bind_rows(cor_long) %>%
  mutate(
    Sig = case_when(
      p.value < .001 ~ "***",
      p.value < .01  ~ "**",
      p.value < .05  ~ "*",
      TRUE           ~ ""
    )
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

cat("\n--- Rubin-pooled implied correlation matrix ---\n")
print(round(R_pool, 3))

cat("\n--- P-values ---\n")
print(round(P_pool, 4))

cat("\n--- Pairwise detail ---\n")
print(as.data.frame(cor_table))

write.csv(round(R_pool, 3), out_main("implied_correlation_matrix_rubin.csv"))
write.csv(round(P_pool, 4), out_main("implied_correlation_pvalues_rubin.csv"))
write.csv(cor_table, out_main("correlation_pairs_rubin.csv"), row.names = FALSE)


####################################################################
## 26. ULKE BAZINDA POOLED KORELASYON (PVREAD - ESCS)
####################################################################

cor_by_country <- person_dat %>%
  group_split(CNT) %>%
  purrr::map_dfr(function(d) {

    cc <- purrr::map(pv_z_names, ~ wtd_cor(d[[.x]], d$ESCS, d$WEIGHT))

    pool_correlation_pv(
      r_vec = purrr::map_dbl(cc, "r"),
      n_vec = purrr::map_dbl(cc, "n")
    ) %>%
      mutate(CNT = as.character(d$CNT[1]), Pair = "PVREAD-ESCS", .before = 1)
  }) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

print(as.data.frame(cor_by_country))

write.csv(cor_by_country,
          out_main("correlation_PVREAD_ESCS_by_country.csv"),
          row.names = FALSE)


####################################################################
## 27. EO GRUP ANALIZLERI
##
## Ayni Model 0 - Model 3 merdiveni her EO grubunda tekrarlanir.
## PV'ler tum ornek uzerinden standardize edildigi icin gruplarin
## katsayilari ayni olcektedir.
####################################################################

dat_eo_formal    <- data_long %>% filter(EO == "Formal")    %>% droplevels()
dat_eo_nonformal <- data_long %>% filter(EO == "Non-Formal") %>% droplevels()

cat("\nFormal EO birey sayisi:",    n_distinct(dat_eo_formal$ID), "\n")
cat("Non-Formal EO birey sayisi:", n_distinct(dat_eo_nonformal$ID), "\n")

eo_formal_res    <- run_model_ladder(dat_eo_formal,    "Formal",    OUT_EO)
eo_nonformal_res <- run_model_ladder(dat_eo_nonformal, "NonFormal", OUT_EO)

saveRDS(
  list(Formal = eo_formal_res, NonFormal = eo_nonformal_res),
  out_eo("EO_results.rds")
)

## --- EO pooled beta tablosu --------------------------------------
beta_table_eo <- bind_rows(
  eo_formal_res$pooled_fixed,
  eo_nonformal_res$pooled_fixed
) %>%
  rename(EnglishOfficial = Group)

print(as.data.frame(beta_table_eo))

write.csv(beta_table_eo, out_eo("pooled_beta_by_EO.csv"), row.names = FALSE)

## --- EO model karsilastirma tablosu ------------------------------
eo_model_comparison <- bind_rows(
  eo_formal_res$comparison,
  eo_nonformal_res$comparison
) %>%
  rename(EnglishOfficial = Group)

print(as.data.frame(eo_model_comparison))

write.csv(eo_model_comparison, out_eo("EO_model_comparison.csv"),
          row.names = FALSE)

## --- EO LRT tablosu ----------------------------------------------
eo_lrt <- bind_rows(
  eo_formal_res$lrt_table,
  eo_nonformal_res$lrt_table
) %>%
  rename(EnglishOfficial = Group)

print(as.data.frame(eo_lrt))

write.csv(eo_lrt, out_eo("EO_lrt.csv"), row.names = FALSE)


####################################################################
## 28. FOREST PLOT (odds ratio)
####################################################################

forest_data <- beta_table_eo %>%
  filter(
    Model == unname(model_labels["m3"]),
    term  != "(Intercept)"
  ) %>%
  mutate(
    term      = factor(term, levels = rev(unique(term))),
    Sig_group = ifelse(p.value < .05, "p < .05", "p >= .05")
  )

forest_plot <- ggplot(
  forest_data,
  aes(x = OR, y = term, xmin = OR_lo, xmax = OR_hi,
      colour = Sig_group, shape = Sig_group)
) +
  geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
  geom_errorbarh(height = 0.15) +
  geom_point(size = 3) +
  facet_wrap(~ EnglishOfficial, scales = "free_x") +
  labs(
    title    = "Model 3: pooled coefficients across 10 READ plausible values",
    subtitle = "Odds ratios with Rubin-pooled 95% confidence intervals (unweighted models)",
    x        = "Odds ratio",
    y        = NULL,
    shape    = "Significance",
    colour   = "Significance"
  ) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold"), legend.position = "bottom")

print(forest_plot)

ggsave(out_eo("EO_model3_forest_plot.png"), forest_plot,
       width = 12, height = 8, dpi = 300)
ggsave(out_eo("EO_model3_forest_plot.pdf"), forest_plot,
       width = 12, height = 8)


####################################################################
## 29. TUM MODELLER ICIN YAKINSAMA KONTROLLERI
####################################################################

all_checks <- bind_rows(
  check_glmer_model(rasch_glmer, "Rasch (descriptive)") %>%
    mutate(Group = "Full sample", .before = 1),
  main_res$checks,
  eo_formal_res$checks,
  eo_nonformal_res$checks
)

print(as.data.frame(all_checks))

write.csv(all_checks, out_main("all_model_checks.csv"), row.names = FALSE)

if (any(all_checks$Singular)) {
  warning(
    "Singular uyum veren modeller var: ",
    paste(unique(all_checks$Model[all_checks$Singular]), collapse = ", ")
  )
}


####################################################################
## 30. CALISMA CIKTILARININ KAYDEDILMESI
####################################################################

analysis_objects <- list(
  item_response_counts            = item_response_counts,
  item_response_counts_by_country = item_response_counts_by_country,
  country_counts                  = country_counts,
  country_counts_by_EO            = country_counts_by_EO,
  raw_country_counts              = raw_country_counts,
  complexity_df                   = complexity_df,
  complexity_mean                 = complexity_mean,
  pv_moments                      = pv_moments,
  item_difficulty                 = item_difficulty,
  complexity_difficulty_cor       = complexity_difficulty_cor,
  item_R2                         = item_R2,
  var_components                  = var_components,
  complexity_effect_table         = complexity_effect_table,
  pooled_report                   = pooled_report,
  model_comparison                = model_comparison,
  lrt_table                       = lrt_table,
  lrt_by_pv                       = main_res$lrt_by_pv,
  R_pool                          = R_pool,
  P_pool                          = P_pool,
  cor_table                       = cor_table,
  cor_by_country                  = cor_by_country,
  beta_table_eo                   = beta_table_eo,
  eo_model_comparison             = eo_model_comparison,
  eo_lrt                          = eo_lrt,
  all_model_checks                = all_checks
)

saveRDS(analysis_objects, out_main("analysis_summary_objects.rds"))


####################################################################
## 31. SESSION INFO
####################################################################

capture.output(sessionInfo(), file = out_main("sessionInfo.txt"))

cat("\n============================================================\n")
cat("Analiz tamamlandi.\n")
cat("Ana ciktilar:", OUT_MAIN, "\n")
cat("EO grup ciktilari:", OUT_EO, "\n")
cat("Model tahmin motoru:", if (HAS_EIRM) "eirm" else "lme4::glmer", "\n")
cat("Homework sinifi:", class(data_long$Homework), "\n")
cat("PV1READ_z sinifi:", class(data_long$PV1READ_z), "\n")
cat("============================================================\n")
