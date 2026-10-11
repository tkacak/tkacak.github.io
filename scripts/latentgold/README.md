# Latent GOLD — faktör karışım modeli (K = 1–5)

`fmm_openmx_enumeration.R` modelinin Latent GOLD 6.x Syntax karşılığı. Dosyalar
`make_latentgold_files.R` ile üretildi; aynı veri (N = 1000, seed 20261011).

| Dosya | İçerik |
|---|---|
| `fmm_data_lg.txt` | Veri. İlk satır değişken adları (`id y1 … y10 trueclass`) |
| `fmm_K1_check_n10.lgs`, `fmm_K1_check_n50.lgs` | **Önce bunları çalıştır** (tanı) |
| `fmm_K2.lgs` | K = 2, sınıflama dosyası yazar |
| `fmm_K1to5_range.lgs` | `Class nominal 1:5`: BLRT (`LLdiff`) ve VLMR tek çalıştırmada |
| `fmm_K1to5_blocks.lgs` | Yedek: range kabul edilmezse K = 1..5 ayrı bloklar |

Hepsini aynı klasöre koy. Sürüm 6.1 ise ilk iki satırı `//LG6.1//` ve `version = 6.1` yap.

## Sıra

1. **Tanı (K = 1, iki dosya).** LG'nin `Log-likelihood (LL)` değeri OpenMx ile eşleşmeli: **−12638.387**.
   - İkisinde de −12638.387 → LG olabilirliği tam hesaplıyor; `nodes` önemsiz.
   - nodes=10'da farklı, nodes=50'de −12638.387 → LG sayısal integrasyon yapıyor; tüm çalıştırmalarda `nodes=50` kalmalı.
   - İkisinde de farklı → model eşdeğer değil; `Npar` = 31 mi ve faktör varyansları çıktıda tam 1 mi, kontrol et.
2. **`fmm_K2.lgs`.** Beklenen: LL = −12636.327, Npar = 34, sınıf oranları ≈ .555/.445,
   sınıf 2 ortalamaları ≈ (0.254, −1.211). Popülasyon değerlerini (.80/.20, 1/1) **bekleme**;
   bu örneklemin ML çözümü budur.
3. **`fmm_K1to5_range.lgs`.** Hata verirse `fmm_K1to5_blocks.lgs`.

## OpenMx referans değerleri

| K | LL | Npar | BIC |
|---|---|---|---|
| 1 | −12638.387 | 31 | 25490.91 |
| 2 | −12636.327 | 34 | 25507.52 |
| 3 | −12632.972 | 37 | 25521.53 |
| 4 | −12631.971 | 40 | 25540.25 |
| 5 | −12629.334 | 43 | 25555.70 |

## Bilinen riskler

- `(1) F1;` faktör varyansını 1'e sabitler; belgelerden doğrulanamadı. Reddedilirse ya da
  varyanslar çıktıda 1 görünmezse: `(1) F1; (1) F2;` yerine `F1; F2;` yaz,
  `y1 <- 1 + (1) F1;` ve `y6 <- 1 + (1) F2;` yap (aynı olabilirlik, aynı Npar).
- VLMR'nin faktör karışım modellerinde (sürekli faktörlü) üretildiği belgelenmedi.
- LG'nin "Entropy R-squared" değeri Mplus/OpenMx entropisinden farklı formülle hesaplanır.
- K = 5'te faktör korelasyonu ≈ .97; OpenMx'te üst sınır .99. LG daha yüksek LL bulursa önce buna bak.
