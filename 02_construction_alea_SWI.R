# ============================================================
# 02_construction_alea_SWI.R
# Recognition-rate series (2002-2023), regime-break tests, and
# construction of the monthly communal SWI (area-weighted).
# Run from the repository root (see README for the data/ layout).
# ============================================================
library(dplyr); library(tidyr); library(lubridate); library(stringr)
library(readxl); library(Matrix)

dir_out <- file.path("output", "etape0")
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Step 1: annual recognition series, seasonal vs annual windows
# ------------------------------------------------------------
base <- readRDS(file.path(dir_out, "base_evenements_couts.rds"))

ev_m <- base %>%
  filter(annee >= 2002, annee <= 2023) %>%
  mutate(m_deb = month(debut), m_fin = month(fin),
         annuelle = m_deb == 1 & m_fin == 12,
         type_dem = ifelse(annuelle, "annuelle", "trimestrielle"))  # "trimestrielle" = seasonal window
cat("Rows 2002-2023:", nrow(ev_m), "| annual-window:", sum(ev_m$annuelle), "\n")

ser2 <- ev_m %>%
  group_by(annee, type_dem) %>%
  summarise(demandes = n(), reconnues = sum(reconnu), .groups = "drop") %>%
  pivot_wider(names_from = type_dem, values_from = c(demandes, reconnues), values_fill = 0) %>%
  mutate(taux_hors_annuelle = reconnues_trimestrielle / demandes_trimestrielle,
         taux_global = (reconnues_trimestrielle + reconnues_annuelle) /
           (demandes_trimestrielle + demandes_annuelle))

print(round(ser2 %>% filter(annee %in% c(2003, 2018, 2022)) %>%
              select(annee, taux_hors_annuelle, taux_global), 3))
saveRDS(ser2, file.path(dir_out, "serie_reconnaissance_hors_annuelle.rds"))

# ------------------------------------------------------------
# Step 2: tests of a break in the recognition regime
# ------------------------------------------------------------
d4 <- ser2 %>% filter(annee >= 2003) %>%
  mutate(rec = reconnues_trimestrielle,
         non = demandes_trimestrielle - reconnues_trimestrielle,
         taux = rec / (rec + non),
         apres17 = annee >= 2017)

m0 <- glm(cbind(rec, non) ~ 1, data = d4, family = quasibinomial)
cat("Null deviance:", round(deviance(m0)), "| dispersion:", round(summary(m0)$dispersion, 1), "\n")

print(d4 %>% group_by(apres17) %>%
        summarise(n_annees = n(), moyenne = round(mean(taux), 3), mediane = round(median(taux), 3),
                  min = round(min(taux), 3), max = round(max(taux), 3)))
print(wilcox.test(taux ~ apres17, data = d4, exact = FALSE))

# Largest deviance gain over candidate break dates c0 in 2008-2020
gain_max <- function(dat) {
  m00 <- glm(cbind(rec, non) ~ 1, data = dat, family = quasibinomial)
  max(sapply(2008:2020, function(c0) {
    m <- glm(cbind(rec, non) ~ I(annee >= c0), data = dat, family = quasibinomial)
    deviance(m00) - deviance(m)
  }))
}
obs <- gain_max(d4)
set.seed(7)
perm <- replicate(499, {          # permutation test on the order of the years
  dd <- d4; ordre <- sample(nrow(dd))
  dd$rec <- d4$rec[ordre]; dd$non <- d4$non[ordre]
  gain_max(dd)
})
cat("Observed peak:", round(obs), "| quantiles 50/95/99 %:", round(quantile(perm, c(.5, .95, .99))),
    "| share >= observed:", mean(perm >= obs), "\n")

# ------------------------------------------------------------
# Step 3: Beta parameters of the recognition rate by regime
# ------------------------------------------------------------
par_reg <- d4 %>% group_by(apres17) %>%
  summarise(n = n(), m = mean(taux), v = var(taux), sd = sd(taux)) %>%
  mutate(k = m * (1 - m) / v - 1, alpha = m * k, beta = (1 - m) * k)
print(round(as.data.frame(par_reg), 3))
saveRDS(par_reg, file.path(dir_out, "parametres_regimes_reconnaissance.rds"))

# ------------------------------------------------------------
# Step 4: SWI inputs (checks)
# ------------------------------------------------------------
f_swi <- file.path("data", "swi", "SWI_2000_2025.csv")
f_map <- file.path("data", "swi", "mapping_SWI_communes_weighted.csv")
f_evt <- file.path("data", "catnat", "CatNat_drought_events_SWI.xlsx")
f_py  <- "build_swi_communal.py"   # script that builds the commune-grid cell mapping
stopifnot(file.exists(f_swi), file.exists(f_map), file.exists(f_evt), file.exists(f_py))

cat("--- SWI by grid cell: first 4 lines ---\n");  cat(readLines(f_swi, n = 4), sep = "\n")
cat("\n--- commune-grid cell mapping: first 4 lines ---\n"); cat(readLines(f_map, n = 4), sep = "\n")

evt <- read_excel(f_evt, n_max = 5)
cat("\nColumns of the events table:\n"); print(names(evt))
print(evt)

swi <- read.csv(f_swi)
map <- read.csv(f_map, colClasses = c(code_insee = "character", num_maille = "integer"))

cat("SWI:", nrow(swi), "rows | grid cells:", n_distinct(swi$NUMERO),
    "| dates:", n_distinct(swi$DATE), "| period:", min(swi$DATE), "->", max(swi$DATE), "\n")
cat("Duplicates (grid cell, date):", sum(duplicated(swi[, c("NUMERO", "DATE")])), "\n")
cat("Missing SWI values:", sum(is.na(swi$SWI_UNIF_MENS)), "\n")
print(summary(swi$SWI_UNIF_MENS))

cat("\nMapping:", nrow(map), "rows | communes:", n_distinct(map$code_insee),
    "| grid cells:", n_distinct(map$num_maille), "\n")
cat("INSEE codes of length 5:", mean(nchar(map$code_insee) == 5), "\n")
s <- tapply(map$weight, map$code_insee, sum)
cat("Sum of weights by commune:\n"); print(summary(s))
cat("Communes whose weights deviate from 1 by more than 1e-6:", sum(abs(s - 1) > 1e-6), "\n")
cat("Mapping grid cells absent from the SWI:",
    length(setdiff(unique(map$num_maille), unique(swi$NUMERO))), "\n")

# ------------------------------------------------------------
# Step 5: communal SWI by area weighting, B = W S (renormalised weights)
# ------------------------------------------------------------
maille_ids <- sort(unique(swi$NUMERO)); dates <- sort(unique(swi$DATE))
S <- matrix(NA_real_, length(maille_ids), length(dates), dimnames = list(maille_ids, dates))
S[cbind(match(swi$NUMERO, maille_ids), match(swi$DATE, dates))] <- swi$SWI_UNIF_MENS
cat("Missing values in S:", sum(is.na(S)), "| dim:", dim(S), "\n")

com_ids <- sort(unique(map$code_insee))
W <- sparseMatrix(i = match(map$code_insee, com_ids), j = match(map$num_maille, maille_ids),
                  x = map$weight, dims = c(length(com_ids), length(maille_ids)),
                  dimnames = list(com_ids, maille_ids))
Wsum <- rowSums(W)

A <- as.matrix(W %*% S)     # weighted sum without renormalisation (earlier version)
B <- A / Wsum               # renormalised version (used in the paper)
cat("dim A:", dim(A), "\n")

print(round(quantile(Wsum, c(0, .01, .05, .5, .95, .99, 1)), 3))
cat("Communes with weight sum < 0.99:", sum(Wsum < 0.99), "| < 0.90:", sum(Wsum < 0.90), "\n")
print(round(head(sort(Wsum), 10), 3))
cat("Gap |A - B|: max", round(max(abs(A - B)), 3), "| mean", round(mean(abs(A - B)), 4), "\n")

# Check on one commune: the 3 months of the 2003 event in 01004
print(round(A["01004", c("200307", "200308", "200309")], 3))
print(round(B["01004", c("200307", "200308", "200309")], 3))
evt <- read_excel(f_evt)
print(evt %>% filter(code_insee == "01004") %>%
        select(`Début d'Evénement`, `Fin d'Evénement`, swi_min, month_swi_min) %>% head(6))

# ------------------------------------------------------------
# Step 6: validation against the reference event-level SWI minimum
# ------------------------------------------------------------
evt <- evt %>% mutate(insee = str_pad(code_insee, 5, pad = "0"),
                      d0 = format(floor_date(as.Date(`Début d'Evénement`), "month"), "%Y%m"),
                      d1 = format(floor_date(as.Date(`Fin d'Evénement`), "month"), "%Y%m"))
cat("Events:", nrow(evt), "| dates outside 2000-2025:", sum(evt$d0 < "200001" | evt$d1 > "202512"), "\n")

min_evt <- function(M) {
  ok <- evt$insee %in% rownames(M)
  out <- rep(NA_real_, nrow(evt))
  for (k in which(ok)) {
    cols <- as.character(dates[dates >= as.integer(evt$d0[k]) & dates <= as.integer(evt$d1[k])])
    if (length(cols)) out[k] <- min(M[evt$insee[k], cols])
  }
  out
}
evt$swi_A <- min_evt(A); evt$swi_B <- min_evt(B)

ref <- evt$swi_min
cat("Events without reference:", sum(is.na(ref)), "| without A computation:", sum(is.na(evt$swi_A)), "\n")
for (v in c("swi_A", "swi_B")) {
  e <- abs(evt[[v]] - ref); e <- e[!is.na(e)]
  cat(v, ": max gap", signif(max(e), 3), "| mean", signif(mean(e), 3),
      "| share < 0.001:", round(mean(e < 0.001), 4), "| share < 0.01:", round(mean(e < 0.01), 4), "\n")
}
print(head(evt[which(abs(evt$swi_A - ref) > 0.01), c("insee", "d0", "d1", "swi_min", "swi_A", "swi_B")], 10))

# ------------------------------------------------------------
# Step 7: outputs
# ------------------------------------------------------------
# 1) Monthly communal SWI, renormalised weights (version B): 34,741 communes x 312 months
saveRDS(B, file.path(dir_out, "SWI_communal_mensuel.rds"))

# 2) Communes whose commune-grid cell weights do not sum to 1 (to monitor)
communes_poids_faibles <- data.frame(insee = com_ids, somme_poids = as.numeric(Wsum)) %>%
  filter(somme_poids < 0.99)
saveRDS(communes_poids_faibles, file.path(dir_out, "SWI_communes_poids_faibles.rds"))

cat("SWI:", round(file.size(file.path(dir_out, "SWI_communal_mensuel.rds")) / 1e6), "MB |",
    "communes to monitor:", nrow(communes_poids_faibles), "\n")
list.files(dir_out)

