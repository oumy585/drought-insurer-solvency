library(dplyr); library(tidyr); library(lubridate)
dir_out <- "C:/Users/OumyDIONE/Documents/output/etape0"
base <- readRDS(file.path(dir_out, "base_evenements_couts.rds"))

ev_m <- base %>%
  filter(annee >= 2002, annee <= 2023) %>%
  mutate(m_deb = month(debut), m_fin = month(fin),
         annuelle = m_deb == 1 & m_fin == 12,
         type_dem = ifelse(annuelle, "annuelle", "trimestrielle"))
cat("Lignes 2002-2023 :", nrow(ev_m), "| annuelles :", sum(ev_m$annuelle), "\n")

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


################################# Etape 2 

d4 <- ser2 %>% filter(annee >= 2003) %>%
  mutate(rec = reconnues_trimestrielle,
         non = demandes_trimestrielle - reconnues_trimestrielle,
         taux = rec / (rec + non),
         apres17 = annee >= 2017)

m0 <- glm(cbind(rec, non) ~ 1, data = d4, family = quasibinomial)
cat("Deviance nulle :", round(deviance(m0)), "| dispersion :", round(summary(m0)$dispersion, 1), "\n")

print(d4 %>% group_by(apres17) %>%
        summarise(n_annees = n(), moyenne = round(mean(taux), 3), mediane = round(median(taux), 3),
                  min = round(min(taux), 3), max = round(max(taux), 3)))
print(wilcox.test(taux ~ apres17, data = d4, exact = FALSE))

gain_max <- function(dat) {
  m00 <- glm(cbind(rec, non) ~ 1, data = dat, family = quasibinomial)
  max(sapply(2008:2020, function(c0) {
    m <- glm(cbind(rec, non) ~ I(annee >= c0), data = dat, family = quasibinomial)
    deviance(m00) - deviance(m)
  }))
}
obs <- gain_max(d4)
set.seed(7)
perm <- replicate(499, {
  dd <- d4; ordre <- sample(nrow(dd))
  dd$rec <- d4$rec[ordre]; dd$non <- d4$non[ordre]
  gain_max(dd)
})
cat("Pic observe :", round(obs), "| quantiles 50/95/99 % :", round(quantile(perm, c(.5, .95, .99))),
    "| part >= observe :", mean(perm >= obs), "\n")

############################step 3 
par_reg <- d4 %>% group_by(apres17) %>%
  summarise(n = n(), m = mean(taux), v = var(taux), sd = sd(taux)) %>%
  mutate(k = m * (1 - m) / v - 1, alpha = m * k, beta = (1 - m) * k)
print(round(as.data.frame(par_reg), 3))
saveRDS(par_reg, file.path(dir_out, "parametres_regimes_reconnaissance.rds"))

######################################## etape with swi 
library(dplyr); library(readxl)

root  <- "C:/Users/OumyDIONE/Downloads/Drought_risk_review"
f_swi <- file.path(root, "DATA/SWI_2000_2025.csv")
f_map <- file.path(root, "DATA/mapping_SWI_communes_weighted.csv")
f_evt <- file.path(root, "DATA/CatNat_drought_events_SWI.xlsx")
f_py  <- file.path(root, "code/build_swi_communal.py")
stopifnot(file.exists(f_swi), file.exists(f_map), file.exists(f_evt), file.exists(f_py))

cat("--- SWI par maille : 4 premieres lignes ---\n");  cat(readLines(f_swi, n = 4), sep = "\n")
cat("\n--- correspondance commune-maille : 4 premieres lignes ---\n"); cat(readLines(f_map, n = 4), sep = "\n")

py <- readLines(f_py, warn = FALSE)
cat("\nLignes du script :", length(py), "\n")
idx <- grep("monthly_swi|weight|groupby|read_csv|to_csv|to_parquet|nearest|fillna|sum\\(|merge|swi",
            py, ignore.case = TRUE)
print(data.frame(ligne = idx, code = trimws(py[idx])), right = FALSE, row.names = FALSE)

evt <- read_excel(f_evt, n_max = 5)
cat("\nColonnes de la table des evenements :\n"); print(names(evt))
print(evt)


##############################################
library(dplyr)
root  <- "C:/Users/OumyDIONE/Downloads/Drought_risk_review"
f_swi <- file.path(root, "DATA/SWI_2000_2025.csv")
f_map <- file.path(root, "DATA/mapping_SWI_communes_weighted.csv")

swi <- read.csv(f_swi)
map <- read.csv(f_map, colClasses = c(code_insee = "character", num_maille = "integer"))

cat("SWI :", nrow(swi), "lignes | mailles :", n_distinct(swi$NUMERO),
    "| dates :", n_distinct(swi$DATE), "| periode :", min(swi$DATE), "->", max(swi$DATE), "\n")
cat("Doublons (maille, date) :", sum(duplicated(swi[, c("NUMERO", "DATE")])), "\n")
cat("Valeurs manquantes de SWI :", sum(is.na(swi$SWI_UNIF_MENS)), "\n")
print(summary(swi$SWI_UNIF_MENS))

cat("\nCorrespondance :", nrow(map), "lignes | communes :", n_distinct(map$code_insee),
    "| mailles :", n_distinct(map$num_maille), "\n")
cat("Codes INSEE de longueur 5 :", mean(nchar(map$code_insee) == 5), "\n")
s <- tapply(map$weight, map$code_insee, sum)
cat("Somme des poids par commune :\n"); print(summary(s))
cat("Communes dont la somme des poids s'ecarte de 1 de plus de 1e-6 :", sum(abs(s - 1) > 1e-6), "\n")
cat("Mailles de la correspondance absentes du SWI :",
    length(setdiff(unique(map$num_maille), unique(swi$NUMERO))), "\n")

#############################################
library(Matrix); library(readxl)

maille_ids <- sort(unique(swi$NUMERO)); dates <- sort(unique(swi$DATE))
S <- matrix(NA_real_, length(maille_ids), length(dates), dimnames = list(maille_ids, dates))
S[cbind(match(swi$NUMERO, maille_ids), match(swi$DATE, dates))] <- swi$SWI_UNIF_MENS
cat("Valeurs manquantes dans S :", sum(is.na(S)), "| dim :", dim(S), "\n")

com_ids <- sort(unique(map$code_insee))
W <- sparseMatrix(i = match(map$code_insee, com_ids), j = match(map$num_maille, maille_ids),
                  x = map$weight, dims = c(length(com_ids), length(maille_ids)),
                  dimnames = list(com_ids, maille_ids))
Wsum <- rowSums(W)

A <- as.matrix(W %*% S)     # version chapitre 2 : somme ponderee sans renormalisation
B <- A / Wsum               # version renormalisee
cat("dim A :", dim(A), "\n")

print(round(quantile(Wsum, c(0, .01, .05, .5, .95, .99, 1)), 3))
cat("Communes avec somme < 0.99 :", sum(Wsum < 0.99), "| < 0.90 :", sum(Wsum < 0.90), "\n")
print(round(head(sort(Wsum), 10), 3))
cat("Ecart |A - B| : max", round(max(abs(A - B)), 3), "| moyen", round(mean(abs(A - B)), 4), "\n")

# Controle sur une commune : les 3 mois de l'evenement 2003 de 01004
print(round(A["01004", c("200307", "200308", "200309")], 3))
print(round(B["01004", c("200307", "200308", "200309")], 3))
evt <- read_excel(file.path(root, "DATA/CatNat_drought_events_SWI.xlsx"))
print(evt %>% filter(code_insee == "01004") %>%
        select(`Début d'Evénement`, `Fin d'Evénement`, swi_min, month_swi_min) %>% head(6))

#########################################################
library(lubridate)
evt <- evt %>% mutate(insee = str_pad(code_insee, 5, pad = "0"),
                      d0 = format(floor_date(as.Date(`Début d'Evénement`), "month"), "%Y%m"),
                      d1 = format(floor_date(as.Date(`Fin d'Evénement`), "month"), "%Y%m"))
cat("Evenements :", nrow(evt), "| dates de debut hors 2000-2025 :", sum(evt$d0 < "200001" | evt$d1 > "202512"), "\n")

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
cat("Evenements sans reference :", sum(is.na(ref)), "| sans calcul A :", sum(is.na(evt$swi_A)), "\n")
for (v in c("swi_A", "swi_B")) {
  e <- abs(evt[[v]] - ref); e <- e[!is.na(e)]
  cat(v, ": ecart max", signif(max(e), 3), "| moyen", signif(mean(e), 3),
      "| part < 0.001 :", round(mean(e < 0.001), 4), "| part < 0.01 :", round(mean(e < 0.01), 4), "\n")
}
print(head(evt[which(abs(evt$swi_A - ref) > 0.01), c("insee", "d0", "d1", "swi_min", "swi_A", "swi_B")], 10))



dir_out <- "C:/Users/OumyDIONE/Documents/output/etape0"

# 1) SWI communal mensuel, poids renormalises (version B) : 34 741 communes x 312 mois
saveRDS(B, file.path(dir_out, "SWI_communal_mensuel.rds"))

# 2) Communes dont les poids commune-maille ne somment pas a 1 (a surveiller)
communes_poids_faibles <- data.frame(insee = com_ids, somme_poids = as.numeric(Wsum)) %>%
  filter(somme_poids < 0.99)
saveRDS(communes_poids_faibles, file.path(dir_out, "SWI_communes_poids_faibles.rds"))

cat("SWI :", round(file.size(file.path(dir_out, "SWI_communal_mensuel.rds")) / 1e6), "Mo |",
    "communes a surveiller :", nrow(communes_poids_faibles), "\n")
list.files(dir_out)

