# ============================================================
# 01_preparation_CatNat_ONRN_couts.R
# Preparation des bases ONRN + CatNat (secheresse) et modele de cout par arrete
# Sorties : base_evenements_couts.rds/.csv, serie_nationale_annuelle.rds/.csv
# Script NON teste sous R (pas d'interpreteur disponible) : relire une fois avant
# execution complete, et signaler tout message d'erreur complet en cas de souci.
# ============================================================
library(readxl); library(dplyr); library(tidyr); library(stringr); library(lubridate)
library(survival)

# ------------------------------------------------------------
# CONFIG (a adapter)
# ------------------------------------------------------------
f_cout       <- "C:/Users/OumyDIONE/Downloads/ONRN_CoutCumu_SECH_9521.xlsx"
f_cout_moyen <- "C:/Users/OumyDIONE/Downloads/ONRN_CoutMoyen_SECH_9521.xlsx"
f_freq       <- "C:/Users/OumyDIONE/Downloads/ONRN_Frequence_SECH_9521.xlsx"
f_sp         <- "C:/Users/OumyDIONE/Downloads/ONRN_SsurP_SECH_9521.xlsx"
f_catnat     <- "D:/travaux recherche/Catnat/traité par comunne/Arretés-CatNat-commune-1.xlsx"
f_carac      <- "C:/Users/OumyDIONE/Downloads/Drought_risk_review/DATA/Data_merging/CatNat_drought_final_analysis.xlsx"

dir_out <- "C:/Users/OumyDIONE/Documents/output/etape0"     # chemin complet, utilise pour toutes les sorties
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(f_cout), file.exists(f_cout_moyen), file.exists(f_freq),
          file.exists(f_sp), file.exists(f_catnat), file.exists(f_carac))

# ============================================================
# 1. LECTURE ONRN (sans la base CatNat)
# ============================================================
annee_min <- 1995; annee_max <- 2021

read_onrn <- function(file, nm) {
  v <- read_excel(file, sheet = 2, col_types = "text")
  names(v) <- c("insee", "commune", "classe")
  r <- read_excel(file, sheet = "Représentativité", col_types = "text")
  names(r) <- c("insee", "commune_r", "repres")
  v %>% select(insee, classe) %>%
    left_join(r %>% select(insee, repres), by = "insee") %>%
    mutate(insee = str_pad(insee, 5, pad = "0")) %>%
    rename(!!nm := classe)
}

# Correction : le fichier "cout moyen" a ses colonnes commune/classe triees en ordre
# inverse (code croissant, valeurs decroissantes). On les remet dans le bon ordre.
read_onrn_moy <- function(file) {
  v <- read_excel(file, sheet = 2, col_types = "text")
  names(v) <- c("insee", "commune", "classe")
  v <- v %>% mutate(insee = str_pad(insee, 5, pad = "0")) %>% arrange(insee)
  tibble(insee = v$insee, commune_moy = rev(v$commune), moy_classe = rev(v$classe))
}

cout <- read_onrn(f_cout, "cout_classe")
freq <- read_onrn(f_freq, "freq_classe") %>% select(-repres)
sp   <- read_onrn(f_sp,   "sp_classe")   %>% select(-repres)

nm_cout <- read_excel(f_cout, sheet = 2, col_types = "text")
names(nm_cout) <- c("insee", "commune", "x")
nm_cout <- nm_cout %>% mutate(insee = str_pad(insee, 5, pad = "0")) %>% select(insee, commune)

moy <- read_onrn_moy(f_cout_moyen)
chk <- moy %>% left_join(nm_cout, by = "insee")
cat("Alignement cout moyen apres correction :",
    round(mean(chk$commune_moy == chk$commune, na.rm = TRUE), 4), "(attendu : 1)\n")

onrn <- cout %>%
  left_join(freq, by = "insee") %>%
  left_join(sp, by = "insee") %>%
  left_join(moy %>% select(insee, moy_classe), by = "insee") %>%
  filter(!str_starts(insee, "97"), !str_starts(insee, "98"))

cat("Communes ONRN :", nrow(onrn), "| doublons :", sum(duplicated(onrn$insee)), "(attendu : 34839 / 0)\n")
print(table(onrn$cout_classe, useNA = "ifany"))

# ============================================================
# 2. LECTURE ET NETTOYAGE DE LA BASE CATNAT
# ============================================================
cn <- read_excel(f_catnat)
names(cn) <- c("insee", "commune", "debut", "fin", "date_arrete", "date_jo",
               "peril", "code_peril", "franchise", "nor", "decision")

# Certaines dates sont stockees comme numeros de serie Excel (253 lignes) plutot
# que comme texte jj/mm/aaaa : parse_date_any() gere les deux formats.
parse_date_any <- function(x) {
  x   <- as.character(x)
  out <- suppressWarnings(dmy(x))
  num <- suppressWarnings(as.numeric(x))
  idx <- is.na(out) & !is.na(num)
  out[idx] <- as.Date(num[idx], origin = "1899-12-30")
  out
}

cn <- cn %>% mutate(
  insee       = str_pad(insee, 5, pad = "0"),
  debut       = parse_date_any(debut), fin = parse_date_any(fin),
  date_arrete = parse_date_any(date_arrete), date_jo = parse_date_any(date_jo))
cat("Dates NA apres parse_date_any : debut", sum(is.na(cn$debut)), "| fin", sum(is.na(cn$fin)),
    "| arrete", sum(is.na(cn$date_arrete)), "| JO", sum(is.na(cn$date_jo)), "(attendu : 0 partout)\n")

sech <- cn %>% filter(code_peril == "SEC") %>%
  mutate(reconnu = as.integer(str_starts(decision, "Reconnue")))
cat("Lignes SEC :", nrow(sech), "\n")
print(table(sech$decision, useNA = "ifany"))

# Deduplication : derniere decision par (commune, debut, fin), triee par date d'arrete/JO
sech_final <- sech %>% arrange(insee, debut, fin, date_arrete, date_jo) %>%
  group_by(insee, debut, fin) %>% slice_tail(n = 1) %>% ungroup()
cat("Lignes SEC :", nrow(sech), "| apres deduplication :", nrow(sech_final), "(attendu : 75488)\n")
print(table(sech_final$reconnu))

# ============================================================
# 3. VERIFICATION DE L'UNITE DE COMPTAGE (recoupement avec la fiche ONRN)
# ============================================================
cible <- c("0" = 22496, "1-3" = 9091, "4-6" = 2265, "7-10" = 836, "11-15" = 140, ">15" = 11)
date_coupure_jo <- as.Date("2024-10-09")

comp <- function(cle, coupure = TRUE) {
  df <- sech_final %>% filter(reconnu == 1)
  if (coupure) df <- df %>% filter(date_jo <= date_coupure_jo)
  n <- df %>% distinct(insee, across(all_of(cle))) %>% count(insee, name = "n_all")
  x <- onrn %>% select(insee) %>% left_join(n, by = "insee") %>%
    mutate(n_all = coalesce(n_all, 0L),
           classe = cut(n_all, c(-1, 0, 3, 6, 10, 15, Inf), labels = names(cible)))
  as.integer(table(x$classe))
}
res_unite <- rbind(cible = cible,
                   evenement         = comp(c("debut", "fin"), coupure = FALSE),
                   evenement_coupure = comp(c("debut", "fin")),
                   par_date_arrete   = comp("date_arrete"),
                   par_nor           = comp("nor"))
colnames(res_unite) <- names(cible)
print(res_unite)   # l'unite retenue est "nor" (arrete), la plus proche de la fiche ONRN sur les classes hautes

# ============================================================
# 4. DIAGNOSTIC DE COHERENCE ONRN / CATNAT (groupes A/B/C/D)
# ============================================================
onrn <- onrn %>% mutate(pas_sinistre = str_starts(cout_classe, "Pas de"))

# Bornes des classes de cout ONRN, avec un point milieu par classe
top_hi <- 3e8
bornes <- tribble(
  ~cout_classe,               ~lo,  ~hi,
  "Entre 0 k€ et 100 k€",      0,    1e5,
  "Entre 100 k€ et 500 k€",    1e5,  5e5,
  "Entre 500 k€ et 2 M€",      5e5,  2e6,
  "Entre 2 M€ et 5 M€",        2e6,  5e6,
  "Entre 5 M€ et 10 M€",       5e6,  1e7,
  "Entre 10 M€ et 50 M€",      1e7,  5e7,
  "Entre 50 M€ et 100 M€",     5e7,  1e8,
  "Supérieur à 100 M€",        1e8,  top_hi
) %>% mutate(mid = (lo + hi) / 2)

onrn <- onrn %>% select(-any_of(c("lo", "hi", "mid"))) %>% left_join(bornes, by = "cout_classe")

# ev : une ligne par (commune, arrete) reconnu, fenetre 1995-2021, JO <= coupure
ev <- sech_final %>%
  filter(reconnu == 1, year(debut) >= annee_min, year(debut) <= annee_max,
         date_jo <= date_coupure_jo) %>%
  arrange(insee, debut) %>%
  distinct(insee, nor, .keep_all = TRUE) %>%
  mutate(annee = year(debut))

nrec <- ev %>% count(insee, name = "N_rec")
d <- onrn %>% left_join(nrec, by = "insee") %>% mutate(N_rec = coalesce(N_rec, 0L))

diag <- d %>% mutate(groupe = case_when(
  !pas_sinistre & N_rec == 0 ~ "A. cout > 0 mais aucune reconnaissance 1995-2021",
  pas_sinistre  & N_rec > 0  ~ "B. reconnue(s) mais aucun sinistre dans la base CCR",
  !pas_sinistre & N_rec > 0  ~ "C. cout > 0 et reconnue(s)",
  TRUE                       ~ "D. ni cout ni reconnaissance")) %>% count(groupe)
print(diag)   # attendu : A=241, B=613, C=11130, D=22855

# ============================================================
# 5. PONT ANNUEL (allocation du cout cumule aux annees d'evenement)
# ============================================================
v_pont <- d %>% filter(!pas_sinistre, N_rec > 0) %>%
  transmute(insee, v_lo = lo / N_rec, v_mid = mid / N_rec, v_hi = hi / N_rec)

pred_annee <- ev %>% inner_join(v_pont, by = "insee") %>%
  group_by(annee) %>%
  summarise(C_lo = sum(v_lo) / 1e9, C_mid = sum(v_mid) / 1e9, C_hi = sum(v_hi) / 1e9,
            n_evenements = n())
print(pred_annee, n = 40)
cat("Total Md EUR bas/milieu/haut :",
    round(c(sum(pred_annee$C_lo), sum(pred_annee$C_mid), sum(pred_annee$C_hi)), 2), "\n")

# Sensibilite : allocation a parts egales (v_pont) vs allocation ponderee par
# l'ampleur nationale de l'annee (v_alt)
nat <- ev %>% count(annee, name = "n_nat")
ev_alt <- ev %>% left_join(nat, by = "annee") %>%
  inner_join(d %>% select(insee, mid), by = "insee") %>%
  filter(!is.na(mid)) %>%
  group_by(insee) %>% mutate(w_alt = n_nat / sum(n_nat), v_alt = mid * w_alt) %>% ungroup()

alt <- ev_alt %>% group_by(annee) %>% summarise(C_alt = sum(v_alt) / 1e9)
cmp_alloc <- pred_annee %>% select(annee, C_mid) %>% left_join(alt, by = "annee")
cat("Spearman entre les deux allocations :",
    round(cor(cmp_alloc$C_mid, cmp_alloc$C_alt, method = "spearman"), 3), "(attendu ~0.988)\n")

# ============================================================
# 6. TABLE COMMUNALE (logements, argile, chapitre 2) ET MODELE DE COUT
# ============================================================
carac <- read_excel(f_carac) %>% mutate(insee = str_pad(as.character(code_insee), 5, pad = "0"))
cat("Table communale :", nrow(carac), "lignes |", n_distinct(carac$insee), "communes\n")

comm <- carac %>% group_by(insee) %>%
  summarise(n_dwellings = first(n_dwellings), CSS = first(CSS),
            old_share = first(old_housing_share), houses = first(share_individual_houses),
            n_evt_ch2 = n()) %>%
  mutate(CSS100 = pmin(CSS, 100) / 100, old_s = old_share / 100, houses_s = houses / 100)

est <- d %>% inner_join(comm, by = "insee") %>% filter(N_rec > 0)
cat("Communes reconnues 1995-2021 dans la table :", nrow(est),
    "| avec cout :", sum(!est$pas_sinistre), "| sans cout :", sum(est$pas_sinistre), "\n")

pos <- est %>% filter(!pas_sinistre) %>%
  mutate(lo_s = ifelse(lo == 0, NA, lo),                                     # censure a gauche
         hi_s = ifelse(cout_classe == "Supérieur à 100 M€", NA, hi))         # censure a droite

# Validation croisee (5 blocs) : modele retenu = "sans_old" (sans old_share, cf. fil de travail)
set.seed(1)
pos$fold <- sample(1:5, nrow(pos), replace = TRUE)
S <- function(rhs) as.formula(paste("Surv(lo_s, hi_s, type = 'interval2') ~", rhs))
forms <- list(
  nul        = S("1"),
  N_rec_seul = S("log(N_rec)"),
  initial    = S("log(n_dwellings) + CSS100 + old_s + houses_s + log(N_rec)"),
  sans_old   = S("log(n_dwellings) + CSS100 + houses_s + log(N_rec)")
)
cv_score <- function(f, data, foldvar = "fold") {
  s <- numeric(nrow(data))
  for (k in 1:5) {
    tr <- data[data[[foldvar]] != k, ]; te <- data[data[[foldvar]] == k, ]
    m  <- survreg(f, data = tr, dist = "lognormal")
    mu <- predict(m, newdata = te, type = "lp"); sg <- m$scale
    plo <- ifelse(is.na(te$lo_s), 0, pnorm((log(te$lo_s) - mu) / sg))
    phi <- ifelse(is.na(te$hi_s), 1, pnorm((log(te$hi_s) - mu) / sg))
    s[data[[foldvar]] == k] <- log(pmax(phi - plo, 1e-12))
  }
  mean(s)
}
print(round(sapply(forms, cv_score, data = pos), 4))   # sans_old attendu proche de -1.12

# Modeles finaux retenus
m_p_f <- glm(!pas_sinistre ~ log(n_dwellings) + CSS100 + houses_s + log(N_rec),
             data = est, family = binomial)
m_c_f <- survreg(forms$sans_old, data = pos, dist = "lognormal")
cat("Sigma (dispersion log) du modele retenu :", round(m_c_f$scale, 3), "\n")

# ============================================================
# 7. BASE FINALE : COUT PREDIT PAR COMMUNE, PUIS PAR EVENEMENT
# ============================================================
comm_pred <- comm %>% mutate(N_rec = 1)
comm_pred$cout_1_arrete_modele <-
  predict(m_p_f, newdata = comm_pred, type = "response") *
  exp(predict(m_c_f, newdata = comm_pred, type = "lp"))       # mediane, EUR 2021

commune_tab <- d %>%
  left_join(v_pont, by = "insee") %>%
  left_join(comm %>% select(insee, n_dwellings, CSS, houses), by = "insee") %>%
  left_join(comm_pred %>% select(insee, cout_1_arrete_modele), by = "insee") %>%
  transmute(insee, cout_classe, cout_lo = lo, cout_hi = hi, cout_mid = mid,
            N_rec, v_lo, v_mid, v_hi, n_dwellings, CSS, houses, cout_1_arrete_modele)

base_evt <- sech_final %>%
  mutate(annee = year(debut)) %>%
  select(insee, nor, debut, fin, date_arrete, date_jo, annee, decision, reconnu) %>%
  left_join(commune_tab, by = "insee")

stopifnot(nrow(base_evt) == nrow(sech_final))
cat("Base finale :", nrow(base_evt), "lignes (attendu 75488)\n")
cat("Lignes avec cout ONRN :", sum(!is.na(base_evt$cout_classe) & !str_starts(base_evt$cout_classe, "Pas de")),
    "| avec cout modele :", sum(!is.na(base_evt$cout_1_arrete_modele)), "\n")

# NOTE : fermer le fichier .csv dans Excel avant de relancer cette ligne, sinon
# write.csv2 echoue avec "Permission denied" (le .rds, lui, s'ecrit toujours).
saveRDS(base_evt, file.path(dir_out, "base_evenements_couts.rds"))
write.csv2(base_evt, file.path(dir_out, "base_evenements_couts.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ============================================================
# 8. SERIE NATIONALE ANNUELLE (part de demandes, part reconnue, ponderees par le cout)
# ============================================================
w <- base_evt %>% filter(!is.na(cout_1_arrete_modele)) %>% distinct(insee, cout_1_arrete_modele)
N_com <- nrow(w)                             # nombre de communes ponderees (scalaire)
W_tot <- sum(w$cout_1_arrete_modele)         # poids total (scalaire)
cat("Communes ponderees :", N_com, "| poids total :", round(W_tot / 1e9, 2), "Md EUR\n")

yr <- base_evt %>%
  filter(!is.na(cout_1_arrete_modele), annee >= 2000, annee <= 2023) %>%
  group_by(insee, annee) %>%
  summarise(R = max(reconnu), poids = first(cout_1_arrete_modele), .groups = "drop")

serie <- yr %>% group_by(annee) %>%
  summarise(communes_demande = n(), communes_reconnues = sum(R),
            a_simple = n() / N_com,      rho_simple = sum(R) / N_com,
            a_pond   = sum(poids) / W_tot, rho_pond = sum(poids * R) / W_tot) %>%
  mutate(r_simple = rho_simple / a_simple, r_pond = rho_pond / a_pond) %>%
  mutate(across(c(a_simple, rho_simple, a_pond, rho_pond, r_simple, r_pond), ~ round(.x, 3)))
print(serie, n = 30)

saveRDS(serie, file.path(dir_out, "serie_nationale_annuelle.rds"))
write.csv2(serie, file.path(dir_out, "serie_nationale_annuelle.csv"), row.names = FALSE)

# ============================================================
# 9. PREMIER REGARD : SWI DES EVENEMENTS DECLARES (indicatif, avant l'etape SWI complete)
# ============================================================
swi_an <- carac %>% mutate(annee = year(event_start)) %>%
  group_by(annee) %>%
  summarise(swi_moy = mean(SWI, na.rm = TRUE), swi_p10 = quantile(SWI, 0.10, na.rm = TRUE), n = n())

d3 <- serie %>% filter(annee >= 2002) %>% left_join(swi_an, by = "annee")
print(d3 %>% select(annee, a_pond, rho_pond, r_pond, swi_moy, swi_p10, n), n = 30)
print(round(cor(d3[, c("a_pond", "rho_pond", "r_pond", "swi_moy", "swi_p10")],
                method = "spearman", use = "complete.obs"), 2))

# Modeles indicatifs sur ce SWI provisoire (le SWI final, sur toutes les communes,
# est construit dans 02_construction_alea_SWI.R -- ceci est un premier regard seulement)
d3 <- d3 %>% mutate(regime_provisoire = as.integer(annee >= 2020))   # hypothese de rupture, non retenue telle quelle

m_a <- glm(a_pond ~ swi_moy, data = d3, family = quasibinomial)
m_r <- glm(r_pond ~ swi_p10, data = d3, family = quasibinomial)
m_r2 <- glm(r_pond ~ swi_p10 + regime_provisoire, data = d3, family = quasibinomial)
print(round(summary(m_a)$coefficients, 3))
print(round(summary(m_r)$coefficients, 3))
print(round(summary(m_r2)$coefficients, 3))

# ============================================================
# 10. SCENARIO D'INCERTITUDE PROVISOIRE SUR r (loi Beta, 22 annees)
#     Remplace par le modele a deux regimes dans 03_reconnaissance_regime.R
# ============================================================
r_vec <- d3$r_pond
m_b <- mean(r_vec); v_b <- var(r_vec)
k_beta <- m_b * (1 - m_b) / v_b - 1
cat("r : moyenne", round(m_b, 3), "| ecart-type", round(sd(r_vec), 3),
    "| Beta(", round(m_b * k_beta, 2), ",", round((1 - m_b) * k_beta, 2), ")\n")
cat("Autocorrelation d'ordre 1 de r :", round(cor(r_vec[-1], r_vec[-length(r_vec)]), 2), "\n")
cat("Correlation de rang r vs volume (a_pond) :", round(cor(d3$r_pond, d3$a_pond, method = "spearman"), 2), "\n")

# Independance a et r : la queue de rho = a*r est-elle affectee ?
set.seed(10); n_sim <- 1e5
a_obs <- d3$a_pond; r_obs <- d3$r_pond
i_a <- sample(seq_along(a_obs), n_sim, replace = TRUE)
rho_pair <- a_obs[i_a] * r_obs[i_a]                                          # couples observes
rho_shuf <- a_obs[i_a] * r_obs[sample(seq_along(r_obs), n_sim, replace = TRUE)]  # independance empirique
rho_beta <- a_obs[i_a] * rbeta(n_sim, m_b * k_beta, (1 - m_b) * k_beta)      # independance, r ~ Beta

qs <- c(0.50, 0.90, 0.95, 0.99)
tab_dep <- rbind(couples_observes  = c(quantile(rho_pair, qs), moyenne = mean(rho_pair)),
                 independance_emp  = c(quantile(rho_shuf, qs), moyenne = mean(rho_shuf)),
                 independance_beta = c(quantile(rho_beta, qs), moyenne = mean(rho_beta)))
print(round(tab_dep, 3))
cat("rho observe, maximum :", round(max(a_obs * r_obs), 3), "\n")

cat("\nTermine. Sorties dans", dir_out, "\n")
cat("Suite du pipeline : 02_construction_alea_SWI.R (SWI communal, indice d'alea)\n")

################################################ Accessoire 
# Partie 1 (logit) et partie 2 (lognormal censure)
cat("n communes reconnues (partie 1) :", nrow(est), "| n avec sinistre (partie 2) :", nrow(pos), "\n")
print(round(coef(summary(m_p_f)), 3))
cat("logLik partie 1 :", round(as.numeric(logLik(m_p_f))), "| AIC :", round(AIC(m_p_f)), "\n")
print(summary(m_c_f))
cat("sigma :", round(m_c_f$scale, 3), "| logLik partie 2 :", round(m_c_f$loglik[2]), "\n")
# Validation croisee des 4 specifications
print(round(sapply(forms, cv_score, data = pos), 4))

