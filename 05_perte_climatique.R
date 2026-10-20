# ============================================================
# 05_perte_climatique.R
# Annual climate loss (drought -> recognition -> cost),
# calibrated on the ONRN cost bridge
# Depends only on the files saved in output/etape0
# Outputs: calibrage_perte_climatique.rds, scenarios_perte_climatique.rds
# ============================================================
library(dplyr); library(tidyr); library(stringr); library(lubridate)

dir_out <- file.path("output", "etape0")
base  <- readRDS(file.path(dir_out, "base_evenements_couts.rds"))
Zj    <- readRDS(file.path(dir_out, "anomalies_communales_jul_dec_2000_2025.rds"))
w     <- readRDS(file.path(dir_out, "poids_cout_communes.rds"))
mod   <- readRDS(file.path(dir_out, "modeles_communaux_depot_reconnaissance.rds"))
panel <- readRDS(file.path(dir_out, "panel_commune_annee_alea.rds"))
cat("Zj:", dim(Zj), "| weighted communes:", nrow(w), "| panel:", nrow(panel), "\n")

# ------------------------------------------------------------
# 1. Pont ONRN : cout annuel C_mid (Md EUR), arretes reconnus 1995-2021
# ------------------------------------------------------------
ev <- base %>%
  filter(reconnu == 1, annee >= 1995, annee <= 2021, date_jo <= as.Date("2024-10-09")) %>%
  arrange(insee, debut) %>%
  distinct(insee, nor, .keep_all = TRUE)
pont <- ev %>% group_by(annee) %>% summarise(C_mid = sum(v_mid, na.rm = TRUE) / 1e9)
print(round(pont %>% filter(annee %in% c(2003, 2011, 2018, 2019, 2020)), 3))

# ------------------------------------------------------------
# 2. Poids de cout, cout observe avec poids (C_obs) et cout predit (C_pred)
#    C_pred = Omega * somme_i omega_i * p^A_it * p^R_it   (Md EUR)
# ------------------------------------------------------------
omega <- setNames(w$omega, w$insee)
Vtot  <- sum(w$cout_1_arrete_modele)            # Omega, en EUR
panel$R0 <- coalesce(panel$R, 0L)               # R est vide quand la commune n'a pas demande
panel$om <- omega[panel$insee]

# (1) cout "observe" avec les poids : reconnaissances observees x cout attendu d'un arrete
C_obs <- tapply(Vtot * panel$om * panel$A * panel$R0, panel$annee, sum) / 1e9

# (2) cout predit par le modele communal (demande x reconnaissance)
p_pred <- predict(mod$m_dep, newdata = panel, type = "response") *
  predict(mod$r_rec, newdata = panel, type = "response")
C_pred <- tapply(Vtot * panel$om * p_pred, panel$annee, sum) / 1e9

cmp <- data.frame(annee = as.integer(names(C_obs)), C_obs = C_obs, C_pred = C_pred) %>%
  left_join(pont, by = "annee") %>%
  mutate(facteur_obs  = C_mid / C_obs,
         facteur_pred = C_mid / C_pred)
print(round(cmp %>% filter(annee <= 2021), 3), row.names = FALSE)

s <- cmp %>% filter(annee <= 2021)
cat("Total 2002-2021 (Md EUR) : pont", round(sum(s$C_mid), 2),
    "| observe avec poids", round(sum(s$C_obs), 2),
    "| predit", round(sum(s$C_pred), 2), "\n")
cat("Correlation de rang avec le pont : observe", round(cor(s$C_mid, s$C_obs, method = "spearman"), 2),
    "| predit", round(cor(s$C_mid, s$C_pred, method = "spearman"), 2), "\n")

# ------------------------------------------------------------
# 3. Facteur de severite lie a l'alea (exploration)
# ------------------------------------------------------------
s <- cmp %>% filter(annee >= 2002, annee <= 2021) %>%
  left_join(readRDS(file.path(dir_out, "H_indice_alea_national.rds")) %>% select(annee, H_pond),
            by = "annee")
s$lm <- log(s$facteur_obs)
print(round(cor(s$lm, s$H_pond, method = "spearman"), 2))
print(round(cor(s$lm, log(s$C_obs), method = "spearman"), 2))

m_sev <- lm(lm ~ log(C_obs), data = s)
print(round(summary(m_sev)$coefficients, 3))
cat("R2 :", round(summary(m_sev)$r.squared, 2), "| sigma :", round(summary(m_sev)$sigma, 2), "\n")

m_sev2 <- lm(lm ~ H_pond, data = s)
print(round(summary(m_sev2)$coefficients, 3))
cat("R2 :", round(summary(m_sev2)$r.squared, 2), "| sigma :", round(summary(m_sev2)$sigma, 2), "\n")

# Perte reconstruite avec la severite modelisee : quelle erreur ?
s$C_rec <- s$C_obs * exp(fitted(m_sev))
cat("Total pont / reconstruit :", round(sum(s$C_mid), 2), round(sum(s$C_rec), 2),
    "| RMSE :", round(sqrt(mean((s$C_mid - s$C_rec)^2)), 3),
    "| ecart-type du pont :", round(sd(s$C_mid), 3), "\n")

# ------------------------------------------------------------
# 4. Facteur constant vs relation log-lineaire, sur C_pred (validation leave-one-out)
# ------------------------------------------------------------
s2 <- cmp %>% filter(annee >= 2002, annee <= 2021)
lp <- log(s2$C_mid); x <- log(s2$C_pred); xo <- log(s2$C_obs)

m_full <- lm(lp ~ x)
print(round(summary(m_full)$coefficients, 3))
cat("R2 :", round(summary(m_full)$r.squared, 2), "| sigma (log) :", round(summary(m_full)$sigma, 2), "\n")

loo_rmse <- function(pred_fun) sqrt(mean(sapply(seq_len(nrow(s2)), function(i) {
  (s2$C_mid[i] - pred_fun(i))^2 })))
p_const_obs  <- function(i) exp(mean(lp[-i] - xo[-i]) + xo[i])
p_slope_obs  <- function(i) { m <- lm(lp[-i] ~ xo[-i]); exp(coef(m)[1] + coef(m)[2] * xo[i]) }
p_const_pred <- function(i) exp(mean(lp[-i] - x[-i]) + x[i])
p_slope_pred <- function(i) { m <- lm(lp[-i] ~ x[-i]); exp(coef(m)[1] + coef(m)[2] * x[i]) }
print(round(c(ecart_type_pont = sd(s2$C_mid),
              obs_constant  = loo_rmse(p_const_obs),  obs_pente  = loo_rmse(p_slope_obs),
              pred_constant = loo_rmse(p_const_pred), pred_pente = loo_rmse(p_slope_pred)), 3))

# ------------------------------------------------------------
# 5. Calibration : facteur d'echelle k (moyenne geometrique) et dispersion sigma
# ------------------------------------------------------------
k_geo <- exp(mean(lp - x)); sig <- sd(lp - x)
cat("k =", round(k_geo, 2), "| sigma (log) =", round(sig, 2), "| n =", length(lp), "\n")
saveRDS(list(k = k_geo, sigma = sig, n = length(lp), Vtot = Vtot,
             note = "L_cat = k * C_pred * exp(eps), eps ~ N(0, sigma^2) ; calibre sur le pont ONRN 2002-2021 (Md EUR, euros 2021)"),
        file.path(dir_out, "calibrage_perte_climatique.rds"))

# ------------------------------------------------------------
# 6. Simulation de la perte de marche selon (delta, regime)
#    delta : decalage additif de z_it ; regime : P_t = 0 (avant) ou 1 (apres 2017)
# ------------------------------------------------------------
set.seed(42)
yrs_obs <- sort(unique(panel$annee))

cp_years <- function(delta, regime) sapply(yrs_obs, function(y) {
  nd <- panel[panel$annee == y, ]
  nd$z <- nd$z + delta
  nd$apres17 <- (regime == "apres")
  p <- predict(mod$m_dep, newdata = nd, type = "response") *
    predict(mod$r_rec, newdata = nd, type = "response")
  Vtot * sum(nd$om * p) / 1e9
})

sim_C <- function(delta, regime, n_sim = 2000, sig_eps = sig) {
  Cp <- cp_years(delta, regime)
  draw <- sample(Cp, n_sim, replace = TRUE)
  k_geo * draw * exp(rnorm(n_sim, 0, sig_eps))
}
res <- do.call(rbind, lapply(c(0, 0.25, 0.5), function(dl) do.call(rbind, lapply(c("avant", "apres"), function(rg) {
  L <- sim_C(dl, rg)
  data.frame(delta = dl, regime = rg, moyenne = mean(L), q50 = quantile(L, .5), q95 = quantile(L, .95),
             q995 = quantile(L, .995), max = max(L))
}))))
res_p <- res %>% mutate(across(where(is.numeric), ~ round(.x, 2)))
rownames(res_p) <- NULL
print(res_p)

# ------------------------------------------------------------
# 7. Sensibilite au bruit : aucun / normal / residus empiriques (reference) / normal ecrete a 1,5 sigma
# ------------------------------------------------------------
Cp <- list(avant = cp_years(0, "avant"), apres = cp_years(0, "apres"))
eps_obs <- (lp - x) - log(k_geo)                 # 20 residus observes, moyenne 0 en log
cat("Residus observes : min", round(exp(min(eps_obs)), 2), "| max", round(exp(max(eps_obs)), 2), "(multiplicateurs)\n")

sim <- function(cp, type, n = 20000) {
  d <- sample(cp, n, replace = TRUE)
  e <- switch(type,
              aucun     = rep(0, n),
              lognormal = rnorm(n, 0, sig),
              residus   = sample(eps_obs, n, TRUE),
              tronque   = pmax(pmin(rnorm(n, 0, sig), 1.5 * sig), -1.5 * sig))   # ecretage a +/- 1,5 sigma
  k_geo * d * exp(e)
}
set.seed(1)
tab <- do.call(rbind, lapply(names(Cp), function(rg) do.call(rbind, lapply(c("aucun", "lognormal", "residus", "tronque"), function(tp) {
  L <- sim(Cp[[rg]], tp)
  data.frame(regime = rg, bruit = tp, moyenne = mean(L), q50 = median(L), q95 = quantile(L, .95),
             q995 = quantile(L, .995), max = max(L))
}))))
print(tab %>% mutate(across(where(is.numeric), ~ round(.x, 2))), row.names = FALSE)

# ------------------------------------------------------------
# 8. Grille complete des scenarios et enregistrement
# ------------------------------------------------------------
grille <- expand.grid(delta = c(0, 0.25, 0.5), regime = c("avant", "apres"), stringsAsFactors = FALSE)
Cp_grille <- lapply(seq_len(nrow(grille)), function(i) cp_years(grille$delta[i], grille$regime[i]))
names(Cp_grille) <- paste0("d", grille$delta, "_", grille$regime)

saveRDS(list(Cp = Cp_grille, annees_alea = yrs_obs, k = k_geo, sigma = sig, eps_obs = eps_obs,
             note = "Cp = perte sans bruit (Md EUR, euros 2021) par annee d'alea observee ; L_cat = k * Cp * exp(eps) ; eps tire dans eps_obs (reference) ou N(0,sigma) ecrete a 1.5 sigma (sensibilite)"),
        file.path(dir_out, "scenarios_perte_climatique.rds"))
print(round(sapply(Cp_grille, function(v) c(moyenne = mean(v), max = max(v))), 2))
