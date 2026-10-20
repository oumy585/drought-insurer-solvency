# ============================================================
# 07_assemblage.R
# Total loss = financial + climate; VaR/ES, grid of configurations,
# tables and figures.
# Depends only on files saved in output/macro and output/etape0.
# Run from the repository root (see README for the data/ layout).
# ============================================================
library(dplyr); library(ggplot2)

dir_macro <- file.path("data", "macro")
dir_out   <- file.path("output", "etape0")

fin   <- readRDS(file.path(dir_macro, "simulation_financiere_reference.rds"))
clim  <- readRDS(file.path(dir_out, "scenarios_perte_climatique.rds"))
mod   <- readRDS(file.path(dir_out, "modeles_communaux_depot_reconnaissance.rds"))
panel <- readRDS(file.path(dir_out, "panel_commune_annee_alea.rds"))
cat("Module financier :", fin$config_sim, "\n")
cat("  n =", length(fin$L_net), "| moyenne =", round(mean(fin$L_net), 1), "| sd =", round(sd(fin$L_net), 1), "\n")
cat("\nModule climatique : scenarios =", names(clim$Cp), "\n")
cat("  annees d'alea :", length(clim$annees_alea), "| k =", round(clim$k, 2), "| sigma =", round(clim$sigma, 2), "\n")

qs <- c(.95, .99, .995)
es <- function(L, p) mean(L[L >= quantile(L, p)])

# Part climatique dans une bande etroite autour du quantile p (VaR), plutot
# qu'un seul tirage ponctuel, pour reduire le bruit Monte Carlo.
part_clim_at_VaR <- function(Lc, Lt, p = .995, window = 0.001) {
  lo <- quantile(Lt, p - window); hi <- quantile(Lt, p + window)
  idx <- Lt >= lo & Lt <= hi
  mean(Lc[idx]) / mean(Lt[idx])
}

# ------------------------------------------------------------
# Theme journal (defini une seule fois, utilise par toutes les figures)
# ------------------------------------------------------------
theme_journal <- theme_minimal(base_size = 11, base_family = "serif") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0, margin = margin(b = 6)),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    plot.caption = element_text(size = 8, color = "grey35", hjust = 0, margin = margin(t = 8)),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8.5, color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95", color = NA)
  )

# ============================================================
# 1. EXPLORATION D'ECHELLE (diagnostic, non retenu comme resultat)
#    Montre pourquoi un Q < 1 sans reference de marche ne suffit pas :
#    le portefeuille financier stylise n'a pas d'echelle monetaire propre.
# ============================================================
Cp_ref <- clim$Cp$d0_apres
moy_clim_brut <- mean(clim$k * Cp_ref)   # Md EUR, avant bruit, echelle nationale
cat("\n[diagnostic] Perte climatique nationale moyenne (apres/delta=0, avant bruit) :",
    round(moy_clim_brut, 3), "Md EUR\n")

ratios_cibles <- c(0.1, 0.25, 0.5, 1, 2)
Q_grille_diag <- ratios_cibles * abs(mean(fin$L_net)) / moy_clim_brut
cat("[diagnostic] Q necessaire pour atteindre ces ratios sans ancrage de taille (irrealiste, Q>1) :\n")
print(data.frame(ratio_cible = ratios_cibles, Q = round(Q_grille_diag, 5)))

taille_actifs_grille <- c(5, 10, 20, 50)          # Md EUR, actifs totaux de l'assureur stylise
lambda_grille <- taille_actifs_grille / 800       # 800 = E+Ms+Ml = 100+300+400 unites du portefeuille stylise

res_echelle <- sapply(lambda_grille, function(l) {
  Lf_eur <- l * fin$L_net
  c(moyenne = mean(Lf_eur), sd = sd(Lf_eur), q995 = quantile(Lf_eur, .995))
})
colnames(res_echelle) <- paste0(taille_actifs_grille, "MdEUR")
cat("\n[diagnostic] Perte financiere selon la taille d'actifs :\n")
print(round(res_echelle, 3))

Q_marche <- c(0.01, 0.02, 0.05, 0.10)
res_ratio <- sapply(taille_actifs_grille, function(ta) {
  Lf_moy <- mean((ta / 800) * fin$L_net)
  (Q_marche * moy_clim_brut) / abs(Lf_moy)
})
rownames(res_ratio) <- paste0("Q=", Q_marche); colnames(res_ratio) <- paste0(taille_actifs_grille, "MdEUR")
cat("\n[diagnostic] Ratio perte climatique / perte financiere (moyennes), Q x taille :\n")
print(round(res_ratio, 3))

# ============================================================
# 2. ANCRAGE DE DELTA SUR UNE SOURCE PUBLIEE (Covea, 2022 ; RCP 8.5, horizon 2050 :
#    +70% de frequence d'eligibilite aux arretes secheresse)
#    Objectif : ne pas presenter delta comme un choix arbitraire.
# ============================================================
base_dep <- mean(predict(mod$m_dep, newdata = panel, type = "response"))

f_hausse <- function(delta) {
  nd <- panel
  nd$z <- nd$z + delta
  mean(predict(mod$m_dep, newdata = nd, type = "response")) / base_dep - 1
}

grille_delta_test <- seq(0, 1.5, by = 0.05)
hausse_pct <- sapply(grille_delta_test, f_hausse)
cat("\n[ancrage] Hausse de la probabilite moyenne de depot selon delta :\n")
print(data.frame(delta = grille_delta_test, hausse_pct = round(hausse_pct * 100, 1)))

delta_covea <- tryCatch(
  uniroot(function(d) f_hausse(d) - 0.70, c(0, 1.5))$root,
  error = function(e) { cat("uniroot a echoue : la hausse cible (70%) n'est peut-etre pas atteinte dans [0, 1.5] ; regarder la table ci-dessus\n"); NA }
)
cat("\n[ancrage] Delta correspondant a +70% de hausse du depot (Covea, RCP 8.5, horizon 2050) :",
    round(delta_covea, 3), "\n")
cat("[ancrage] Comparaison : delta=0.5 produit une hausse de", round(f_hausse(0.5) * 100, 1), "% du depot moyen\n")

# ============================================================
# 3. SIMULATION DE REFERENCE : Q = 0.05, actifs = 20 MdEUR, scenario apres/delta=0
# ============================================================
Q_ref <- 0.05
taille_ref <- 20
lambda_ref <- taille_ref / 800

n_sim <- 10000
set.seed(5000)

Lf_eur <- lambda_ref * sample(fin$L_net, n_sim, replace = TRUE)
Cp_tirage  <- sample(clim$Cp$d0_apres, n_sim, replace = TRUE)
eps_tirage <- sample(clim$eps_obs, n_sim, replace = TRUE)
Lc_marche  <- clim$k * Cp_tirage * exp(eps_tirage)
Lc_eur     <- Q_ref * Lc_marche

L_total <- Lf_eur + Lc_eur

cat("\n[reference] Perte financiere  : moyenne", round(mean(Lf_eur), 3), "| sd", round(sd(Lf_eur), 3), "\n")
cat("[reference] Perte climatique  : moyenne", round(mean(Lc_eur), 3), "| sd", round(sd(Lc_eur), 3), "\n")
cat("[reference] Perte totale      : moyenne", round(mean(L_total), 3), "| sd", round(sd(L_total), 3), "\n")

cat("\nVaR (Md EUR) :\n")
print(round(rbind(financier = quantile(Lf_eur, qs), climatique = quantile(Lc_eur, qs),
                  total = quantile(L_total, qs)), 3))
cat("\nES (Md EUR) :\n")
print(round(c(financier = es(Lf_eur, .995), climatique = es(Lc_eur, .995), total = es(L_total, .995)), 3))

seuil <- quantile(L_total, .995)
idx_queue <- L_total >= seuil
cat("\nContribution moyenne dans la queue (au-dela du seuil ES99.5%) :\n")
print(round(c(financier = mean(Lf_eur[idx_queue]), climatique = mean(Lc_eur[idx_queue]),
              somme = mean(Lf_eur[idx_queue]) + mean(Lc_eur[idx_queue])), 3))

# ============================================================
# 4. GRILLE COMPLETE : Q x taille x (regime avant/apres) x delta
#    Calcule part_clim_ES995 ET part_clim_VaR995 dans la MEME boucle.
#    n_grille tirages par configuration (la queue a 0.5 % ne contient que
#    n_grille * 0.005 tirages : 10 000 tirages donnent une part bruitee).
# ============================================================
grille_Q      <- c(0.02, 0.05, 0.10)
grille_taille <- c(5, 10, 20, 50)
grille_scenario <- names(clim$Cp)   # d0_avant, d0.25_avant, d0.5_avant, d0_apres, d0.25_apres, d0.5_apres
n_grille <- 200000

res_grille <- expand.grid(Q = grille_Q, taille = grille_taille, scenario = grille_scenario,
                          stringsAsFactors = FALSE)
res_grille$regime <- ifelse(grepl("apres", res_grille$scenario), "2017+", "Before 2017")
res_grille$delta  <- as.numeric(gsub("d|_avant|_apres", "", res_grille$scenario))
res_grille$part_clim_ES995  <- NA
res_grille$part_clim_VaR995 <- NA

for (i in seq_len(nrow(res_grille))) {
  lam <- res_grille$taille[i] / 800
  set.seed(7000 + i)
  Lf_i  <- lam * sample(fin$L_net, n_grille, replace = TRUE)
  Cp_i  <- sample(clim$Cp[[res_grille$scenario[i]]], n_grille, replace = TRUE)
  eps_i <- sample(clim$eps_obs, n_grille, replace = TRUE)
  Lc_i  <- res_grille$Q[i] * clim$k * Cp_i * exp(eps_i)
  Lt_i  <- Lf_i + Lc_i
  
  idxq <- Lt_i >= quantile(Lt_i, .995)
  res_grille$part_clim_ES995[i]  <- mean(Lc_i[idxq]) / mean(Lt_i[idxq])
  res_grille$part_clim_VaR995[i] <- part_clim_at_VaR(Lc_i, Lt_i, p = .995)
  
  
  VaR_tot <- quantile(Lt_i, .995);  VaR_fin <- quantile(Lf_i, .995)
  ES_tot  <- mean(Lt_i[Lt_i >= VaR_tot]);  ES_fin <- mean(Lf_i[Lf_i >= VaR_fin])
  res_grille$dVaR[i] <- VaR_tot - VaR_fin
  res_grille$dES[i]  <- ES_tot  - ES_fin
}

print(res_grille %>%
        mutate(across(c(part_clim_ES995, part_clim_VaR995), ~ round(.x, 3))) %>%
        arrange(desc(part_clim_ES995)))

# ============================================================
# 5. ENREGISTREMENT (grille principale)
# ============================================================
saveRDS(list(res_grille = res_grille, L_total_ref = L_total, Lf_eur_ref = Lf_eur, Lc_eur_ref = Lc_eur,
             Q_ref = Q_ref, taille_ref = taille_ref, lambda_ref = lambda_ref, n_grille = n_grille,
             delta_covea_ancrage = delta_covea, hausse_pct_delta = data.frame(delta = grille_delta_test, hausse_pct = hausse_pct),
             note = paste("Grille Q x taille x regime x delta ; contribution du climat a l'ES99.5% et a la VaR99.5% du total.",
                          "Reference : Q=0.05, taille=20 MdEUR, scenario d0_apres.",
                          "Delta ancre sur Covea (2022), RCP 8.5, horizon 2050 : +70% de frequence d'eligibilite secheresse.",
                          "Grille simulee avec n_grille = 200000 tirages par configuration.")),
        file.path(dir_out, "assemblage_perte_totale.rds"))

# ============================================================
# 6. FIGURES
# ============================================================
## Figure 1 -- Hazard index
H <- readRDS(file.path(dir_out, "H_indice_alea_national.rds"))
p1 <- ggplot(H, aes(annee, H_pond)) +
  geom_col(aes(fill = H_pond > quantile(H_pond, .75, na.rm = TRUE)), width = 0.7, show.legend = FALSE) +
  scale_fill_manual(values = c("grey70", "#B2182B")) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  labs(x = NULL, y = "Hazard index (SWI anomaly, Jul-Dec)",
       title = "Drought hazard index, France, 2000-2025") +
  theme_journal
print(p1)
ggsave(file.path(dir_out, "fig1_hazard_index.png"), p1, width = 7, height = 3.8, dpi = 300)

## Figure 2 -- Recognition rate by hazard intensity and period
dd <- readRDS(file.path(dir_out, "table_annuelle_reconnaissance_modele.rds"))
p2 <- ggplot(dd, aes(jul_dec, taux, color = apres17)) +
  geom_point(size = 2.2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.6) +
  scale_color_manual(values = c("grey50", "#2166AC"), labels = c("Before 2017", "2017+")) +
  labs(x = "Hazard index (SWI anomaly, Jul-Dec)", y = "Recognition rate",
       title = "CatNat recognition rate by hazard intensity and period") +
  theme_journal
print(p2)
ggsave(file.path(dir_out, "fig2_recognition_regime.png"), p2, width = 6.5, height = 4.5, dpi = 300)

## Figure 3 -- Simulated annual loss distributions (reference configuration)
df_plot <- rbind(data.frame(loss_val = Lf_eur, type = "Financial"),
                 data.frame(loss_val = Lc_eur, type = "Climate"),
                 data.frame(loss_val = L_total, type = "Total"))
df_plot$type <- factor(df_plot$type, levels = c("Financial", "Climate", "Total"))

vlines <- data.frame(type = factor(c("Financial", "Climate", "Total"), levels(df_plot$type)),
                     var995 = c(quantile(Lf_eur, .995), quantile(Lc_eur, .995), quantile(L_total, .995)))

p3 <- ggplot(df_plot, aes(loss_val, fill = type, color = type)) +
  geom_density(alpha = 0.25, linewidth = 0.7) +
  geom_vline(data = vlines, aes(xintercept = var995, color = type), linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~type, scales = "free", ncol = 1) +
  labs(x = "Annual loss (EUR bn)", y = "Density", title = "Simulated annual loss distributions") +
  theme_journal + theme(legend.position = "none")
print(p3)
ggsave(file.path(dir_out, "fig3_loss_distributions.png"), p3, width = 6.5, height = 7.5, dpi = 300)

## Figure 4a -- Climate contribution, ES 99.5%
res_grille_fig <- res_grille %>%
  mutate(regime = factor(regime, levels = c("Before 2017", "2017+")))

p4_es <- ggplot(res_grille_fig, aes(x = factor(taille), y = factor(Q), fill = part_clim_ES995)) +
  geom_tile(color = "black", linewidth = 0.25) +
  geom_text(aes(label = scales::percent(part_clim_ES995, accuracy = 1),
                color = ifelse(part_clim_ES995 > 0.55, "white", "black")),
            size = 2.4, family = "serif") +
  scale_color_identity() +
  scale_fill_gradient(low = "#FFFFFF", high = "#212121", limits = c(0, 1),
                      labels = scales::percent_format(accuracy = 1), name = "Climate\nshare") +
  facet_grid(regime ~ delta, labeller = labeller(delta = function(x) paste0("delta = ", x))) +
  labs(x = "Insurer balance-sheet scale (EUR bn)", y = "Market share (Q)",
       title = "Climate contribution to tail risk (ES 99.5%)") +
  theme_journal
print(p4_es)
ggsave(file.path(dir_out, "fig4_full_grid_ES995.png"), p4_es, width = 9, height = 5, dpi = 300)

## Figure 4b -- Climate contribution, VaR 99.5%
p4_var <- ggplot(res_grille_fig, aes(x = factor(taille), y = factor(Q), fill = part_clim_VaR995)) +
  geom_tile(color = "black", linewidth = 0.25) +
  geom_text(aes(label = scales::percent(part_clim_VaR995, accuracy = 1),
                color = ifelse(part_clim_VaR995 > 0.55, "white", "black")),
            size = 2.4, family = "serif") +
  scale_color_identity() +
  scale_fill_gradient(low = "#FFFFFF", high = "#212121", limits = c(0, 1),
                      labels = scales::percent_format(accuracy = 1), name = "Climate\nshare") +
  facet_grid(regime ~ delta, labeller = labeller(delta = function(x) paste0("delta = ", x))) +
  labs(x = "Insurer balance-sheet scale (EUR bn)", y = "Market share (Q)",
       title = "Climate contribution at VaR 99.5% threshold") +
  theme_journal
print(p4_var)
ggsave(file.path(dir_out, "fig4_full_grid_VaR995.png"), p4_var, width = 9, height = 5, dpi = 300)

cat("\nFigures terminees. Fichiers de sortie dans", dir_out, "\n")

# ============================================================
# 7. ROBUSTESSE : independance financier/climat, bootstrap
# ============================================================
esg <- readRDS(file.path(dir_macro, "ESG_data_2000_2024_CAC40corrige.rds"))

an_cac <- as.integer(substr(esg$time, 1, 4))
cac_annuel <- tapply(esg$CAC40, an_cac, sum)
mrg <- merge(data.frame(annee = as.integer(names(cac_annuel)), cac = as.numeric(cac_annuel)),
             H[, c("annee", "H_pond")], by = "annee")
cat("\n[independance] Correlation de rang H_pond / rendement CAC40 annuel :",
    round(cor(mrg$H_pond, mrg$cac, method = "spearman"), 2), "| n =", nrow(mrg), "\n")

infl_annuel <- tapply(esg$Inflation, an_cac, sum)
mrg2 <- merge(mrg, data.frame(annee = as.integer(names(infl_annuel)), infl = as.numeric(infl_annuel)), by = "annee")
cat("[independance] Correlation de rang H_pond / inflation annuelle :",
    round(cor(mrg2$H_pond, mrg2$infl, method = "spearman"), 2), "\n")

set.seed(9000)
n_boot <- 300
n_sim_boot <- 3000

boot_q995 <- boot_es995 <- boot_part_clim <- numeric(n_boot)

for (b in seq_len(n_boot)) {
  annees_b  <- sample(clim$annees_alea, length(clim$annees_alea), replace = TRUE)
  Cp_pool_b <- clim$Cp$d0_apres[match(annees_b, clim$annees_alea)]
  
  Lf_b  <- lambda_ref * sample(fin$L_net, n_sim_boot, replace = TRUE)
  Cp_b  <- sample(Cp_pool_b, n_sim_boot, replace = TRUE)
  eps_b <- sample(clim$eps_obs, n_sim_boot, replace = TRUE)
  Lc_b  <- Q_ref * clim$k * Cp_b * exp(eps_b)
  Lt_b  <- Lf_b + Lc_b
  
  boot_q995[b]  <- quantile(Lt_b, .995)
  boot_es995[b] <- es(Lt_b, .995)
  idxq_b <- Lt_b >= quantile(Lt_b, .995)
  boot_part_clim[b] <- mean(Lc_b[idxq_b]) / mean(Lt_b[idxq_b])
}

cat("\nVaR99.5% : point central", round(quantile(L_total, .995), 3),
    "| IC bootstrap 95% :", round(quantile(boot_q995, c(.025, .975)), 3), "\n")
cat("ES99.5%  : point central", round(es(L_total, .995), 3),
    "| IC bootstrap 95% :", round(quantile(boot_es995, c(.025, .975)), 3), "\n")
cat("Part climat dans ES99.5% : point central", round(mean(Lc_eur[idx_queue]) / mean(L_total[idx_queue]), 3),
    "| IC bootstrap 95% :", round(quantile(boot_part_clim, c(.025, .975)), 3), "\n")

saveRDS(list(Q_ref = Q_ref, taille_ref = taille_ref,
             VaR995_central = as.numeric(quantile(L_total, .995)), VaR995_IC = as.numeric(quantile(boot_q995, c(.025, .975))),
             ES995_central = es(L_total, .995), ES995_IC = as.numeric(quantile(boot_es995, c(.025, .975))),
             part_clim_ES995_central = mean(Lc_eur[idx_queue]) / mean(L_total[idx_queue]),
             part_clim_ES995_IC = as.numeric(quantile(boot_part_clim, c(.025, .975))),
             corr_H_cac40 = round(cor(mrg$H_pond, mrg$cac, method = "spearman"), 2),
             corr_H_inflation = round(cor(mrg2$H_pond, mrg2$infl, method = "spearman"), 2),
             n_boot = n_boot, n_sim_boot = n_sim_boot,
             note = "Bootstrap sur les 22 annees d'alea observees (avec remise), 300 repliques. Test d'independance financier/climat : correlation de rang sur 25 ans (CAC40, inflation)."),
        file.path(dir_out, "robustesse_assemblage.rds"))

# ============================================================
# 8. SENSIBILITE A LA COMPOSITION DU PORTEFEUILLE STYLISE
#    Compare la part climatique a 3 compositions d'actifs differentes,
#    a taille totale egale (800 unites), passif de reference fixe
#    (Mp = 500, Dp = 6, Ap = 100), pour Q et scenario fixes.
#    Meme perte nette que la grille principale (perte_net2 dans 06_macro_esg.R).
#    Memes tirages pour les trois compositions (graine par taille), n_comp grand.
# ============================================================
a3r <- fin$a
stopifnot(all(c("cac", "infl", "d_court", "d_long") %in% colnames(a3r)))

perte_compo <- function(a, E, Ms, Ds, Ml, Dl, Mp, Dp, Api, Ap) {
  -(E * (exp(a[, "cac"] / 100) - 1) - (Ds * Ms) * (a[, "d_court"] / 100) -
      (Dl * Ml - Mp * Dp) * (a[, "d_long"] / 100) - (Api - Ap) * (a[, "infl"] / 100))
}

# reference     : composition de la grille principale (E=100, Ms=300, Ml=400)
# conservatrice : moins d'actions, plus d'obligations courtes
# actions       : plus d'actions
compositions <- list(
  reference     = list(E = 100, Ms = 300, Ds = 2, Ml = 400, Dl = 8, Mp = 500, Dp = 6, Api = 200, Ap = 100),
  conservatrice = list(E = 50,  Ms = 400, Ds = 2, Ml = 350, Dl = 8, Mp = 500, Dp = 6, Api = 200, Ap = 100),
  actions       = list(E = 150, Ms = 250, Ds = 2, Ml = 400, Dl = 8, Mp = 500, Dp = 6, Api = 200, Ap = 100)
)

# Verification : la composition de reference doit reproduire la perte nette fin$L_net
cp0 <- compositions$reference
stopifnot(isTRUE(all.equal(
  perte_compo(a3r, cp0$E, cp0$Ms, cp0$Ds, cp0$Ml, cp0$Dl, cp0$Mp, cp0$Dp, cp0$Api, cp0$Ap),
  as.numeric(fin$L_net))))

Q_test <- 0.05
tailles_test <- c(5, 50)   # les deux extremes de la grille
scenario_test <- "d0_apres"
n_comp <- 200000

res_compo <- expand.grid(composition = names(compositions), taille = tailles_test,
                         stringsAsFactors = FALSE)
res_compo$part_clim_ES995 <- NA

for (i in seq_len(nrow(res_compo))) {
  cp  <- compositions[[res_compo$composition[i]]]
  lam <- res_compo$taille[i] / 800
  set.seed(9000 + res_compo$taille[i])     # meme graine pour les 3 compositions
  
  Lf_i <- lam * perte_compo(a3r[sample(nrow(a3r), n_comp, replace = TRUE), ],
                            cp$E, cp$Ms, cp$Ds, cp$Ml, cp$Dl, cp$Mp, cp$Dp, cp$Api, cp$Ap)
  Cp_i  <- sample(clim$Cp[[scenario_test]], n_comp, replace = TRUE)
  eps_i <- sample(clim$eps_obs, n_comp, replace = TRUE)
  Lc_i  <- Q_test * clim$k * Cp_i * exp(eps_i)
  
  Lt_i <- Lf_i + Lc_i
  idxq <- Lt_i >= quantile(Lt_i, .995)
  res_compo$part_clim_ES995[i] <- mean(Lc_i[idxq]) / mean(Lt_i[idxq])
}

cat("\n[sensibilite composition] Part climatique selon la composition du portefeuille :\n")
print(res_compo %>% mutate(part_clim_ES995 = round(part_clim_ES995, 3)) %>%
        arrange(taille, composition))

cat("\n[controle] Grille principale, meme cellule (Q=0.05, d0_apres) :\n")
print(res_grille %>% filter(Q == Q_test, scenario == scenario_test, taille %in% tailles_test) %>%
        select(taille, part_clim_ES995) %>% mutate(part_clim_ES995 = round(part_clim_ES995, 3)))

saveRDS(res_compo, file.path(dir_out, "sensibilite_composition_portefeuille.rds"))

# Ecart-type, quantile 99.5 % et sensibilite nette de la perte financiere
# (800 unites) selon la composition : colonnes du tableau d'annexe.
sd_compo <- sapply(compositions, function(cp) {
  L <- perte_compo(a3r, cp$E, cp$Ms, cp$Ds, cp$Ml, cp$Dl, cp$Mp, cp$Dp, cp$Api, cp$Ap)
  c(ecart_type = sd(L), q995 = unname(quantile(L, .995)),
    sens_taux_long_nette = cp$Dl * cp$Ml - cp$Mp * cp$Dp)
})
cat("\n[sensibilite composition] Perte financiere (800 unites) selon la composition :\n")
print(round(t(sd_compo), 1))

# ============================================================
# 9. SENSIBILITE A LA LOI DU RESIDU (annexe)
#    Compare trois lois pour eps : residus empiriques (reference),
#    normale N(0, sigma_eps^2), normale ecretee a +/- 1.5 sigma_eps.
#    Perte climatique de MARCHE (avant part de marche Q), Md EUR.
# ============================================================
sig_eps <- sd(clim$eps_obs)
borne   <- 1.5 * sig_eps

tirer_eps <- function(n, bruit = c("empirique", "normal", "ecrete")) {
  bruit <- match.arg(bruit)
  switch(bruit,
         empirique = sample(clim$eps_obs, n, replace = TRUE),
         normal    = rnorm(n, 0, sig_eps),
         ecrete    = pmax(pmin(rnorm(n, 0, sig_eps), borne), -borne))
}

n_sens <- 1e6    # grand n : la queue de la normale est instable avec 1e4 tirages
modes  <- c("empirique", "normal", "ecrete")
scen   <- c(avant = "d0_avant", apres = "d0_apres")

res_bruit <- expand.grid(regime = names(scen), bruit = modes, stringsAsFactors = FALSE)
res_bruit$moyenne <- res_bruit$q995 <- res_bruit$q999 <- res_bruit$maximum <- NA
res_bruit$part_clim_ES995 <- NA   # configuration de reference (Q_ref, taille_ref)

for (i in seq_len(nrow(res_bruit))) {
  set.seed(11000 + i)
  Cp_i <- sample(clim$Cp[[scen[res_bruit$regime[i]]]], n_sens, replace = TRUE)
  Lc_m <- clim$k * Cp_i * exp(tirer_eps(n_sens, res_bruit$bruit[i]))
  res_bruit$moyenne[i] <- mean(Lc_m)
  res_bruit$q995[i]    <- quantile(Lc_m, .995)
  res_bruit$q999[i]    <- quantile(Lc_m, .999)
  res_bruit$maximum[i] <- max(Lc_m)
  
  Lf_i <- lambda_ref * sample(fin$L_net, n_sens, replace = TRUE)
  Lt_i <- Lf_i + Q_ref * Lc_m
  idxq <- Lt_i >= quantile(Lt_i, .995)
  res_bruit$part_clim_ES995[i] <- mean(Q_ref * Lc_m[idxq]) / mean(Lt_i[idxq])
}

cat("\n[sensibilite loi du residu] Perte climatique de marche (Md EUR), delta = 0 :\n")
print(res_bruit %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
        arrange(regime, match(bruit, modes)))

saveRDS(res_bruit, file.path(dir_out, "sensibilite_loi_residu.rds"))
write.csv(res_bruit, file.path(dir_out, "sensibilite_loi_residu.csv"), row.names = FALSE)

cat("\nScript 07_assemblage.R termine. Fichiers de sortie dans", dir_out, "\n")


######################################## Accesoires 
library(ggplot2)

dd <- readRDS(file.path(dir_out, "table_annuelle_reconnaissance_modele.rds"))
dd$haut <- dd$jul_dec > quantile(dd$jul_dec, .75)   # dernier quartile

# contrôle : doit donner 0.84 et 1.03 (texte de la 4.2)
cat(round(mean(dd$jul_dec[dd$annee %in% 2003:2016]), 2),
    round(mean(dd$jul_dec[dd$annee %in% 2017:2023]), 2), "\n")

p_haz <- ggplot(dd, aes(annee, jul_dec, fill = haut)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey30") +
  geom_text(data = subset(dd, annee == 2022), aes(label = annee),
            vjust = -0.5, size = 3, fontface = "bold", colour = "black") +
  scale_fill_manual(values = c("grey65", "#2166AC")) +
  scale_x_continuous(breaks = c(2005, 2010, 2015, 2020)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(x = "Year", y = "National hazard index") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92"),
    panel.border       = element_rect(colour = "grey70"),
    axis.title         = element_text(face = "bold"),
    axis.text          = element_text(colour = "black")
  )

print(p_haz)
ggsave(file.path(dir_out, "fig_hazard_index.png"), p_haz,
       width = 6.3, height = 2.8, dpi = 300)

