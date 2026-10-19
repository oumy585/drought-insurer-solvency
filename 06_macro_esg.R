
# ============================================================
# 06_macro_esg.R : ESG macro-financier (Module 1), sans climat
# Point de depart : fichier corrige (CAC40 remis dans le bon sens)
# ============================================================
library(dplyr); library(ahead)
dir_macro <- "C:/Users/OumyDIONE/Documents/output/macro"
esg <- readRDS(file.path(dir_macro, "ESG_data_2000_2024_CAC40corrige.rds"))
vars <- c("CAC40", "IPI_growth", "Inflation", "ShortRate", "LongRate")
stopifnot(!is.unsorted(esg$time), nrow(esg) == 300, !anyNA(esg[vars]))

Y <- ts(as.matrix(esg[vars]), start = c(2000, 1), frequency = 12)
cat("Y :", dim(Y), "|", paste(colnames(Y), collapse = ", "), "\n")

cat("\nArguments de ridge2f :\n"); print(names(formals(ahead::ridge2f)))

set.seed(1)
fit <- ahead::ridge2f(Y, h = 12, lags = 5, nb_hidden = 20, lambda_1 = 2.05, lambda_2 = 7.1,
                      nodes_sim = "sobol", activ = "relu",
                      type_pi = "movingblockbootstrap", block_length = 12, B = 5, seed = 1)
cat("\nElements de la sortie :", names(fit), "\n")
str(fit, max.level = 1, list.len = 12)
if (!is.null(fit$sims)) {
  cat("\nClasse de sims :", class(fit$sims), "| longueur :", length(fit$sims), "\n")
  s1 <- if (is.list(fit$sims)) fit$sims[[1]] else fit$sims
  cat("Premiere simulation : classe", class(s1), "| dim", dim(s1), "\n"); print(round(head(s1, 3), 2))
}
print(round(fit$mean, 2))

############################################################################### verification 
# 1) Les chemins simules sont-ils "prevision centrale + bloc de residus" ?
R  <- as.matrix(fit$residuals)
for (b in seq_along(fit$sims)) {
  d1 <- as.matrix(fit$sims[[b]]) - as.matrix(fit$mean)
  ecart <- sapply(1:(nrow(R) - 11), function(s) max(abs(R[s:(s + 11), ] - d1)))
  cat("Simulation", b, ": meilleur bloc de residus, debut", which.min(ecart), "| ecart maximal", signif(min(ecart), 3), "\n")
}

# 2) Moyenne sur 12 mois de la prevision centrale contre l'historique
v3 <- c("CAC40", "IPI_growth", "Inflation")
print(round(rbind(prevision_12m = colSums(as.matrix(fit$mean))[v3],
                  historique_12m = 12 * colMeans(as.matrix(Y))[v3]), 2))
print(round(rbind(variance = apply(as.matrix(Y), 2, var), loocv = fit$loocv_per_series,
                  rapport = fit$loocv_per_series / apply(as.matrix(Y), 2, var)), 3))

########################################################## Dispersion  des chemins  selon l'horizon 
set.seed(2)
fit2 <- ahead::ridge2f(Y, h = 12, lags = 5, nb_hidden = 20, lambda_1 = 2.05, lambda_2 = 7.1,
                       nodes_sim = "sobol", activ = "relu", type_pi = "movingblockbootstrap",
                       block_length = 12, B = 200, seed = 2, show_progress = FALSE)
S <- simplify2array(lapply(fit2$sims, as.matrix))            # 12 x 5 x 200
sd_h <- apply(S, c(1, 2), sd); colnames(sd_h) <- colnames(Y)
cat("Ecart-type des chemins simules, par horizon (mois 1, 3, 6, 9, 12) :\n"); print(round(sd_h[c(1, 3, 6, 9, 12), ], 2))

cat("\nEcart-type des residus (un mois) :\n"); print(round(apply(as.matrix(fit2$residuals), 2, sd), 2))
d12 <- function(v) sd(diff(as.numeric(v), lag = 12))
cat("\nEcart-type historique des variations sur 12 mois (niveaux de taux) :\n")
print(round(c(ShortRate = d12(Y[, "ShortRate"]), LongRate = d12(Y[, "LongRate"])), 2))
cat("\nMoyenne des chemins a h = 12 contre prevision centrale :\n")
print(round(rbind(sims = apply(S[12, , ], 1, mean), centrale = as.numeric(fit2$mean[12, ])), 2))

#########################################################################################

set.seed(3)
fit3 <- ahead::ridge2f(Y, h = 12, lags = 5, nb_hidden = 20, lambda_1 = 2.05, lambda_2 = 7.1,
                       nodes_sim = "sobol", activ = "relu", type_pi = "movingblockbootstrap",
                       block_length = 12, B = 1000, seed = 3, show_progress = FALSE)
S3 <- simplify2array(lapply(fit3$sims, as.matrix)); B3 <- dim(S3)[3]
x0 <- as.numeric(tail(Y, 1))                                    # dernier etat observe (2024-12)
cat("Dernier etat observe :", round(x0, 2), "\n")

nm <- c("cac_12m", "infl_12m", "d_court", "d_long")
sim_a <- t(sapply(1:B3, function(b) c(sum(S3[, 1, b]), sum(S3[, 3, b]), S3[12, 4, b] - x0[4], S3[12, 5, b] - x0[5])))
Ym <- as.matrix(Y)
his_a <- t(sapply(1:(nrow(Ym) - 12), function(t) c(sum(Ym[(t + 1):(t + 12), 1]), sum(Ym[(t + 1):(t + 12), 3]),
                                                   Ym[t + 12, 4] - Ym[t, 4], Ym[t + 12, 5] - Ym[t, 5])))
colnames(sim_a) <- colnames(his_a) <- nm

f <- function(M) round(rbind(moyenne = colMeans(M), ecart_type = apply(M, 2, sd),
                             q05 = apply(M, 2, quantile, .05), q95 = apply(M, 2, quantile, .95)), 2)
cat("\nSIMULE (12 mois a partir de 2024-12) :\n"); print(f(sim_a))
cat("\nHISTORIQUE (288 fenetres glissantes de 12 mois) :\n"); print(f(his_a))
cat("\nCorrelations simulees :\n"); print(round(cor(sim_a), 2))
cat("\nCorrelations historiques :\n"); print(round(cor(his_a), 2))

perte <- function(a) -(100 * (exp(a[, 1] / 100) - 1) - 600 * (a[, 3] / 100) - 3200 * (a[, 4] / 100) - 200 * (a[, 2] / 100))
Ls <- perte(sim_a); Lh <- perte(his_a)
cat("\nPerte annuelle (unites du portefeuille stylise) :\n")
print(round(rbind(simulee = c(moy = mean(Ls), sd = sd(Ls), q95 = quantile(Ls, .95), q99 = quantile(Ls, .99), q995 = quantile(Ls, .995)),
                  historique = c(mean(Lh), sd(Lh), quantile(Lh, .95), quantile(Lh, .99), quantile(Lh, .995))), 1))


################################################### Robustesse du similateur 
cfg <- data.frame(nom = c("actuel", "lam1_10^2.05", "peu_regularise", "lags1", "lags3", "bloc6", "bloc24"),
                  lags = c(5, 5, 5, 1, 3, 5, 5), l1 = c(2.05, 10^2.05, 0.1, 2.05, 2.05, 2.05, 2.05),
                  l2 = c(7.1, 7.1, 0.1, 7.1, 7.1, 7.1, 7.1), bloc = c(12, 12, 12, 12, 12, 6, 24))
run_cfg <- function(i) {
  fx <- ahead::ridge2f(Y, h = 12, lags = cfg$lags[i], nb_hidden = 20, lambda_1 = cfg$l1[i], lambda_2 = cfg$l2[i],
                       nodes_sim = "sobol", activ = "relu", type_pi = "movingblockbootstrap",
                       block_length = cfg$bloc[i], B = 500, seed = 10 + i, show_progress = FALSE)
  Sx <- simplify2array(lapply(fx$sims, as.matrix)); Bx <- dim(Sx)[3]
  a <- t(sapply(1:Bx, function(b) c(sum(Sx[, 1, b]), sum(Sx[, 3, b]), Sx[12, 4, b] - x0[4], Sx[12, 5, b] - x0[5])))
  L <- perte(a)
  c(sd_cac = sd(a[, 1]), sd_infl = sd(a[, 2]), mu_court = mean(a[, 3]), sd_court = sd(a[, 3]),
    mu_long = mean(a[, 4]), sd_long = sd(a[, 4]), perte_sd = sd(L), q95 = unname(quantile(L, .95)), q995 = unname(quantile(L, .995)))
}
res_cfg <- t(sapply(seq_len(nrow(cfg)), run_cfg)); rownames(res_cfg) <- cfg$nom
print(round(res_cfg, 2))

#################################################fixer le simulateur, puis mesurer l'erreur de simulation
set.seed(2026)
fit_ref <- ahead::ridge2f(Y, h = 12, lags = 5, nb_hidden = 20, lambda_1 = 2.05, lambda_2 = 7.1,
                          nodes_sim = "sobol", activ = "relu", type_pi = "movingblockbootstrap",
                          block_length = 12, B = 10000, seed = 2026, show_progress = FALSE)
S_ref <- simplify2array(lapply(fit_ref$sims, as.matrix)); B_ref <- dim(S_ref)[3]
a_ref <- cbind(cac = apply(S_ref[, 1, ], 2, sum), infl = apply(S_ref[, 3, ], 2, sum),
               d_court = S_ref[12, 4, ] - x0[4], d_long = S_ref[12, 5, ] - x0[5])
L_ref <- perte(a_ref)
qs <- c(.95, .99, .995)
cat("B =", B_ref, "| perte : moyenne", round(mean(L_ref), 1), "| sd", round(sd(L_ref), 1), "\n")
print(round(quantile(L_ref, qs), 1))
es <- function(L, p) mean(L[L >= quantile(L, p)])
cat("ES 95 / 99 / 99.5 :", round(c(es(L_ref, .95), es(L_ref, .99), es(L_ref, .995)), 1), "\n")

set.seed(7)
bq <- t(replicate(500, { l <- sample(L_ref, replace = TRUE); c(quantile(l, qs), es(l, .995)) }))
colnames(bq) <- c("q95", "q99", "q995", "ES995")
cat("Intervalle a 95 % sur les quantiles (bootstrap des chemins) :\n")
print(round(apply(bq, 2, quantile, c(.025, .5, .975)), 1))
saveRDS(list(a = a_ref, L = L_ref, x0 = x0, perte_formule = "L = -(100*(exp(cac/100)-1) - 600*d_court/100 - 3200*d_long/100 - 200*infl/100)",
             config = "ridge2f lags=5 nb_hidden=20 lambda_1=2.05 lambda_2=7.1 sobol relu, movingblockbootstrap block=12, B=10000, seed=2026"),
        file.path(dir_macro, "simulation_financiere_reference.rds"))

########################## Test 
cat("Pertes distinctes (arrondi a 3 decimales) :", length(unique(round(L_ref, 3))), "sur", length(L_ref), "\n")
print(round(sort(L_ref, decreasing = TRUE)[1:15], 2))
cat("Repetitions de la pire valeur :", sum(round(L_ref, 3) == round(max(L_ref), 3)), "\n")


set.seed(2027)
fit_ref2 <- ahead::ridge2f(Y, h = 12, lags = 5, nb_hidden = 20, lambda_1 = 2.05, lambda_2 = 7.1,
                           nodes_sim = "sobol", activ = "relu", type_pi = "movingblockbootstrap",
                           block_length = 3, B = 10000, seed = 2027, show_progress = FALSE)
S2 <- simplify2array(lapply(fit_ref2$sims, as.matrix))
a2 <- cbind(cac = apply(S2[, 1, ], 2, sum), infl = apply(S2[, 3, ], 2, sum),
            d_court = S2[12, 4, ] - x0[4], d_long = S2[12, 5, ] - x0[5])
L2 <- perte(a2)

cat("Pertes distinctes :", length(unique(round(L2, 3))), "sur", length(L2), "\n")
cat("Perte : moyenne", round(mean(L2), 1), "| sd", round(sd(L2), 1), "\n")
print(round(quantile(L2, qs), 1))
cat("ES 95/99/99.5 :", round(c(es(L2, .95), es(L2, .99), es(L2, .995)), 1), "\n")

set.seed(8)
bq2 <- t(replicate(500, { l <- sample(L2, replace = TRUE); c(quantile(l, qs), es(l, .995)) }))
colnames(bq2) <- c("q95", "q99", "q995", "ES995")
print(round(apply(bq2, 2, quantile, c(.025, .5, .975)), 1))

#############################################
for (bl in c(3, 6, 9, 12)) {
  set.seed(3000 + bl)
  fx <- ahead::ridge2f(Y, h = 12, lags = 5, nb_hidden = 20, lambda_1 = 2.05, lambda_2 = 7.1,
                       nodes_sim = "sobol", activ = "relu", type_pi = "movingblockbootstrap",
                       block_length = bl, B = 5000, seed = 3000 + bl, show_progress = FALSE)
  Sx <- simplify2array(lapply(fx$sims, as.matrix))
  ax <- cbind(cac = apply(Sx[, 1, ], 2, sum), infl = apply(Sx[, 3, ], 2, sum),
              d_court = Sx[12, 4, ] - x0[4], d_long = Sx[12, 5, ] - x0[5])
  Lx <- perte(ax)
  cat("bloc =", bl, "| distinctes :", length(unique(round(Lx, 3))), "/ 5000",
      "| VaR99.5 :", round(quantile(Lx, .995), 1), "| ES99.5 :", round(es(Lx, .995), 1), "\n")
}

saveRDS(list(a = a2, L = L2, x0 = x0, perte_formule = "L = -(100*(exp(cac/100)-1) - 600*d_court/100 - 3200*d_long/100 - 200*infl/100)",
             config = "ridge2f lags=5 nb_hidden=20 lambda_1=2.05 lambda_2=7.1 sobol relu, movingblockbootstrap block=3, B=10000, seed=2027",
             note = "Reference retenue : block_length=3 (block_length=12 rejette car seulement 284/10000 pertes distinctes)"),
        file.path(dir_macro, "simulation_financiere_reference.rds"))
list.files(dir_macro)

##########################################
# Verifier l'hypothese 2 : la variance de la perte suit-elle le nombre de blocs par chemin ?
res_bl <- data.frame(bloc = c(3, 6, 9, 12), nb_blocs_approx = 12 / c(3, 6, 9, 12),
                     sd_perte = NA, var_ratio_attendu = NA)
for (i in seq_along(res_bl$bloc)) {
  bl <- res_bl$bloc[i]
  set.seed(4000 + bl)
  fx <- ahead::ridge2f(Y, h = 12, lags = 5, nb_hidden = 20, lambda_1 = 2.05, lambda_2 = 7.1,
                       nodes_sim = "sobol", activ = "relu", type_pi = "movingblockbootstrap",
                       block_length = bl, B = 5000, seed = 4000 + bl, show_progress = FALSE)
  Sx <- simplify2array(lapply(fx$sims, as.matrix))
  ax <- cbind(cac = apply(Sx[, 1, ], 2, sum), infl = apply(Sx[, 3, ], 2, sum),
              d_court = Sx[12, 4, ] - x0[4], d_long = Sx[12, 5, ] - x0[5])
  res_bl$sd_perte[i] <- sd(perte(ax))
}
print(res_bl)
cat("\nSi effet purement mecanique (moins de blocs = moins de lissage), sd_perte devrait croitre",
    "a peu pres comme sqrt(nb_blocs_max / nb_blocs) : verifier la forme\n")


################################################################################ reference 
set.seed(2028)
fit_ref3 <- ahead::ridge2f(Y, h = 12, lags = 5, nb_hidden = 20, lambda_1 = 2.05, lambda_2 = 7.1,
                           nodes_sim = "sobol", activ = "relu", type_pi = "movingblockbootstrap",
                           block_length = 6, B = 10000, seed = 2028, show_progress = FALSE)
S3r <- simplify2array(lapply(fit_ref3$sims, as.matrix))
a3r <- cbind(cac = apply(S3r[, 1, ], 2, sum), infl = apply(S3r[, 3, ], 2, sum),
             d_court = S3r[12, 4, ] - x0[4], d_long = S3r[12, 5, ] - x0[5])
L3r <- perte(a3r)

cat("Pertes distinctes :", length(unique(round(L3r, 3))), "sur", length(L3r), "\n")
cat("Perte : moyenne", round(mean(L3r), 1), "| sd", round(sd(L3r), 1), "\n")
print(round(quantile(L3r, qs), 1))
cat("ES 95/99/99.5 :", round(c(es(L3r, .95), es(L3r, .99), es(L3r, .995)), 1), "\n")

set.seed(9)
bq3 <- t(replicate(500, { l <- sample(L3r, replace = TRUE); c(quantile(l, qs), es(l, .995)) }))
colnames(bq3) <- c("q95", "q99", "q995", "ES995")
print(round(apply(bq3, 2, quantile, c(.025, .5, .975)), 1))

saveRDS(list(a = a3r, L = L3r, x0 = x0,
             perte_formule = "L = -(100*(exp(cac/100)-1) - 600*d_court/100 - 3200*d_long/100 - 200*infl/100)",
             config = "ridge2f lags=5 nb_hidden=20 lambda_1=2.05 lambda_2=7.1 sobol relu, movingblockbootstrap block=6, B=10000, seed=2028",
             note = "Reference retenue : block_length=6 (compromis diversite de scenarios / dependance temporelle ; block_length=12 rejete, 284/10000 pertes distinctes ; block_length=3 et 9 en sensibilite)"),
        file.path(dir_macro, "simulation_financiere_reference.rds"))
list.files(dir_macro)


######################################################################### Passif stylisé 
perte_net <- function(a, Mp, Dp, Ap) {
  -(100 * (exp(a[, 1] / 100) - 1) - 600 * (a[, 3] / 100) - 3200 * (a[, 4] / 100) - 200 * (a[, 2] / 100)) +
    (-Mp * Dp * (a[, 4] / 100) - Ap * (a[, 2] / 100))
}

grille <- expand.grid(Mp = c(0, 300, 500, 700), Dp = c(4, 6, 8), Ap = c(0, 100, 200))
res_p <- t(apply(grille, 1, function(r) {
  L <- perte_net(a3r, r["Mp"], r["Dp"], r["Ap"])
  c(sens_taux_long_nette = 3200 - r["Mp"] * r["Dp"], sd = sd(L), q995 = quantile(L, .995), es995 = es(L, .995))
}))
print(round(cbind(grille, res_p), 1))


####################################### enregistrement de la fonction de perte 
Mp_ref <- 500; Dp_ref <- 6; Ap_ref <- 100
L_net <- perte_net(a3r, Mp_ref, Dp_ref, Ap_ref)

cat("Perte nette : moyenne", round(mean(L_net), 1), "| sd", round(sd(L_net), 1), "\n")
print(round(quantile(L_net, qs), 1))
cat("ES 95/99/99.5 :", round(c(es(L_net, .95), es(L_net, .99), es(L_net, .995)), 1), "\n")

saveRDS(list(a = a3r, L_brut = L3r, L_net = L_net, x0 = x0,
             portefeuille = list(E = 100, Ms = 300, Ml = 400, Ds = 2, Dl = 8, Api = 200,
                                 Mp = Mp_ref, Dp = Dp_ref, Ap = Ap_ref),
             perte_formule_nette = "L = -(100*(exp(cac/100)-1) - 600*d_court/100 - 3200*d_long/100 - 200*infl/100) - Mp*Dp*d_long/100 - Ap*infl/100",
             config_sim = "ridge2f lags=5 nb_hidden=20 lambda_1=2.05 lambda_2=7.1 sobol relu, movingblockbootstrap block=6, B=10000, seed=2028",
             grille_sensibilite_passif = grille, res_sensibilite_passif = res_p,
             note = "Passif stylise retenu : Mp=500, Dp=6, Ap=100 (sensibilite nette taux long = 200, positive et moderee)"),
        file.path(dir_macro, "simulation_financiere_reference.rds"))
list.files(dir_macro)

cat("Moyenne de la variation de taux long simulee sur 12 mois :", round(mean(a3r[, "d_long"]), 3), "\n")
cat("Dernier taux long observe :", round(x0["LongRate"], 2), "\n")
cat("Contribution moyenne du terme de couverture (-Mp*Dp*d_long/100) :",
    round(mean(-Mp_ref * Dp_ref * a3r[, "d_long"] / 100), 2), "\n")
cat("Contribution moyenne du terme actif seul (-3200*d_long/100) :",
    round(mean(-3200 * a3r[, "d_long"] / 100), 2), "\n")

############################################## correction perte 
perte_net2 <- function(a, Mp, Dp, Ap) {
  -(100 * (exp(a[, 1] / 100) - 1) - 600 * (a[, 3] / 100) - (3200 - Mp * Dp) * (a[, 4] / 100) - (200 - Ap) * (a[, 2] / 100))
}
L_net2 <- perte_net2(a3r, Mp_ref, Dp_ref, Ap_ref)
cat("Sensibilite nette taux long :", 3200 - Mp_ref * Dp_ref, "| inflation :", 200 - Ap_ref, "\n")
cat("Perte nette corrigee : moyenne", round(mean(L_net2), 1), "| sd", round(sd(L_net2), 1), "\n")
print(round(quantile(L_net2, qs), 1))
cat("ES 95/99/99.5 :", round(c(es(L_net2, .95), es(L_net2, .99), es(L_net2, .995)), 1), "\n")


contrib_eq   <- -100 * (exp(a3r[, "cac"] / 100) - 1)
contrib_ct   <- 600 * (a3r[, "d_court"] / 100)
contrib_lt   <- (3200 - Mp_ref * Dp_ref) * (a3r[, "d_long"] / 100)
contrib_infl <- (200 - Ap_ref) * (a3r[, "infl"] / 100)

cat("Contributions moyennes a la perte nette :\n")
print(round(c(equity = mean(contrib_eq), taux_court = mean(contrib_ct),
              taux_long_net = mean(contrib_lt), inflation_net = mean(contrib_infl),
              somme = mean(contrib_eq) + mean(contrib_ct) + mean(contrib_lt) + mean(contrib_infl)), 2))

cat("\nMoyennes des variables sous-jacentes :\n")
print(round(c(cac_12m = mean(a3r[, "cac"]), d_court = mean(a3r[, "d_court"]),
              d_long = mean(a3r[, "d_long"]), infl_12m = mean(a3r[, "infl"])), 3))

