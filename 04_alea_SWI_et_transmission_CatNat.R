# ============================================================
# 04_alea_SWI_et_transmission_CatNat.R
# Drought hazard index (communal SWI)
# Depends only on: SWI_communal_mensuel.rds
# ============================================================
library(dplyr); library(tidyr); library(stringr)

dir_out <- file.path("output", "etape0")
print(list.files(dir_out))

f_swi_rds <- file.path(dir_out, "SWI_communal_mensuel.rds")
if (!file.exists(f_swi_rds)) stop("SWI_communal_mensuel.rds not found: run 02_construction_alea_SWI.R first")

B <- readRDS(f_swi_rds)
cat("Communal SWI:", nrow(B), "communes x", ncol(B), "months |",
    colnames(B)[1], "->", colnames(B)[ncol(B)], "| NA:", sum(is.na(B)), "\n")
#################################################################### Monthly swi 
nat <- data.frame(mois = colnames(B),
                  moyenne  = colMeans(B),
                  p10      = apply(B, 2, quantile, probs = 0.10),
                  part_sec = colMeans(B < 0.10)) %>%
  mutate(annee = as.integer(substr(mois, 1, 4)), m = as.integer(substr(mois, 5, 6)))

cat("SWI moyen national (mois) :\n"); print(round(summary(nat$moyenne), 3))

# Cycle saisonnier : moyenne par mois calendaire, sur 2000-2025
print(round(tapply(nat$moyenne, nat$m, mean), 3))

# Grandes annees : SWI moyen de avril a decembre
grandes <- c(2003, 2011, 2018, 2019, 2020, 2022)
print(nat %>% filter(annee %in% grandes, m %in% 4:12) %>%
        select(annee, m, moyenne) %>%
        mutate(moyenne = round(moyenne, 3)) %>%
        pivot_wider(names_from = m, values_from = moyenne), n = 10)

# Minimum annuel du SWI moyen national et mois correspondant
print(nat %>% group_by(annee) %>%
        summarise(min_moy = round(min(moyenne), 3), mois_min = m[which.min(moyenne)],
                  max_part_sec = round(max(part_sec), 3)) %>% as.data.frame(), row.names = FALSE)



############################################
serie <- readRDS(file.path(dir_out, "serie_nationale_annuelle.rds"))
ser2  <- readRDS(file.path(dir_out, "serie_reconnaissance_hors_annuelle.rds"))

an <- as.integer(substr(colnames(B), 1, 4)); mm <- as.integer(substr(colnames(B), 5, 6))
fen <- list(T1 = 1:3, T2 = 4:6, T3 = 7:9, T4 = 10:12, avr_sep = 4:9, jul_dec = 7:12)

int_fen <- function(ms) sapply(2000:2025, function(y) {
  cols <- which(an == y & mm %in% ms)
  mean(apply(B[, cols, drop = FALSE], 1, min))         # moyenne nationale des minima communaux
})
H <- sapply(fen, int_fen); rownames(H) <- 2000:2025
print(round(H, 3))

d <- serie %>% filter(annee >= 2002, annee <= 2023) %>% select(annee, a_pond, rho_pond) %>%
  left_join(ser2 %>% select(annee, taux = taux_hors_annuelle), by = "annee") %>%
  left_join(data.frame(annee = 2000:2025, H), by = "annee")

av <- d$annee %in% 2003:2016; ap <- d$annee >= 2017
cors <- sapply(names(fen), function(f) c(
  a_pond       = cor(d$a_pond,   d[[f]], method = "spearman"),
  rho_pond     = cor(d$rho_pond, d[[f]], method = "spearman"),
  taux_avant17 = cor(d$taux[av], d[[f]][av], method = "spearman"),
  taux_apres17 = cor(d$taux[ap], d[[f]][ap], method = "spearman")))
print(round(cors, 2))

#################################################
an  <- as.integer(substr(colnames(B), 1, 4)); mm <- as.integer(substr(colnames(B), 5, 6))
ref <- an <= 2024

Z <- B
for (mo in 1:12) {
  cols <- which(mm == mo); cr <- which(mm == mo & ref)
  mu <- rowMeans(B[, cr, drop = FALSE])
  sdv <- apply(B[, cr, drop = FALSE], 1, sd)
  sdv[sdv < 1e-6] <- NA                                   # communes sans variation ce mois-la
  Z[, cols] <- (B[, cols, drop = FALSE] - mu) / sdv
}
Z <- -Z                                                    # sens : valeur elevee = plus sec
cat("Anomalies : NA =", sum(is.na(Z)), "| communes avec au moins un NA :", sum(rowSums(is.na(Z)) > 0), "\n")
print(round(summary(as.vector(Z[, an == 2022 & mm == 8])), 2))    # 2022, aout
print(round(summary(as.vector(Z[, an == 2013 & mm == 8])), 2))    # annee calme, aout
print(round(apply(Z[, an == 2003 & mm %in% 7:9], 2, mean), 2))     # 2003, moyenne nationale

#################################################################
fen <- list(T2 = 4:6, T3 = 7:9, T4 = 10:12, avr_sep = 4:9, jul_dec = 7:12)

int_z <- function(ms) sapply(2000:2025, function(y) {
  cols <- which(an == y & mm %in% ms)
  mean(apply(Z[, cols, drop = FALSE], 1, max))          # moyenne nationale du maximum communal de l'anomalie
})
HZ <- sapply(fen, int_z); rownames(HZ) <- 2000:2025
print(round(HZ, 2))

d2 <- d %>% select(annee, a_pond, rho_pond, taux) %>% left_join(data.frame(annee = 2000:2025, HZ), by = "annee")
av <- d2$annee %in% 2003:2016; ap <- d2$annee >= 2017
cz <- sapply(names(fen), function(f) c(
  a_pond = cor(d2$a_pond, d2[[f]], method = "spearman"),
  rho_pond = cor(d2$rho_pond, d2[[f]], method = "spearman"),
  taux_avant17 = cor(d2$taux[av], d2[[f]][av], method = "spearman"),
  taux_apres17 = cor(d2$taux[ap], d2[[f]][ap], method = "spearman")))
print(round(cz, 2))

#####################################################le saut de niveau tient-il à intensité égale ?
dd <- d2 %>% filter(annee >= 2003) %>%
  left_join(ser2 %>% select(annee, reconnues_trimestrielle, demandes_trimestrielle), by = "annee") %>%
  mutate(rec = reconnues_trimestrielle, non = demandes_trimestrielle - reconnues_trimestrielle,
         apres17 = annee >= 2017)

print(dd %>% group_by(apres17) %>% summarise(n = n(), intensite_moy = round(mean(jul_dec), 2)))

m_a <- glm(cbind(rec, non) ~ jul_dec,           data = dd, family = quasibinomial)
m_b <- glm(cbind(rec, non) ~ apres17 + jul_dec, data = dd, family = quasibinomial)
print(round(summary(m_a)$coefficients, 3)); print(round(summary(m_b)$coefficients, 3))

print(dd %>% filter(jul_dec >= 0.9, jul_dec <= 1.6) %>%
        select(annee, jul_dec, taux, apres17) %>% arrange(jul_dec) %>%
        mutate(across(c(jul_dec, taux), ~ round(.x, 3))))


################################ Test à l'echelle de l'année
dd <- dd %>% mutate(lg = qlogis(taux))

m2 <- lm(lg ~ apres17 + jul_dec, data = dd)
print(round(summary(m2)$coefficients, 3))
cat("R2 sans / avec la periode :", round(summary(lm(lg ~ jul_dec, data = dd))$r.squared, 2),
    round(summary(m2)$r.squared, 2), "\n")

# Balayage de la date de rupture, intensite en controle
sc <- do.call(rbind, lapply(2008:2020, function(c0) {
  m <- lm(lg ~ I(annee >= c0) + jul_dec, data = dd)
  data.frame(rupture = c0, coef = round(coef(m)[2], 2), rss = round(sum(resid(m)^2), 2))
}))
print(sc, row.names = FALSE)

# Stabilite : on retire chaque annee
loo <- sapply(dd$annee, function(y) coef(lm(lg ~ apres17 + jul_dec, data = dd[dd$annee != y, ]))["apres17TRUE"])
cat("Coefficient de la periode en retirant chaque annee : min", round(min(loo), 2),
    "| max", round(max(loo), 2), "| annee la plus influente :", dd$annee[which.min(loo)], "\n")

###############################################################################
# 1) Enregistrer le modele de reconnaissance (niveau annuel)
res_m2 <- list(coef = coef(m2), sigma = summary(m2)$sigma, df = m2$df.residual,
               vcov = vcov(m2), n = nrow(dd), fenetre = "jul_dec",
               definition = "logit(taux hors annuelle) ~ apres17 + anomalie SWI jul-dec (moyenne nationale des maxima communaux)")
saveRDS(res_m2, file.path(dir_out, "modele_reconnaissance_annuel.rds"))
saveRDS(data.frame(annee = 2000:2025, HZ)[, c("annee", "jul_dec", "avr_sep", "T3")],
        file.path(dir_out, "indice_alea_national_brut.rds"))
cat("sigma residuel (logit) :", round(summary(m2)$sigma, 2), "| df :", m2$df.residual, "\n")
print(round(summary(m2)$sigma, 3))

# 2) Predit vs observe, par annee
dd$pred <- plogis(fitted(m2))
print(dd %>% select(annee, jul_dec, taux, pred) %>% mutate(across(c(jul_dec, taux, pred), ~ round(.x, 3))), n = 25)

saveRDS(d2, file.path(dir_out, "table_annuelle_alea_reconnaissance.rds"))
saveRDS(dd, file.path(dir_out, "table_annuelle_reconnaissance_modele.rds"))
list.files(dir_out)

###################################################
da <- d2 %>% select(annee, a_pond, rho_pond, jul_dec, T3) %>% filter(annee >= 2003) %>%
  mutate(apres17 = annee >= 2017, lga = qlogis(a_pond))

ma1 <- lm(lga ~ jul_dec, data = da)
ma2 <- lm(lga ~ apres17 + jul_dec, data = da)
print(round(summary(ma1)$coefficients, 3)); print(round(summary(ma2)$coefficients, 3))
cat("R2 :", round(summary(ma1)$r.squared, 2), round(summary(ma2)$r.squared, 2),
    "| sigma :", round(summary(ma1)$sigma, 2), "\n")

da$pred <- plogis(fitted(ma1))
print(da %>% select(annee, jul_dec, a_pond, pred) %>%
        mutate(across(c(jul_dec, a_pond, pred), ~ round(.x, 3))), n = 25)

###########################################################################
# Variante 1 : part des communes dont l'anomalie de juil-dec depasse un seuil (largeur de la secheresse)
part_sec <- function(seuil) sapply(2000:2025, function(y) {
  cols <- which(an == y & mm %in% 7:12)
  mean(apply(Z[, cols, drop = FALSE], 1, max) > seuil)
})
P <- data.frame(annee = 2000:2025, p1 = part_sec(1), p15 = part_sec(1.5), p2 = part_sec(2))
print(round(P, 3), row.names = FALSE)

# Variante 2 : les deux indices dans le meme modele des demandes
da2 <- da %>% left_join(P, by = "annee")
for (v in c("jul_dec", "p1", "p15", "p2")) {
  m <- lm(reformulate(c("apres17", v), response = "lga"), data = da2)
  cat(v, ": R2", round(summary(m)$r.squared, 2), "| coef indice", round(coef(m)[3], 2),
      "| coef apres17", round(coef(m)[2], 2), "\n")
}

###################################################################
library(lubridate)
yrs <- 2002:2023
Zmax <- sapply(yrs, function(y) apply(Z[, which(an == y & mm %in% 7:12), drop = FALSE], 1, max))
dimnames(Zmax) <- list(rownames(Z), yrs)

base <- readRDS(file.path(dir_out, "base_evenements_couts.rds"))

com <- base %>% filter(!is.na(cout_1_arrete_modele)) %>%
  distinct(insee, n_dwellings, CSS, houses)

ev_c <- base %>% filter(annee %in% yrs, !(month(debut) == 1 & month(fin) == 12)) %>%
  group_by(insee, annee) %>% summarise(A = 1L, R = max(reconnu), .groups = "drop")

panel <- expand.grid(insee = com$insee, annee = yrs, stringsAsFactors = FALSE) %>%
  left_join(com, by = "insee") %>% left_join(ev_c, by = c("insee", "annee")) %>%
  mutate(A = coalesce(A, 0L),
         z = Zmax[cbind(match(insee, rownames(Zmax)), match(as.character(annee), colnames(Zmax)))],
         apres17 = annee >= 2017, CSSs = pmin(CSS, 100) / 100, lw = log(n_dwellings))
cat("Panel :", nrow(panel), "| communes :", n_distinct(panel$insee), "| z manquants :", sum(is.na(panel$z)),
    "| part de dépôts :", round(mean(panel$A), 3), "\n")

#####################################################
saveRDS(panel, file.path(dir_out, "panel_commune_annee_alea.rds"))
cat("Panel enregistre :", round(file.size(file.path(dir_out, "panel_commune_annee_alea.rds")) / 1e6, 1), "Mo\n")

m1 <- glm(A ~ z + apres17 + lw + CSSs, data = panel, family = binomial)
m2 <- glm(A ~ z + factor(annee) + lw + CSSs, data = panel, family = binomial)
print(round(summary(m1)$coefficients[c("z", "apres17TRUE", "lw", "CSSs"), 1:2], 3))
print(round(summary(m2)$coefficients[c("z", "lw", "CSSs"), 1:2], 3))

pa <- panel %>% filter(A == 1)
r1 <- glm(R ~ z + apres17 + lw + CSSs, data = pa, family = binomial)
r2 <- glm(R ~ z + factor(annee) + lw + CSSs, data = pa, family = binomial)
cat("\nReconnaissance | demande :", nrow(pa), "lignes | taux :", round(mean(pa$R), 3), "\n")
print(round(summary(r1)$coefficients[c("z", "apres17TRUE", "lw", "CSSs"), 1:2], 3))
print(round(summary(r2)$coefficients[c("z", "lw", "CSSs"), 1:2], 3))

############################################################
panel <- panel %>% mutate(zc = cut(z, c(-Inf, 0, 1, 2, Inf), labels = c("<0", "0-1", "1-2", ">2")),
                          periode = ifelse(apres17, "2017+", "avant"))
print(panel %>% group_by(periode, zc) %>% summarise(n = n(), depot = round(mean(A), 3), .groups = "drop"))
print(panel %>% filter(A == 1) %>% group_by(periode, zc) %>%
        summarise(n = n(), reconnu = round(mean(R), 3), .groups = "drop"))

m3 <- glm(A ~ z * apres17 + lw + CSSs, data = panel, family = binomial)
r3 <- glm(R ~ z * apres17 + lw + CSSs, data = pa, family = binomial)
print(round(coef(summary(m3))[c("z", "apres17TRUE", "z:apres17TRUE"), 1:2], 3))
print(round(coef(summary(r3))[c("z", "apres17TRUE", "z:apres17TRUE"), 1:2], 3))

if (requireNamespace("sandwich", quietly = TRUE) && requireNamespace("lmtest", quietly = TRUE)) {
  print(round(lmtest::coeftest(m2, vcov = sandwich::vcovCL(m2, cluster = panel$insee))[c("z", "lw", "CSSs"), 1:2], 3))
  print(round(lmtest::coeftest(r2, vcov = sandwich::vcovCL(r2, cluster = pa$insee))[c("z", "lw", "CSSs"), 1:2], 3))
} else cat("Installer : install.packages(c('sandwich', 'lmtest'))\n")

################################################################
print(panel %>% filter(zc == ">2", !apres17) %>% group_by(annee) %>%
        summarise(n = n(), depot = round(mean(A), 3), .groups = "drop"))

ev <- base %>%
  filter(annee %in% 2002:2023, !(month(debut) == 1 & month(fin) == 12), insee %in% rownames(Z)) %>%
  mutate(i0 = (year(debut) - 2000) * 12 + month(debut), i1 = (year(fin) - 2000) * 12 + month(fin)) %>%
  filter(i0 >= 1, i1 <= 312, i1 >= i0)
ri <- match(ev$insee, rownames(Z))
ev$z_ev <- mapply(function(r, a, b) max(Z[r, a:b]), ri, ev$i0, ev$i1)
ev$z_jd <- Zmax[cbind(match(ev$insee, rownames(Zmax)), match(as.character(ev$annee), colnames(Zmax)))]
ev <- ev %>% left_join(com %>% transmute(insee, lw = log(n_dwellings), CSSs = pmin(CSS, 100) / 100), by = "insee") %>%
  mutate(apres17 = annee >= 2017) %>% filter(!is.na(lw))
cat("Evenements :", nrow(ev), "| correlation z_ev / z_jd :", round(cor(ev$z_ev, ev$z_jd), 2), "\n")

r_jd <- glm(reconnu ~ z_jd * apres17 + lw + CSSs, data = ev, family = binomial)
r_ev <- glm(reconnu ~ z_ev * apres17 + lw + CSSs, data = ev, family = binomial)
print(round(coef(summary(r_jd))[c("z_jd", "apres17TRUE", "z_jd:apres17TRUE"), 1:2], 3))
print(round(coef(summary(r_ev))[c("z_ev", "apres17TRUE", "z_ev:apres17TRUE"), 1:2], 3))
cat("Deviance :", round(deviance(r_jd)), round(deviance(r_ev)), "\n")

ev$zc <- cut(ev$z_ev, c(-Inf, 0, 1, 2, Inf), labels = c("<0", "0-1", "1-2", ">2"))
print(ev %>% group_by(apres17, zc) %>% summarise(n = n(), reconnu = round(mean(reconnu), 3), .groups = "drop"))


###########################################
ev$dur <- ev$i1 - ev$i0 + 1
print(ev %>% group_by(dur) %>% summarise(n = n(), z = round(mean(z_ev), 2), reconnu = round(mean(reconnu), 3)) %>% head(12))

r_ev2 <- glm(reconnu ~ z_ev * apres17 + dur + lw + CSSs, data = ev, family = binomial)
print(round(coef(summary(r_ev2))[c("z_ev", "apres17TRUE", "dur", "z_ev:apres17TRUE"), 1:2], 3))

# Fenetre fixe : juillet-septembre uniquement (duree constante = 3 mois)
e3 <- ev %>% filter(month(debut) == 7, month(fin) == 9)
cat("Juil-sept :", nrow(e3), "evenements\n")
r3 <- glm(reconnu ~ z_ev * apres17 + lw + CSSs, data = e3, family = binomial)
print(round(coef(summary(r3))[c("z_ev", "apres17TRUE", "z_ev:apres17TRUE"), 1:2], 3))
e3$zc <- cut(e3$z_ev, c(-Inf, 0, 1, 2, Inf), labels = c("<0", "0-1", "1-2", ">2"))
print(e3 %>% group_by(apres17, zc) %>% summarise(n = n(), reconnu = round(mean(reconnu), 3), .groups = "drop"))

#####################################
e3 <- e3 %>% mutate(annee = year(debut))
print(e3 %>% filter(z_ev >= 0.5, z_ev < 1.5) %>% group_by(annee) %>%
        summarise(n = n(), z = round(mean(z_ev), 2), reconnu = round(mean(reconnu), 3), .groups = "drop"), n = 30)

ev <- ev %>% mutate(aligne = month(debut) %in% c(1, 4, 7, 10) & month(fin) %in% c(3, 6, 9, 12))
print(ev %>% group_by(aligne) %>% summarise(n = n(), z = round(mean(z_ev), 2), reconnu = round(mean(reconnu), 3)))
ea <- ev %>% filter(aligne)
ra <- glm(reconnu ~ z_ev * apres17 + factor(pmin(dur, 12)) + lw + CSSs, data = ea, family = binomial)
print(round(coef(summary(ra))[c("z_ev", "apres17TRUE", "z_ev:apres17TRUE"), 1:2], 3))

##############################################
yrs_all <- 2000:2025
Zj <- sapply(yrs_all, function(y) apply(Z[, which(an == y & mm %in% 7:12), drop = FALSE], 1, max))
dimnames(Zj) <- list(rownames(Z), yrs_all)

w <- base %>% filter(!is.na(cout_1_arrete_modele), insee %in% rownames(Z)) %>%
  distinct(insee, cout_1_arrete_modele)
w$omega <- w$cout_1_arrete_modele / sum(w$cout_1_arrete_modele)
Zw <- Zj[w$insee, , drop = FALSE]

H_indice <- data.frame(annee = yrs_all,
                       H_simple  = colMeans(Zj),            # 34 741 communes, sans ponderation
                       H_pond    = colSums(Zw * w$omega),   # ponderee par le cout attendu d'un arret
                       part_sup1 = colMeans(Zj > 1))        # part des communes avec anomalie > 1

cat("Communes ponderees :", nrow(w), "| controle H_simple vs jul_dec de HZ :",
    signif(max(abs(H_indice$H_simple - HZ[, "jul_dec"])), 2), "\n")
cat("Correlation de rang H_simple / H_pond :",
    round(cor(H_indice$H_simple, H_indice$H_pond, method = "spearman"), 3), "\n")
print(round(H_indice, 3), row.names = FALSE)
saveRDS(H_indice, file.path(dir_out, "H_indice_alea_national.rds"))

########################################
H_indice$H_s15 <- colMeans(Zw)     # moyenne simple sur les memes 15 764 communes que la ponderation
cat("Ecart moyen |H_simple - H_s15| :", round(mean(abs(H_indice$H_simple - H_indice$H_s15)), 3),
    "| |H_s15 - H_pond| :", round(mean(abs(H_indice$H_s15 - H_indice$H_pond)), 3), "\n")
print(round(H_indice[H_indice$annee %in% c(2003, 2007, 2012, 2016, 2018, 2022),
                     c("annee", "H_simple", "H_s15", "H_pond")], 3), row.names = FALSE)

h <- H_indice$H_pond[H_indice$annee <= 2024]      # reference 2000-2024 ; 2025 hors reference
cat("\nH_pond 2000-2024 : moyenne", round(mean(h), 2), "| ecart-type", round(sd(h), 2),
    "| asymetrie", round(mean((h - mean(h))^3) / sd(h)^3, 2), "\n")
print(round(quantile(h, c(0, .1, .25, .5, .75, .9, 1)), 2))
print(shapiro.test(h))
cat("Tendance (rang) H_pond vs annee :", round(cor(2000:2024, h, method = "spearman"), 2), "\n")
print(round(summary(lm(h ~ I(2000:2024)))$coefficients, 3))


#######################################################
part_sup <- function(delta, seuil = 1) {
  sapply(as.character(2000:2024), function(y) mean(Zj[, y] + delta > seuil))
}
tab_d <- sapply(c(0, 0.25, 0.5, 1), function(dl) part_sup(dl))
colnames(tab_d) <- paste0("delta_", c(0, 0.25, 0.5, 1))
print(round(tab_d[as.character(c(2003, 2013, 2018, 2022, 2024)), ], 3))
cat("Moyenne sur 25 ans :", round(colMeans(tab_d), 3), "\n")
cat("H_simple, moyenne 2000-2024 selon delta :", round(mean(colMeans(Zj[, as.character(2000:2024)])) + c(0, .25, .5, 1), 2), "\n")

####################################################### aggreagtion la reconnaisance de la commune à la nation 
m_dep <- glm(A ~ z * apres17 + lw + CSSs, data = panel, family = binomial)                    # depot
r_rec <- glm(R ~ z * apres17 + lw + CSSs, data = panel[panel$A == 1, ], family = binomial)    # reconnaissance | depot
print(round(coef(summary(m_dep))[c("z", "apres17TRUE", "z:apres17TRUE"), 1], 3))   # attendu : 1.242 0.016 0.766
print(round(coef(summary(r_rec))[c("z", "apres17TRUE", "z:apres17TRUE"), 1], 3))   # attendu : 0.631 -3.111 3.517

omega <- setNames(w$omega, w$insee)
panel$R0 <- coalesce(panel$R, 0L)
panel$om <- omega[panel$insee]

rho_year <- function(delta = 0, regime = "obs") {
  nd <- panel
  nd$z <- nd$z + delta
  if (regime == "avant") nd$apres17 <- FALSE
  if (regime == "apres") nd$apres17 <- TRUE
  p <- predict(m_dep, newdata = nd, type = "response") * predict(r_rec, newdata = nd, type = "response")
  tapply(nd$om * p, nd$annee, sum)
}

obs   <- tapply(panel$om * panel$A * panel$R0, panel$annee, sum)
pred0 <- rho_year(0, "obs")
cmp <- data.frame(annee = names(obs), obs = round(obs, 3), pred = round(pred0, 3))
print(cmp, row.names = FALSE)
cat("Correlation de rang :", round(cor(obs, pred0, method = "spearman"), 2),
    "| RMSE :", round(sqrt(mean((obs - pred0)^2)), 3), "| ecart-type de obs :", round(sd(obs), 3), "\n")

sc <- do.call(rbind, lapply(c(0, 0.25, 0.5, 1), function(dl) rbind(
  data.frame(delta = dl, regime = "avant 2017", rho = rho_year(dl, "avant")),
  data.frame(delta = dl, regime = "apres 2017", rho = rho_year(dl, "apres")))))
print(sc %>% group_by(delta, regime) %>%
        summarise(moyenne = round(mean(rho), 3), mediane = round(median(rho), 3), max = round(max(rho), 3),
                  .groups = "drop"))

################################# enregistrement 
saveRDS(list(m_dep = m_dep, r_rec = r_rec, omega = omega,
             note = "A ~ z*apres17 + lw + CSSs (depot) ; R|A ~ z*apres17 + lw + CSSs (reconnaissance) ; z = anomalie SWI jul-dec ; omega = poids de cout attendu"),
        file.path(dir_out, "modeles_communaux_depot_reconnaissance.rds"))
cat("Enregistre :", file.exists(file.path(dir_out, "modeles_communaux_depot_reconnaissance.rds")), "\n")


saveRDS(Zj, file.path(dir_out, "anomalies_communales_jul_dec_2000_2025.rds"))    # communes x annees, sens : haut = sec
saveRDS(w,  file.path(dir_out, "poids_cout_communes.rds"))                       # insee, cout_1_arrete_modele, omega
saveRDS(list(m_dep = m_dep, r_rec = r_rec, omega = omega,
             note = "A ~ z*apres17 + lw + CSSs (depot) ; R|A ~ z*apres17 + lw + CSSs (reconnaissance) ; z = anomalie SWI jul-dec ; omega = poids de cout attendu"),
        file.path(dir_out, "modeles_communaux_depot_reconnaissance.rds"))
cat(sapply(c("anomalies_communales_jul_dec_2000_2025.rds", "poids_cout_communes.rds",
             "modeles_communaux_depot_reconnaissance.rds", "panel_commune_annee_alea.rds",
             "H_indice_alea_national.rds"),
           function(f) paste(f, file.exists(file.path(dir_out, f)))), sep = "\n")


################################################### Graphe 
library(ggplot2)

# Theme publication minimaliste & epure
theme_journal <- theme_classic(base_size = 10, base_family = "serif") +
  theme(
    plot.title    = element_text(size = 11, face = "bold", hjust = 0, margin = margin(b = 6)),
    plot.caption  = element_text(size = 8, color = "grey30", hjust = 0, margin = margin(t = 8)),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    axis.title    = element_text(size = 9),
    axis.text     = element_text(size = 8, color = "black"),
    legend.position = "none",
    plot.margin   = margin(8, 12, 6, 8)
  )

# Data prep
H_fig <- H_indice %>%
  filter(annee <= 2024) %>%
  mutate(highlight = H_pond >= quantile(H_pond, 0.75))

top_years <- H_fig %>% slice_max(H_pond, n = 3) %>% pull(annee)

# Plot
p1 <- ggplot(H_fig, aes(x = annee, y = H_pond, fill = highlight)) +
  geom_col(width = 0.7) +
  geom_text(
    data = filter(H_fig, annee %in% top_years),
    aes(label = annee), vjust = -0.5, size = 2.8,
    family = "serif", fontface = "bold", color = "#B2182B"
  ) +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#B2182B")) +
  scale_x_continuous(breaks = seq(2000, 2024, by = 4), expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) + # Evite les labels coupes au sommet
  labs(
    x = NULL, 
    y = "Hazard index (weighted SWI anomaly, Jul.–Dec.)",
    title = "National drought hazard index, France, 2000–2024",
    caption = paste(
      "Notes: Standardized negative soil-wetness-index (SWI) anomaly (July–December max), weighted by expected per-decree cost.",
      "Reference period: 2000–2024. Top quartile highlighted; top 3 years labelled.", sep = "\n"
    )
  ) +
  theme_journal

print(p1)
ggsave(file.path(dir_out, "fig1_hazard_index.png"), p1, width = 6.5, height = 3.8, dpi = 300)

########################################################################################################## figure 2 
library(ggplot2)

# ------------------------------------------------------------
# Figure 2 -- Recognition rate by hazard intensity and regime
# ------------------------------------------------------------
dd_fig <- dd %>%
  mutate(regime = factor(
    ifelse(apres17, "2017–present", "Before 2017"),
    levels = c("Before 2017", "2017–present")
  ))

p2 <- ggplot(dd_fig, aes(x = jul_dec, y = taux, color = regime, shape = regime)) +
  geom_point(size = 2.2, alpha = 0.85) +
  geom_smooth(
    method = "lm", 
    se = TRUE, 
    linewidth = 0.5, 
    alpha = 0.12, 
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("Before 2017" = "grey40", "2017–present" = "#2166AC")) +
  scale_shape_manual(values = c("Before 2017" = 16, "2017–present" = 17)) +
  scale_y_continuous(
    limits = c(0, 1), 
    breaks = seq(0, 1, 0.25),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x = "Hazard index (weighted SWI anomaly, Jul.–Dec.)",
    y = "Recognition rate (non-annual claims)",
    title = "CatNat drought recognition rate by hazard intensity and regime",
    color = NULL, 
    shape = NULL,
    caption = paste(
      "Notes: Annual recognition rate among non-annual-window drought claims against national hazard intensity.",
      "Lines show separate linear fits by regime (regime shift coef. = 1.54, s.e. = 0.35, R² = 0.70; see Section X).",
      sep = "\n"
    )
  ) +
  theme_journal +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.margin = margin(b = 4, t = 0),
    legend.key.height = unit(0.8, "lines")
  )

print(p2)
ggsave(file.path(dir_out, "fig2_recognition_regime.png"), p2, width = 6.5, height = 4.2, dpi = 300)


####################################################################
####################################################################
# Figure: Structural Break Scan (Journal Ready - High Contrast)
####################################################################
library(ggplot2)

# Donnees
sc_fig <- data.frame(
  rupture = 2008:2020,
  coef    = c(0.23, 0.08, 0.29, 0.56, 0.70, 0.54, 0.75, 0.85, 1.37, 1.54, 1.45, 1.33, 1.34),
  rss     = c(20.45, 20.60, 20.26, 19.08, 18.21, 19.11, 17.76, 17.30, 11.84, 10.03, 11.87, 13.90, 14.81)
)

# Palette Journal Elegante
col_main  <- "#1E293B" # Bleu nuit
col_break <- "#DC2626" # Rouge vif académique
col_bg    <- "#F8FAFC" # Gris très clair pour le fond de zone

theme_journal_modern <- theme_minimal(base_size = 10, base_family = "serif") +
  theme(
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.line          = element_line(color = "black", linewidth = 0.5),
    axis.ticks         = element_line(color = "black", linewidth = 0.4),
    axis.title.x       = element_text(size = 9.5, face = "bold", margin = margin(t = 8)),
    axis.title.y       = element_text(size = 9.5, face = "bold", margin = margin(r = 8)),
    axis.text          = element_text(size = 9, color = "black"),
    plot.margin        = margin(10, 15, 10, 10)
  )

p_rupture <- ggplot(sc_fig, aes(x = rupture, y = rss)) +
  # 1. Zone d'ombre légère sur la période post-rupture (2017-2020)
  annotate("rect", xmin = 2016.5, xmax = 2020.5, ymin = -Inf, ymax = Inf,
           fill = "#F1F5F9", alpha = 0.7) +
  
  # 2. Ligne verticale rouge pour marquer la rupture
  geom_vline(xintercept = 2017, linetype = "dashed", color = col_break, linewidth = 0.6) +
  
  # 3. Courbe principale
  geom_line(color = col_main, linewidth = 0.8) +
  
  # 4. Points normaux
  geom_point(size = 2.2, shape = 21, fill = "white", color = col_main, stroke = 0.8) +
  
  # 5. Point minimum (Rouge vif)
  geom_point(data = subset(sc_fig, rupture == 2017), 
             size = 3.2, shape = 21, fill = col_break, color = "black", stroke = 0.8) +
  
  # 6. Label élégant et dégagé (à gauche de la ligne)
  annotate("label", x = 2015.2, y = 10.03, 
           label = "Optimal Break\n(2017, RSS = 10.03)", 
           family = "serif", size = 3, fontface = "bold",
           color = col_break, fill = "white", label.size = 0.3) +
  
  # Configuration des axes
  scale_x_continuous(breaks = seq(2008, 2020, by = 2), expand = c(0.03, 0.03)) +
  scale_y_continuous(limits = c(8, 22), breaks = seq(8, 22, by = 2)) +
  labs(
    x = "Candidate break year",
    y = "Residual Sum of Squares (RSS)"
  ) +
  theme_journal_modern

print(p_rupture)

# Export aux dimensions COMPACTES (4.2 x 3.0 inches) -> Évite l'effet aplati
ggsave(
  filename = file.path(dir_out, "fig_rupture_scan.pdf"),
  plot     = p_rupture,
  width    = 4.5, 
  height   = 3.2, 
  units    = "in",
  dpi      = 300
)

#########################################################################
##################################################################### HARMONIZATION
library(ggplot2)

# ------------------------------------------------------------------------------
# 0. Thème commun et dimensions unifiées
# ------------------------------------------------------------------------------
theme_journal_unified <- theme_bw(base_size = 10, base_family = "serif") +
  theme(
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title.x     = element_text(size = 9.5, face = "bold", margin = margin(t = 6)),
    axis.title.y     = element_text(size = 9.5, face = "bold", margin = margin(r = 6)),
    axis.text        = element_text(size = 8.5, color = "black"),
    axis.ticks       = element_line(color = "black", linewidth = 0.4),
    legend.position  = "top",
    legend.justification = "left",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 8.5),
    plot.margin      = margin(8, 12, 6, 8)
  )

# Palette unifiée pour l'ensemble du papier
col_blue  <- "#2166AC"
col_red   <- "#B2182B"
col_dark  <- "#1E293B"
w_fig     <- 6.5
h_fig     <- 4.0

# ------------------------------------------------------------------------------
# Figure 1 : National Drought Hazard Index
# ------------------------------------------------------------------------------
H_fig <- H_indice %>%
  filter(annee <= 2024) %>%
  mutate(highlight = H_pond >= quantile(H_pond, 0.75))

top_years <- H_fig %>% slice_max(H_pond, n = 3) %>% pull(annee)

p1 <- ggplot(H_fig, aes(x = annee, y = H_pond, fill = highlight)) +
  geom_col(width = 0.72) +
  geom_text(
    data = filter(H_fig, annee %in% top_years),
    aes(label = annee), vjust = -0.5, size = 2.8,
    family = "serif", fontface = "bold", color = col_red
  ) +
  scale_fill_manual(values = c("FALSE" = "grey75", "TRUE" = col_red)) +
  scale_x_continuous(breaks = seq(2000, 2024, by = 4), expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    x = "Year",
    y = "Hazard index (weighted SWI anomaly, Jul.–Dec.)"
  ) +
  theme_journal_unified +
  theme(legend.position = "none")

print(p1)
ggsave(file.path(dir_out, "fig1_hazard_index.png"), p1, width = w_fig, height = h_fig, dpi = 300)

# ------------------------------------------------------------------------------
# Figure 2 : Recognition Rate by Hazard Intensity and Regime
# ------------------------------------------------------------------------------
dd_fig <- dd %>%
  mutate(regime = factor(
    ifelse(apres17, "2017–present", "Before 2017"),
    levels = c("Before 2017", "2017–present")
  ))

p2 <- ggplot(dd_fig, aes(x = jul_dec, y = taux, color = regime, shape = regime)) +
  geom_point(size = 2.3, alpha = 0.85) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.6,
    alpha = 0.12,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("Before 2017" = "grey40", "2017–present" = col_blue)) +
  scale_shape_manual(values = c("Before 2017" = 16, "2017–present" = 17)) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x = "National hazard index (SWI anomaly, Jul.–Dec.)",
    y = "Recognition rate (non-annual claims)"
  ) +
  theme_journal_unified

print(p2)
ggsave(file.path(dir_out, "fig2_recognition_regime.png"), p2, width = w_fig, height = h_fig, dpi = 300)

# ------------------------------------------------------------------------------
# Figure 3 : Structural Break Scan
# ------------------------------------------------------------------------------
sc_fig <- data.frame(
  rupture = 2008:2020,
  coef    = c(0.23, 0.08, 0.29, 0.56, 0.70, 0.54, 0.75, 0.85, 1.37, 1.54, 1.45, 1.33, 1.34),
  rss     = c(20.45, 20.60, 20.26, 19.08, 18.21, 19.11, 17.76, 17.30, 11.84, 10.03, 11.87, 13.90, 14.81)
)

p3 <- ggplot(sc_fig, aes(x = rupture, y = rss)) +
  annotate("rect", xmin = 2016.5, xmax = 2020.5, ymin = -Inf, ymax = Inf,
           fill = "#F1F5F9", alpha = 0.8) +
  geom_vline(xintercept = 2017, linetype = "dashed", color = col_red, linewidth = 0.5) +
  geom_line(color = col_dark, linewidth = 0.7) +
  geom_point(size = 2.0, shape = 21, fill = "white", color = col_dark, stroke = 0.7) +
  geom_point(data = subset(sc_fig, rupture == 2017),
             size = 3.0, shape = 21, fill = col_red, color = "black", stroke = 0.7) +
  annotate("label", x = 2015.3, y = 10.03,
           label = "Optimal Break\n(2017, RSS = 10.03)",
           family = "serif", size = 2.8, fontface = "bold",
           color = col_red, fill = "white", label.size = 0.3) +
  scale_x_continuous(breaks = seq(2008, 2020, by = 2), expand = expansion(mult = c(0.03, 0.03))) +
  scale_y_continuous(limits = c(8, 22), breaks = seq(8, 22, by = 2)) +
  labs(
    x = "Candidate break year",
    y = "Residual Sum of Squares (RSS)"
  ) +
  theme_journal_unified

print(p3)
ggsave(file.path(dir_out, "fig3_rupture_scan.png"), p3, width = w_fig, height = h_fig, dpi = 300)

#########################
library(patchwork)

# (b) : remplacer l'étiquette (6e couche de p3) par une étiquette déplacée
p3$layers[[6]] <- NULL
p3 <- p3 +
  annotate("label", x = 2014.6, y = 9.2,
           label = "Minimum RSS\n(2017, RSS = 10.03)",
           family = "serif", size = 2.8, fontface = "bold",
           color = col_red, fill = "white")

# (a) : légende à l'intérieur, titres d'axes courts
p2c <- p2 +
  labs(x = "National hazard index",
       y = "Recognition rate") +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", color = NA))

fig_rupture <- (p2c | p3) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(family = "serif", face = "bold", size = 10))

print(fig_rupture)
ggsave(file.path(dir_out, "fig_rupture_combined.png"), fig_rupture,
       width = 6.5, height = 3.0, dpi = 300)

############################################################ Accessoires
print(round(coef(summary(m_dep)), 3))
print(round(coef(summary(r_rec)), 3))

cat("N demande :", nobs(m_dep), "| N reconnaissance :", nobs(r_rec), "\n")
cat("Log-vraisemblance :", round(as.numeric(logLik(m_dep))), round(as.numeric(logLik(r_rec))), "\n")
cat("AIC :", round(AIC(m_dep)), round(AIC(r_rec)), "\n")
