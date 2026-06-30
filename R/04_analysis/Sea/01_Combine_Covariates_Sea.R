


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3: WORKING DIRECTORY
# ══════════════════════════════════════════════════════════════════════════════

# ── Working directory ─────────────────────────────────────────────────────────
# Automatically detects the repository root based on the script location.
# No hardcoded path — works on any machine.
current_file <- rstudioapi::getActiveDocumentContext()$path
path_parts   <- unlist(strsplit(normalizePath(current_file, winslash = "/"), "/"))
target_index <- which(path_parts == "phd-chapter3-Atlantic-Salmon")
base_path    <- paste(path_parts[1:target_index], collapse = "/")
setwd(base_path)


# |========================|
# |     LOAD LIBRARIES     |
# |========================|

library(fst)
library(dplyr)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(car)
library(ggpubr)


# ─────────────────────────────────────────────
# 1. DeltaT_SST_River
#    Garder : river, cohort (= year), delta_T_sst_minus_river
#    Retirer les NA dans delta_T_sst_minus_river
# ─────────────────────────────────────────────

deltaT_stj_1sw <- read_fst(file.path(base_path, "R/03_clean_data/Sea/DeltaT_SST_River_stj_1SW.fst")) %>%
  dplyr::select(river, cohort = year, delta_T_sst_minus_river) %>%
  filter(!is.na(delta_T_sst_minus_river))

deltaT_stj_2sw <- read_fst(file.path(base_path, "R/03_clean_data/Sea/DeltaT_SST_River_stj_2SW.fst")) %>%
  dplyr::select(river, cohort = year, delta_T_sst_minus_river) %>%
  filter(!is.na(delta_T_sst_minus_river))

deltaT_tri_1sw <- read_fst(file.path(base_path, "R/03_clean_data/Sea/DeltaT_SST_River_tri_1SW.fst")) %>%
  dplyr::select(river, cohort = year, delta_T_sst_minus_river) %>%
  filter(!is.na(delta_T_sst_minus_river))

deltaT_tri_2sw <- read_fst(file.path(base_path, "R/03_clean_data/Sea/DeltaT_SST_River_tri_2SW.fst")) %>%
  dplyr::select(river, cohort = year, delta_T_sst_minus_river) %>%
  filter(!is.na(delta_T_sst_minus_river))


# ─────────────────────────────────────────────
# 4. plankton_index
#    Garder : toutes les colonnes sauf sw_type
#    Retirer les NA dans plankton_index
# ─────────────────────────────────────────────

food_stj_1sw <- read_fst(file.path(base_path, "R/03_clean_data/Sea/Food_Index_stj_1SW.fst")) %>%
  dplyr::select(-sw_type) %>%
  filter(!is.na(plankton_index))

food_stj_2sw <- read_fst(file.path(base_path, "R/03_clean_data/Sea/Food_Index_stj_2SW.fst")) %>%
  dplyr::select(-sw_type) %>%
  filter(!is.na(plankton_index), !is.na(plankton_index))

food_tri_1sw <- read_fst(file.path(base_path, "R/03_clean_data/Sea/Food_Index_tri_1SW.fst")) %>%
  dplyr::select(-sw_type) %>%
  filter(!is.na(plankton_index))

food_tri_2sw <- read_fst(file.path(base_path, "R/03_clean_data/Sea/Food_Index_tri_2SW.fst")) %>%
  dplyr::select(-sw_type) %>%
  filter(!is.na(plankton_index), !is.na(plankton_index))


##############################################################################################################################################################################

# ─────────────────────────────────────────────
# 5. TH_annual
#    Garder : toutes les colonnes + ajouter river
#    Retirer les NA dans toutes les colonnes sauf cohort
# ─────────────────────────────────────────────

th_stj_1sw <- read_fst(file.path(base_path, "R/03_clean_data/Sea/TH_annual_stj_1sw.fst")) %>%
  mutate(river = "stj") %>%
  filter(if_all(-cohort, ~ !is.na(.)))

th_stj_2sw <- read_fst(file.path(base_path, "R/03_clean_data/Sea/TH_annual_stj_2sw.fst")) %>%
  mutate(river = "stj") %>%
  filter(if_all(-cohort, ~ !is.na(.)))

th_tri_1sw <- read_fst(file.path(base_path, "R/03_clean_data/Sea/TH_annual_tri_1sw.fst")) %>%
  mutate(river = "tri") %>%
  filter(if_all(-cohort, ~ !is.na(.)))

th_tri_2sw <- read_fst(file.path(base_path, "R/03_clean_data/Sea/TH_annual_tri_2sw.fst")) %>%
  mutate(river = "tri") %>%
  filter(if_all(-cohort, ~ !is.na(.)))


# ─────────────────────────────────────────────
# 6. Smolt_Total_Length
#    Garder : toutes les colonnes, year -> cohort
#    Retirer les NA dans smolt_tl
#    (même jeu pour stj et tri, on filtre par river au moment du join)
# ─────────────────────────────────────────────

smolt_tl <- read_fst(file.path(base_path, "03_Clean_Data/Sea/Smolt_Total_Length.fst")) %>%
  rename(cohort = year) %>%
  filter(!is.na(smolt_tl))


# ─────────────────────────────────────────────
# 7. Smolt_Weight
#    Garder : toutes les colonnes, year -> cohort
#    Retirer les NA dans smolt_weight
# ─────────────────────────────────────────────

smolt_weight <- read_fst(file.path(base_path, "03_Clean_Data/Sea/Smolt_Weight.fst")) %>%
  rename(cohort = year) %>%
  filter(!is.na(smolt_weight))


# ─────────────────────────────────────────────
# 8. Smolt_Fulton_K
#    Garder : toutes les colonnes, year -> cohort
#    Retirer les NA dans smolt_k
# ─────────────────────────────────────────────

smolt_k <- read_fst(file.path(base_path, "03_Clean_Data/Sea/Smolt_Fulton_K.fst")) %>%
  rename(cohort = year) %>%
  filter(!is.na(smolt_k))


# ─────────────────────────────────────────────
# 9. ExploitationRate_Sea
#    Garder : toutes les colonnes 
#    Retirer les NA dans er1NAC_lag1, er_WG_lag1 ou er2NAC_lag2
#    (même jeu pour stj et tri, on joint par cohort)
# ─────────────────────────────────────────────

exploit <- read_fst(file.path(base_path, "03_Clean_Data/Sea/ExploitationRate_Sea.fst")) %>%
  filter(!is.na(er1NAC_lag1) & !is.na(er_WG_lag1) & !is.na(er2NAC_lag2))


# ─────────────────────────────────────────────
# 9. Density_Lab
#    Garder : toutes les colonnes 
#    (même jeu pour stj et tri, on joint par cohort)
# ─────────────────────────────────────────────

density_lab <- read_fst(file.path(base_path, "03_Clean_Data/Sea/Density_Lab.fst"))


# ─────────────────────────────────────────────
# 10. SalmonSurvival_Sea  (variable réponse — aucune ligne retirée)
#     river : "Saint-Jean" -> "stj", "Trinite" -> "tri"
#     year -> cohort
# ─────────────────────────────────────────────

survival <- read_fst(file.path(base_path, "03_Clean_Data/Sea/SalmonSurvival_Sea.fst")) %>%
  mutate(river = ifelse(river == "Saint-Jean", "stj", "tri")) %>%
  rename(cohort = year)


# ═════════════════════════════════════════════
# COMBINAISON FINALE : 4 jeux de données
# full_join sur cohort (et river quand pertinent)
# pour conserver toutes les années disponibles
# ═════════════════════════════════════════════

combine_data <- function(river_id, sw_label,
                         deltaT, capelin, plankton, food, th,
                         survival_data, exploit_data,
                         smolt_tl_data, smolt_weight_data, 
                         smolt_k_data, density_data) {
  
  # Survie (variable réponse) : on commence par là pour avoir toutes les années
  base <- survival_data %>% filter(river == river_id)
  
  # Variables communes aux deux rivières (exploit, smolts) -> filtrées par river
  ex  <- exploit_data                                  # pas de colonne river dans le fichier source
  stl <- smolt_tl_data     %>% filter(river == river_id)
  sw  <- smolt_weight_data %>% filter(river == river_id)
  sk  <- smolt_k_data      %>% filter(river == river_id)
  
  base %>%
    full_join(deltaT,   by = c("cohort", "river")) %>%
    full_join(capelin,  by = c("cohort", "river")) %>%
    full_join(plankton, by = c("cohort", "river")) %>%
    full_join(food,     by = c("cohort", "river")) %>%
    full_join(th,       by = c("cohort", "river")) %>%
    full_join(stl,      by = c("cohort", "river")) %>%
    full_join(sw,       by = c("cohort", "river")) %>%
    full_join(sk,       by = c("cohort", "river")) %>%
    full_join(ex,       by = "cohort") %>%
    full_join(density_data, by = "cohort") %>% 
    arrange(cohort)
}

data_stj_1sw <- combine_data(
  river_id = "stj", sw_label = "1sw",
  deltaT   = deltaT_stj_1sw,
  capelin  = capelin_stj_1sw,
  plankton = plankton_stj_1sw,
  food     = food_stj_1sw,
  th       = th_stj_1sw,
  survival_data    = survival %>% dplyr::select(cohort, river, return_1sw),
  exploit_data     = exploit,
  smolt_tl_data    = smolt_tl,
  smolt_weight_data = smolt_weight,
  smolt_k_data     = smolt_k,
  density_data      = density_lab %>% dplyr::select(YEAR, pfa_1sw) %>% rename(cohort = YEAR)
) 

data_stj_2sw <- combine_data(
  river_id = "stj", sw_label = "2sw",
  deltaT   = deltaT_stj_2sw,
  capelin  = capelin_stj_2sw,
  plankton = plankton_stj_2sw,
  food     = food_stj_2sw,
  th       = th_stj_2sw,
  survival_data    = survival %>% dplyr::select(cohort, river, return_2sw),
  exploit_data     = exploit,
  smolt_tl_data    = smolt_tl,
  smolt_weight_data = smolt_weight,
  smolt_k_data     = smolt_k,
  density_data      = density_lab %>% dplyr::select(YEAR, pfa_2sw) %>% dplyr::rename(cohort = YEAR)
) 

data_tri_1sw <- combine_data(
  river_id = "tri", sw_label = "1sw",
  deltaT   = deltaT_tri_1sw,
  capelin  = capelin_tri_1sw,
  plankton = plankton_tri_1sw,
  food     = food_tri_1sw,
  th       = th_tri_1sw,
  survival_data    = survival %>% dplyr::select(cohort, river, return_1sw),
  exploit_data     = exploit,
  smolt_tl_data    = smolt_tl,
  smolt_weight_data = smolt_weight,
  smolt_k_data     = smolt_k,
  density_data      = density_lab %>% dplyr::select(YEAR, pfa_1sw) %>% rename(cohort = YEAR)
) 

data_tri_2sw <- combine_data(
  river_id = "tri", sw_label = "2sw",
  deltaT   = deltaT_tri_2sw,
  capelin  = capelin_tri_2sw,
  plankton = plankton_tri_2sw,
  food     = food_tri_2sw,
  th       = th_tri_2sw,
  survival_data    = survival %>% dplyr::select(cohort, river, return_2sw),
  exploit_data     = exploit,
  smolt_tl_data    = smolt_tl,
  smolt_weight_data = smolt_weight,
  smolt_k_data     = smolt_k,
  density_data      = density_lab %>% dplyr::select(YEAR, pfa_2sw) %>% rename(cohort = YEAR)
) 


# ══════════════════════════════════════════════════════════════════════════════
# CRÉATION D'UN INDICE SYNTHÉTIQUE POUR LES TAUX D'EXPLOITATION (2SW)
# ══════════════════════════════════════════════════════════════════════════════
#
# er_WG_lag1 et er2NAC_lag2 mesurent le même phénomène (pression de pêche)
# dans deux zones géographiques différentes mais biologiquement liées pour
# les saumons 2SW. Un indice synthétique tiré de l'axe 1 d'une PCA est
# construit pour résoudre leur collinéarité mutuelle tout en conservant
# l'information biologique des deux variables.
#
# Étapes :
#   1. Corrélation entre les deux variables (brute et log)
#   2. PCA sur échelle log
#   3. Vérification de l'interprétabilité biologique du score
#   4. Ajout de er_index dans les jeux de données finaux
#
# ══════════════════════════════════════════════════════════════════════════════


# ── 1. Corrélation entre er_WG_lag1 et er2NAC_lag2 ───────────────────────────
#
# Condition nécessaire : les deux variables doivent être fortement et
# positivement corrélées pour que leur variance commune soit biologiquement
# interprétable sur un seul axe.
#
# Règle d'interprétation :
#   r ≥ 0.70  → corrélation suffisante, indice PCA justifié
#   r < 0.70  → les deux variables captent des dynamiques trop distinctes;
#               un indice commun effacerait de l'information importante
#
# Les taux d'exploitation sont également testés sur échelle log, car les
# distributions brutes sont fortement asymétriques. La version log est
# celle utilisée pour la PCA (voir étape 2).

# STJ
cor_stj_raw <- cor.test(data_stj_2sw$er_WG_lag1, data_stj_2sw$er2NAC_lag2, use = "complete.obs")
cor_stj_log <- cor.test(log(data_stj_2sw$er_WG_lag1), log(data_stj_2sw$er2NAC_lag2), use = "complete.obs")

# TRI
cor_tri_raw <- cor.test(data_tri_2sw$er_WG_lag1, data_tri_2sw$er2NAC_lag2, use = "complete.obs")
cor_tri_log <- cor.test(log(data_tri_2sw$er_WG_lag1), log(data_tri_2sw$er2NAC_lag2), use = "complete.obs")

# Résultats — corrélations brutes et log
sprintf("STJ — r raw = %.3f  (p = %.4f)", cor_stj_raw$estimate, cor_stj_raw$p.value)
sprintf("STJ — r log = %.3f  (p = %.4f)", cor_stj_log$estimate, cor_stj_log$p.value)
sprintf("TRI — r raw = %.3f  (p = %.4f)", cor_tri_raw$estimate, cor_tri_raw$p.value)
sprintf("TRI — r log = %.3f  (p = %.4f)", cor_tri_log$estimate, cor_tri_log$p.value)

# Figures côte à côte (brute vs log) par rivière
make_cor_plot <- function(x, y, r, p, color, title, xlab, ylab) {
  label_text <- sprintf("r = %.3f\np = %.4f", r, p)
  ggplot(data.frame(x = x, y = y), aes(x = x, y = y)) +
    geom_point(color = color) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
             label = label_text, size = 4, lineheight = 0.9) +
    labs(title = title, x = xlab, y = ylab) +
    theme_minimal()
}

p_stj_raw <- make_cor_plot(
  x = data_stj_2sw$er_WG_lag1, y = data_stj_2sw$er2NAC_lag2,
  r = cor_stj_raw$estimate, p = cor_stj_raw$p.value,
  color = "steelblue", title = "STJ 2SW — Échelle brute",
  xlab = "er_WG_lag1", ylab = "er2NAC_lag2"
)
p_stj_log <- make_cor_plot(
  x = log(data_stj_2sw$er_WG_lag1), y = log(data_stj_2sw$er2NAC_lag2),
  r = cor_stj_log$estimate, p = cor_stj_log$p.value,
  color = "steelblue", title = "STJ 2SW — Échelle log",
  xlab = "log(er_WG_lag1)", ylab = "log(er2NAC_lag2)"
)
p_tri_raw <- make_cor_plot(
  x = data_tri_2sw$er_WG_lag1, y = data_tri_2sw$er2NAC_lag2,
  r = cor_tri_raw$estimate, p = cor_tri_raw$p.value,
  color = "tomato", title = "TRI 2SW — Échelle brute",
  xlab = "er_WG_lag1", ylab = "er2NAC_lag2"
)
p_tri_log <- make_cor_plot(
  x = log(data_tri_2sw$er_WG_lag1), y = log(data_tri_2sw$er2NAC_lag2),
  r = cor_tri_log$estimate, p = cor_tri_log$p.value,
  color = "tomato", title = "TRI 2SW — Échelle log",
  xlab = "log(er_WG_lag1)", ylab = "log(er2NAC_lag2)"
)

ggarrange(p_stj_raw, p_stj_log, ncol = 2, nrow = 1)
ggarrange(p_tri_raw, p_tri_log, ncol = 2, nrow = 1)


# ── 2. PCA sur er_WG_lag1 et er2NAC_lag2 (échelle log) ───────────────────────
#
# Les taux d'exploitation sont transformés en log avant la PCA pour deux
# raisons. Premièrement, les effets de la mortalité par pêche sont
# multiplicatifs : une augmentation du taux de 2% à 4% n'est pas équivalente
# à une augmentation de 20% à 22%, même si la différence absolue est
# identique. L'échelle log reflète donc mieux la nature biologique du
# phénomène. Deuxièmement, les distributions brutes sont fortement
# asymétriques, avec une masse de valeurs concentrées près de zéro et
# quelques observations extrêmes. Sans transformation, l'axe 1 serait dominé
# par ces valeurs extrêmes. La transformation log homogénéise la distribution
# et assure que chaque observation contribue de façon comparable à la
# structure de covariance.
#
# La transformation log améliore substantiellement la distribution, mais ne
# la rend pas parfaitement symétrique. Quelques points restent écartés de la
# relation linéaire principale, ce qui introduit une incertitude résiduelle
# dans le score PCA. Cet indice doit donc être interprété comme un gradient
# approximatif de pression de pêche globale, et non comme une mesure précise.
#
# Ce qu'on vérifie :
#   a) % de variance expliquée par Dim 1
#      → Objectif : ≥ 70–80 %. En dessous de 60 %, le score est trop imprécis.
#   b) Loadings des deux variables sur Dim 1
#      → Les deux doivent pointer dans le même sens (même signe) et avoir
#        une magnitude similaire.
er_vars_log_stj <- data.frame(
  er_WG_lag1  = log(data_stj_2sw$er_WG_lag1),
  er2NAC_lag2 = log(data_stj_2sw$er2NAC_lag2)
)
er_vars_log_stj <- na.omit(er_vars_log_stj)
er_vars_log_tri <- data.frame(
  er_WG_lag1  = log(data_tri_2sw$er_WG_lag1),
  er2NAC_lag2 = log(data_tri_2sw$er2NAC_lag2)
)
er_vars_log_tri <- na.omit(er_vars_log_tri)
pca_stj <- PCA(er_vars_log_stj, scale.unit = TRUE, graph = FALSE)
pca_tri <- PCA(er_vars_log_tri, scale.unit = TRUE, graph = FALSE)
# Variance expliquée par Dim 1 et loadings — STJ
sprintf("STJ — Variance expliquée par Dim 1 : %.1f %%", pca_stj$eig[1, 2])
round(pca_stj$var$coord[, 1, drop = FALSE], 3)
# Variance expliquée par Dim 1 et loadings — TRI
sprintf("TRI — Variance expliquée par Dim 1 : %.1f %%", pca_tri$eig[1, 2])
round(pca_tri$var$coord[, 1, drop = FALSE], 3)
# Biplots
fviz_pca_var(pca_stj, col.var = "steelblue", repel = TRUE,
             title = "STJ 2SW — PCA log(er_WG_lag1) / log(er2NAC_lag2)")
fviz_pca_var(pca_tri, col.var = "tomato", repel = TRUE,
             title = "TRI 2SW — PCA log(er_WG_lag1) / log(er2NAC_lag2)")
# ── 3. Interprétabilité biologique du score ───────────────────────────────────
#
# Le score PCA (coordonnée individuelle sur Dim 1) est une variable centrée
# réduite : sa moyenne est 0 et son écart-type est 1. Un score élevé doit
# correspondre à une pression de pêche élevée pour que l'indice soit
# biologiquement interprétable dans le modèle.
#
# On vérifie ceci en calculant la corrélation du score avec chacune des deux
# variables originales. Si les deux corrélations sont positives et fortes,
# le score est orienté dans le bon sens. Si elles sont négatives, le signe
# du score est automatiquement inversé.

data_stj_2sw$er_index <- NA_real_
data_tri_2sw$er_index <- NA_real_
data_stj_2sw$er_index[as.integer(rownames(er_vars_log_stj))] <- pca_stj$ind$coord[, 1]
data_tri_2sw$er_index[as.integer(rownames(er_vars_log_tri))] <- pca_tri$ind$coord[, 1]
cor_stj_wg  <- cor(data_stj_2sw$er_index, data_stj_2sw$er_WG_lag1,  use = "complete.obs")
cor_stj_nac <- cor(data_stj_2sw$er_index, data_stj_2sw$er2NAC_lag2, use = "complete.obs")
cor_tri_wg  <- cor(data_tri_2sw$er_index, data_tri_2sw$er_WG_lag1,  use = "complete.obs")
cor_tri_nac <- cor(data_tri_2sw$er_index, data_tri_2sw$er2NAC_lag2, use = "complete.obs")

# Correction automatique de l'orientation si nécessaire
if (cor_stj_wg < 0) {
  data_stj_2sw$er_index <- data_stj_2sw$er_index * -1
  cor_stj_wg            <- cor_stj_wg  * -1
  cor_stj_nac           <- cor_stj_nac * -1
  # STJ — Score inversé automatiquement (er_index * -1)
} else {
  # STJ — Score correctement orienté, aucun changement
}

if (cor_tri_wg < 0) {
  data_tri_2sw$er_index <- data_tri_2sw$er_index * -1
  cor_tri_wg            <- cor_tri_wg  * -1
  cor_tri_nac           <- cor_tri_nac * -1
  # TRI — Score inversé automatiquement (er_index * -1)
} else {
  # TRI — Score correctement orienté, aucun changement
}

# Corrélations finales du score avec les variables originales — STJ
sprintf("STJ — r(er_index, er_WG_lag1)  = %.3f", cor_stj_wg)
sprintf("STJ — r(er_index, er2NAC_lag2) = %.3f", cor_stj_nac)

# Corrélations finales du score avec les variables originales — TRI
sprintf("TRI — r(er_index, er_WG_lag1)  = %.3f", cor_tri_wg)
sprintf("TRI — r(er_index, er2NAC_lag2) = %.3f", cor_tri_nac)



# |=========================|
# |     EXPORT DATASETS     |
# |=========================|

write_fst(data_stj_1sw, path = file.path(base_path, "04_Analysis/Sea", "PreSelection_DataSet_1sw_stj.fst"), compress = 100)
write_fst(data_stj_2sw, path = file.path(base_path, "04_Analysis/Sea", "PreSelection_DataSet_2sw_stj.fst"), compress = 100)
write_fst(data_tri_1sw, path = file.path(base_path, "04_Analysis/Sea", "PreSelection_DataSet_1sw_tri.fst"), compress = 100)
write_fst(data_tri_2sw, path = file.path(base_path, "04_Analysis/Sea", "PreSelection_DataSet_2sw_tri.fst"), compress = 100)


# |=========================================|
# |     CREATE A FIGURE TO VIZUALIZE NA     |
# |=========================================|

to_long <- function(data) {
  data <- data %>%
    group_by(river) %>%
    arrange(cohort, .by_group = TRUE) %>%
    complete(cohort = seq(min(cohort), max(cohort), by = 1)) %>%
    ungroup() %>%
    dplyr::select(-river) %>%
    pivot_longer(cols = -cohort,
                 names_to = "variable",
                 values_to = "value") %>%
    mutate(presence = ifelse(is.na(value), "Missing", "Available"))
  return(data)
}

data_stj_1sw_long <- to_long(data_stj_1sw) # St. Jean River (1SW)
data_stj_2sw_long <- to_long(data_stj_2sw) # St. Jean River (2SW)
data_tri_1sw_long <- to_long(data_tri_1sw) # Trinite River (1SW)
data_tri_2sw_long <- to_long(data_tri_2sw) # Trinite River (2SW)

# Fonction pour ordonner les variables (return_*sw en premier, reste alphabétique)
order_vars <- function(long_data, sw_label) {
  response_var <- paste0("return_", sw_label)
  all_vars     <- unique(long_data$variable)
  other_vars   <- sort(setdiff(all_vars, response_var))
  c(other_vars, response_var)   # ggplot inverse l'axe y → response_var apparaît en haut
}

# Thème commun
theme_availability <- theme_minimal(base_size = 12) +
  theme(
    panel.grid       = element_blank(),
    panel.background = element_rect(fill = "#FAFAFA", color = NA),
    plot.background  = element_rect(fill = "#FAFAFA", color = NA),
    axis.text.y      = element_text(size = 8, color = "#3A3A3A"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 7.5, color = "#3A3A3A"),
    axis.title.x     = element_text(size = 10, margin = margin(t = 6), color = "#3A3A3A"),
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 13, color = "#1A1A2E"),
    legend.position  = "bottom",
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 9),
    plot.margin      = margin(10, 15, 10, 10)
  )

# Palette : vert sauge pour Available, gris très clair pour Missing
fill_values <- c("Available" = "#4A9E8E", "Missing" = "#EDEDEE")

# ── St. Jean River (1SW) ──────────────────────────────────────────────────────
var_order_stj_1sw <- order_vars(data_stj_1sw_long, "1sw")

ggplot(data_stj_1sw_long,
       aes(x = cohort,
           y = factor(variable, levels = var_order_stj_1sw),
           fill = presence)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_manual(values = fill_values, name = "Data availability") +
  scale_x_continuous(expand = c(0, 0),
                     breaks = seq(min(data_stj_1sw_long$cohort),
                                  max(data_stj_1sw_long$cohort), 1)) +
  labs(x = "Cohort", y = NULL, title = "St. Jean River (1SW)") +
  theme_availability

# ── St. Jean River (2SW) ──────────────────────────────────────────────────────
var_order_stj_2sw <- order_vars(data_stj_2sw_long, "2sw")

ggplot(data_stj_2sw_long,
       aes(x = cohort,
           y = factor(variable, levels = var_order_stj_2sw),
           fill = presence)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_manual(values = fill_values, name = "Data availability") +
  scale_x_continuous(expand = c(0, 0),
                     breaks = seq(min(data_stj_2sw_long$cohort),
                                  max(data_stj_2sw_long$cohort), 1)) +
  labs(x = "Cohort", y = NULL, title = "St. Jean River (2SW)") +
  theme_availability

# ── Trinite River (1SW) ───────────────────────────────────────────────────────
var_order_tri_1sw <- order_vars(data_tri_1sw_long, "1sw")

ggplot(data_tri_1sw_long,
       aes(x = cohort,
           y = factor(variable, levels = var_order_tri_1sw),
           fill = presence)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_manual(values = fill_values, name = "Data availability") +
  scale_x_continuous(expand = c(0, 0),
                     breaks = seq(min(data_tri_1sw_long$cohort),
                                  max(data_tri_1sw_long$cohort), 1)) +
  labs(x = "Cohort", y = NULL, title = "Trinite River (1SW)") +
  theme_availability

# ── Trinite River (2SW) ───────────────────────────────────────────────────────
var_order_tri_2sw <- order_vars(data_tri_2sw_long, "2sw")

ggplot(data_tri_2sw_long,
       aes(x = cohort,
           y = factor(variable, levels = var_order_tri_2sw),
           fill = presence)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_manual(values = fill_values, name = "Data availability") +
  scale_x_continuous(expand = c(0, 0),
                     breaks = seq(min(data_tri_2sw_long$cohort),
                                  max(data_tri_2sw_long$cohort), 1)) +
  labs(x = "Cohort", y = NULL, title = "Trinite River (2SW)") +
  theme_availability




