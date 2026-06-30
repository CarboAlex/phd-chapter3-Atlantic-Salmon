




# |=================|
# |    BASE PATH    |
# |=================|

base_path <- "C:/Users/carbo/Documents/Alex/Ecole/Doctorat/Chapitre_3"


# |========================|
# |     LOAD LIBRARIES     |
# |========================|

library(fst)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(FactoMineR)
library(ggrepel)
library(tidyverse)
library(mgcv)
library(car)
library(DHARMa)
library(betareg)


# |======================================================|
# |     LOAD PCA AND CORRELATION (HEATMAP) FUNCTIONS     |
# |======================================================|

plot_cor <- function(data, title_ = ""){
  
  # Correlation matrix
  cor_mat <- cor(
    data %>% dplyr::select(where(is.numeric), -cohort),
    use = "pairwise.complete.obs"
  )
  
  # Convert to long data.frame
  cor_df <- as.data.frame(cor_mat) %>%
    rownames_to_column("var1") %>%
    pivot_longer(-var1, names_to = "var2", values_to = "Cor")
  
  # Keep only the top triangle
  cor_df <- cor_df %>%
    mutate(
      var1 = factor(var1, levels = colnames(cor_mat)),
      var2 = factor(var2, levels = colnames(cor_mat))
    ) %>%
    filter(as.numeric(var1) < as.numeric(var2))
  
  # Plot
  p <- ggplot(cor_df, aes(var1, var2, fill = Cor)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(Cor, 2)), size = 3, color = "black") +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, limits = c(-1, 1)
    ) +
    scale_y_discrete(limits = rev) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(title = title_) +
    coord_fixed()
  
  return(p)
}


# |=======================================|
# |     CREATE A FUNCTION TO PLOT PCA     |
# |=======================================|

plot_pca <- function(data, highlight_var = "surv_egg_smolt", arrow_scale = 3, title_ = "") {
  
  data <- data %>%
    dplyr::select(where(is.numeric), -cohort)
  
  # 1. Supprimer les colonnes avec une seule valeur
  keep <- sapply(data, function(x) length(unique(x)) != 1)
  data <- data[, keep]
  
  # 2. PCA
  pca_res <- PCA(data, scale.unit = TRUE, graph = FALSE)
  
  # 3. Coordonnées des individus
  ind_df <- as.data.frame(pca_res$ind$coord)
  ind_df$ind <- rownames(ind_df)
  
  # 4. Coordonnées des variables
  var_df <- as.data.frame(pca_res$var$coord)
  var_df$var <- rownames(var_df)
  
  # 5. Couleur des variables
  var_df$col <- ifelse(var_df$var == highlight_var, "red", "black")
  
  # 6. Appliquer un facteur de scale pour rendre les flèches plus visibles
  var_df_scaled <- var_df %>%
    mutate(
      Dim.1 = Dim.1 * arrow_scale,
      Dim.2 = Dim.2 * arrow_scale
    )
  
  # 7. Plot ggplot avec ggrepel
  ggplot() +
    # individus
    geom_point(data = ind_df, aes(x = Dim.1, y = Dim.2), color = "grey30") +
    # flèches variables
    geom_segment(data = var_df_scaled,
                 aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2, color = col),
                 arrow = arrow(length = unit(0.3, "cm")), size = 0.7) +
    # labels variables avec ggrepel
    geom_text_repel(data = var_df_scaled,
                    aes(x = Dim.1, y = Dim.2, label = var, color = col),
                    size = 4,
                    box.padding = 0.8,
                    point.padding = 0.8,
                    max.overlaps = Inf,
                    segment.color = "grey50") +
    scale_color_identity() +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = paste0("Dim 1 (", round(pca_res$eig[1,2],1), "%)"),
         y = paste0("Dim 2 (", round(pca_res$eig[2,2],1), "%)"),
         title = title_) +
    coord_fixed()
}

# |=====================|
# |     IMPORT DATA     |
# |=====================|

stj_1sw <- read_fst(file.path(base_path, "04_Analysis/Sea/PreSelection_DataSet_1sw_stj.fst"))
stj_2sw <- read_fst(file.path(base_path, "04_Analysis/Sea/PreSelection_DataSet_2sw_stj.fst"))
tri_1sw <- read_fst(file.path(base_path, "04_Analysis/Sea/PreSelection_DataSet_1sw_tri.fst"))
tri_2sw <- read_fst(file.path(base_path, "04_Analysis/Sea/PreSelection_DataSet_2sw_tri.fst"))


# Make sure return is between 0 and 1 (not between 0 and 100)
if(any(stj_1sw$return_1sw > 1)){
  stj_1sw <- stj_1sw %>%
    mutate(return_1sw = return_1sw / 100)
}
if(any(stj_2sw$return_2sw > 1)){
  stj_2sw <- stj_2sw %>%
    mutate(return_2sw = return_2sw / 100)
}
if(any(tri_1sw$return_1sw > 1)){
  tri_1sw <- tri_1sw %>%
    mutate(return_1sw = return_1sw / 100)
}
if(any(tri_2sw$return_2sw > 1)){
  tri_2sw <- tri_2sw %>%
    mutate(return_2sw = return_2sw / 100)
}


# |=====================================================|
# |     REMOVE COLUMNS THAT WERE DECIDED TO NOT USE     |
# |=====================================================|

# ---------- Define vectors of columns to remove ----------

# Were used to create the food index, and only the food index variable will be used
food <- c(
  "Capelin_Original_tp1",
  "Capelin_Imputed_tp1",
  "Capelin_Original_tp2",
  "Capelin_Imputed_tp2",
  "Plankton_tp1"
)

# This is an artefact from another script and was not planned to use
polygons <- c(
  "StJean_BaieGaspe_TH",
  "Trinite_30km_TH"
)

# In smolts condition variables, there are less missing data in K and since this metric is representative of the smolt condition (maybe even more than the 2 others) only this one is kept
smolt_condition <- c(
  "smolt_tl",
  "smolt_weight"
)


# Combine all columns to remove in one vector
cols_to_remove <- c(
  food,
  polygons,
  smolt_condition
)

# ---------- Apply filtering and column removal ----------

clean_dataframe <- function(df) {
  # Keep only rows with return data
  if(any(colnames(df) == "return_1sw")){
    df <- df %>%
      dplyr::filter(!is.na(return_1sw))
  } else{
    if(any(colnames(df) == "return_2sw")){
      df <- df %>%
        dplyr::filter(!is.na(return_2sw))
    } else{
      stop("Neither of the column return_1sw or return_2sw are available in the data set.")
    }
  }
  
  df <- df %>%
    dplyr::select(where(~ !all(is.na(.)))) %>%     # Remove columns that are 100% NA
    dplyr::select(where(~ !all(. == 0))) %>%       # Remove columns that are 100% 0
    dplyr::select(-any_of(cols_to_remove)) %>%     # Remove selected columns if present
    na.omit()                                      # Remove any remaining rows with NA
  
  return(df)
}

# Apply to stj and tri (and remove the exploitation rate that should not be in the data set)
# Apply to stj and tri (and remove the exploitation rate that should not be in the data set)
stj_1sw <- clean_dataframe(stj_1sw) %>% dplyr::select(-er_WG_lag1, -er2NAC_lag2)
stj_2sw <- clean_dataframe(stj_2sw) %>% dplyr::select(-er1NAC_lag1)  
tri_1sw <- clean_dataframe(tri_1sw) %>% dplyr::select(-er_WG_lag1, -er2NAC_lag2)
tri_2sw <- clean_dataframe(tri_2sw) %>% dplyr::select(-er1NAC_lag1)  


# |====================================|
# |     CREATE GROUPS OF VARIABLES     |
# |====================================|

# Temperatures
stj_1sw_temp <- stj_1sw %>%
  dplyr::select(river, cohort, return_1sw, delta_T_sst_minus_river,  StJeanCorridor_TH, BelleIsle_TH, LabradorEntrance_TH, LabradorSea_TH) %>%
  na.omit()

stj_2sw_temp <- stj_2sw %>%
  dplyr::select(river, cohort, return_2sw, delta_T_sst_minus_river,  StJeanCorridor_TH, BelleIsle_TH, LabradorEntrance_TH, LabradorSea_TH) %>%
  na.omit()

tri_1sw_temp <- tri_1sw %>%
  dplyr::select(river, cohort, return_1sw, delta_T_sst_minus_river, TriniteCorridor_TH, BelleIsle_TH, LabradorEntrance_TH, LabradorSea_TH) %>%
  na.omit()

tri_2sw_temp <- tri_2sw %>%
  dplyr::select(river, cohort, return_2sw, delta_T_sst_minus_river, TriniteCorridor_TH, BelleIsle_TH, LabradorEntrance_TH, LabradorSea_TH) %>%
  na.omit()

# Exploitation rate (only for 2sw since there is only one variable in 1sw)
stj_2sw_exp <- stj_2sw %>%
  dplyr::select(river, cohort, return_2sw, er_WG_lag1) %>%
  na.omit()

tri_2sw_exp <- tri_2sw %>%
  dplyr::select(river, cohort, return_2sw, er_WG_lag1) %>%
  na.omit()


# |==============================|
# |     EXPLORE CORRELATIONS     |
# |==============================|

# Temperature - STJ
plot_pca(stj_1sw_temp, arrow_scale = 5, highlight_var = "return_1sw", title_ = "Temperature St. Jean")
plot_cor(data = stj_1sw_temp, title_ = "Temperature St. Jean")
# No variables should be removed a priori strictly on the basis of correlations

plot_pca(stj_2sw_temp, arrow_scale = 5, highlight_var = "return_2sw", title_ = "Temperature St. Jean")
plot_cor(data = stj_2sw_temp, title_ = "Temperature St. Jean")
# No variables should be removed a priori strictly on the basis of correlations

# Temperature - TRI
plot_pca(tri_1sw_temp, arrow_scale = 5, highlight_var = "return_1sw", title_ = "Temperature Trinite")
plot_cor(data = tri_1sw_temp, title_ = "Temperature Trinite")
# No variables should be removed a priori strictly on the basis of correlations

plot_pca(tri_2sw_temp, arrow_scale = 5, highlight_var = "return_2sw", title_ = "Temperature Trinite")
plot_cor(data = tri_2sw_temp, title_ = "Temperature Trinite")
# No variables should be removed a priori strictly on the basis of correlations

# Exploitation rate - STJ
plot_pca(stj_2sw_exp, arrow_scale = 5, highlight_var = "return_2sw", title_ = "Exploitation rate St. Jean")
plot_cor(data = stj_2sw_exp, title_ = "Exploitation rate St. Jean")
# er_WG_lag1 and er2NAC_lag2 are highly correlated. It is likely that one of them will have to be removed.
# It is more likely that er2NAC_lag2 will be choose to be remove since 2sw salmon spend less time in that area and are more likely to be affected by greenland fisheries 

# Exploitation rate - TRI
plot_pca(tri_2sw_exp, arrow_scale = 5, highlight_var = "return_2sw", title_ = "Exploitation rate Trinite")
plot_cor(data = tri_2sw_exp, title_ = "Exploitation rate Trinite")
# Again, er_WG_lag1 and er2NAC_lag2 are highly correlated. It is likely that one of them will have to be removed.
# It is more likely that er2NAC_lag2 will be choose to be remove since 2sw salmon spend less time in that area and are more likely to be affected by greenland fisheries 

# |============================================|
# |     EXPLORE ALL VARIABLES CORRELATIONS     |
# |============================================|

# -------------------------------------------------------- #
#                     ST. JEAN — 1SW                       #
# -------------------------------------------------------- #
plot_pca(stj_1sw, arrow_scale = 5, highlight_var = "return_1sw", title_ = "All variables - St. Jean (1sw)")
plot_cor(data = stj_1sw, title_ = "All variables - St. Jean (1sw)")
# er1NAC_lag1 <-> food_index : r = 0.76 — ces variables devront être à surveiller mais ne seront pas retirées pour l'instant.
# er1NAC_lag1 <-> plankton_index : r = 0.75 — ces variables devront être à surveiller mais ne seront pas retirées pour l'instant.
# LabradorSea_TH <-> food_index : r = -0.70 — ces variables devront être à surveiller mais ne seront pas retirées pour l'instant.
# LabradorSea_TH <-> plankton_index : r = -0.70 — ces variables devront être à surveiller mais ne seront pas retirées pour l'instant.
# À noter que la forte corrélation entre food_index et plankton_index est attendue et normale —
# une seule de ces deux variables sera considérée dans les modèles. Aussi, plankton_t ne sera sans doute pas utilisée. 


# -------------------------------------------------------- #
#                     ST. JEAN — 2SW                       #
# -------------------------------------------------------- #
plot_pca(stj_2sw, arrow_scale = 5, highlight_var = "return_2sw", title_ = "All variables - St. Jean (2sw)")
plot_cor(data = stj_2sw, title_ = "All variables - St. Jean (2sw)")
# er2NAC_lag2 <-> food_index : r = 0.86 — ces variables devront être à surveiller mais ne seront pas retirées pour l'instant.
# er2NAC_lag2 <-> plankton_index : r = 0.85 — ces variables devront être à surveiller mais ne seront pas retirées pour l'instant.
# er2NAC_lag2 <-> er_WG_lag1 : r = 0.87 — ces variables devront être à surveiller mais ne seront pas retirées pour l'instant.
# er_WG_lag1 <-> food_index : r = 0.80 — ces variables devront être à surveiller mais ne seront pas retirées pour l'instant.
# er_WG_lag1 <-> plankton_index : r = 0.76 — ces variables devront être à surveiller mais ne seront pas retirées pour l'instant.
# LabradorSea_TH <-> food_index : r = -0.71 — ces variables devront être à surveiller mais ne seront pas retirées pour l'instant.
# LabradorSea_TH <-> plankton_index : r = -0.69 — ces variables devront être à surveiller mais ne seront pas retirées pour l'instant.
# er2NAC_lag2 <-> pfa_2sw : r = 0.71 — ces variables devront être à surveiller mais ne seront pas retirées pour l'instant.
# er2NAC_lag2 <-> LabradorSea_TH : r = -0.70 — ces variables devront être à surveiller mais ne seront pas retirées pour l'instant.
# À noter que la forte corrélation entre food_index et plankton_index est attendue et normale —
# une seule de ces deux variables sera considérée dans les modèles.

# -------------------------------------------------------- #
#                      TRINITE — 1SW                       #
# -------------------------------------------------------- #
plot_pca(tri_1sw, arrow_scale = 5, highlight_var = "return_1sw", title_ = "All variables - Trinite (1sw)")
plot_cor(data = tri_1sw, title_ = "All variables - Trinite (1sw)")
# er1NAC_lag1 <-> food_index : r = 0.81 — ces variables devront être à surveiller mais ne seront pas retirées pour l'instant.
# pfa_1sw <-> er1NAC_lag1 : r = 0.82 — ces variables devront être à surveiller mais ne seront pas retirées pour l'instant.
# À noter que la forte corrélation entre food_index et plankton_index est attendue et normale —
# une seule de ces deux variables sera considérée dans les modèles.

# -------------------------------------------------------- #
#                      TRINITE — 2SW                       #
# -------------------------------------------------------- #
plot_pca(tri_2sw, arrow_scale = 5, highlight_var = "return_2sw", title_ = "All variables - Trinite (2sw)")
plot_cor(data = tri_2sw, title_ = "All variables - Trinite (2sw)")
# er_WG_lag1 <-> food_index : r = 0.85 — ces variables devront être à surveiller mais ne seront pas retirées pour l'instant.
# er_WG_lag1 <-> plankton_index : r = 0.72 — ces variables devront être à surveiller mais ne seront pas retirées pour l'instant.
# er2NAC_lag2 <-> food_index : r = 0.84 — ces variables devront être à surveiller mais ne seront pas retirées pour l'instant.
# er2NAC_lag2 <-> plankton_index : r = 0.73 — ces variables devront être à surveiller mais ne seront pas retirées pour l'instant.
# er2NAC_lag2 <-> er_WG_lag1 : r = 0.82 — ces variables devront être à surveiller mais ne seront pas retirées pour l'instant.
# er2NAC_lag2 <-> LabradorSea_TH : r = -0.71 — ces variables devront être à surveiller mais ne seront pas retirées pour l'instant.
# pfa_2sw <-> er_WG_lag1 : r = 0.71 — ces variables devront être à surveiller mais ne seront pas retirées pour l'instant.
# pfa_2sw <-> er2NAC_lag2 : r = 0.84 — ces variables devront être à surveiller mais ne seront pas retirées pour l'instant.
# À noter que la forte corrélation entre food_index et plankton_index est attendue et normale —
# une seule de ces deux variables sera considérée dans les modèles.


# |============================================================|
# |     RÉSUMÉ DES VARIABLES À SURVEILLER (PRÉ-SÉLECTION)      |
# |============================================================|

# --- Variables alimentaires (food_index, plankton_index) ---
# La forte corrélation entre ces deux variables est attendue et normale : elles
# mesurent toutes deux la disponibilité alimentaire en mer. Une seule sera retenue
# par modèle. Les VIF seront testés séparément pour chaque version.

# --- er2NAC_lag2 (retrait fortement envisagé dans les modèles 2SW) ---
# Corrélée avec er_WG_lag1, food_index, plankton_index et pfa_2sw dans les deux rivières.
# Argument biologique appuyant le retrait : les saumons 2SW passent moins de
# temps dans la zone NAC et on s'attend à ce qu'ils soient davantage affectés par les pêcheries du
# Groenland (er_WG_lag1 est donc préféré). Ces variables devront être à
# surveiller mais ne seront pas retirées avant l'étape du VIF.

# --- er1NAC_lag1 (à surveiller dans les modèles 1SW) ---
# Corrélée avec food_index et plankton_index dans les deux rivières.
# Corrélée avec pfa_1sw dans TRI.
# Ces variables devront être à surveiller mais ne seront pas retirées pour l'instant.

# --- er_WG_lag1 (à surveiller dans les modèles 2SW) ---
# Corrélée avec food_index, plankton_index et er2NAC_lag2 dans les deux rivières.
# Ces variables devront être à surveiller mais ne seront pas retirées pour l'instant.

# --- LabradorSea_TH (à surveiller dans tous les modèles) ---
# Corrélée avec food_index et plankton_index dans les deux rivières et les deux stratégies.
# Ces variables devront être à surveiller mais ne seront pas retirées pour l'instant.

# --- Aucune variable ne sera retirée sur la seule base des corrélations ---
# Les décisions de retrait seront prises à l'étape du VIF, en tenant compte
# à la fois des corrélations ET de la justification biologique.


# |====================================================|
# |     SELECT VARIABLES BASED ON VIF: TEMPERATURE     |
# |====================================================|

sort(vif(lm(return_1sw ~ cohort + delta_T_sst_minus_river + StJeanCorridor_TH  + BelleIsle_TH + LabradorEntrance_TH + LabradorSea_TH, data = stj_1sw)))
sort(vif(lm(return_2sw ~ cohort + delta_T_sst_minus_river + StJeanCorridor_TH  + BelleIsle_TH + LabradorEntrance_TH + LabradorSea_TH, data = stj_2sw)))
sort(vif(lm(return_1sw ~ cohort + delta_T_sst_minus_river + TriniteCorridor_TH + BelleIsle_TH + LabradorEntrance_TH + LabradorSea_TH, data = tri_1sw)))
sort(vif(lm(return_2sw ~ cohort + delta_T_sst_minus_river + TriniteCorridor_TH + BelleIsle_TH + LabradorEntrance_TH + LabradorSea_TH, data = tri_2sw)))
# Aucune variable n'a un VIF > 10. Aucune n'est retirée.


# |=====================================================================|
# |     SELECT VARIABLES BASED ON VIF: EXPLOITATION RATE (2SW ONLY)     |
# |=====================================================================|

sort(vif(lm(return_2sw ~ cohort + er_WG_lag1 + er2NAC_lag2, data = stj_2sw)))
sort(vif(lm(return_2sw ~ cohort + er_WG_lag1 + er2NAC_lag2, data = tri_2sw)))
# Aucune variable n'a un VIF > 10. Aucune n'est retirée.


# |======================================================|
# |     SELECT VARIABLES BASED ON VIF: ALL VARIABLES     |
# |======================================================|

# -------------------------------------------------------- #
#                     ST. JEAN — 1SW                       #
# -------------------------------------------------------- #

sort(vif(lm(return_1sw ~ cohort + delta_T_sst_minus_river + plankton_index +
              StJeanCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
              LabradorSea_TH + smolt_k + er1NAC_lag1 + pfa_1sw, data = stj_1sw)))
# Aucune variable n'a un VIF > 10. Aucune n'est retirée.


# -------------------------------------------------------- #
#                      TRINITE — 1SW                       #
# -------------------------------------------------------- #

sort(vif(lm(return_1sw ~ cohort + delta_T_sst_minus_river + plankton_index +
              TriniteCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
              LabradorSea_TH + smolt_k + er1NAC_lag1 + pfa_1sw, data = tri_1sw)))
# er1NAC_lag1 a un VIF de 11.42
# pfa_1sw a un VIF de 6.52
# Bien que la règle habituelle consiste à retirer la variable affichant le VIF le plus
# élevé, nous avons plutôt exclu pfa_1sw pour des raisons biologiques et méthodologiques.
# Premièrement, er1NAC_lag1 est un prédicteur à fort ancrage théorique : la mortalité 
# par la pêche en mer représente un forçage externe direct sur les cohortes constitue
# une hypothèse centrale de notre cadre d'analyse. Deuxièmement, la colinéarité entre
# ces deux variables (qui décrivent toutes deux la pression de mortalité en mer) suggère
# qu'elles encodent une information largement redondante; conserver er1NAC_lag1 permet
# ainsi de représenter ce processus avec la variable dont l'interprétation causale est
# la plus directe. Troisièmement, l'exclusion de pfa_1sw a résolu efficacement le 
# problème de colinéarité : le VIF de er1NAC_lag1 est tombé à 3.41 et aucun prédicteur
# ne dépassait 3.5 dans le modèle résultant, contre un VIF maximal de 2.32 pour
# plankton_index lorsque c'est er1NAC_lag1 qui est exclu. Les deux solutions étant 
# donc acceptables sur le plan statistique, mais la première étant biologiquement préférable.


sort(vif(lm(return_1sw ~ cohort + delta_T_sst_minus_river + plankton_index +
              TriniteCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
              LabradorSea_TH + smolt_k + er1NAC_lag1, data = tri_1sw))) # Sans pfa_1sw
sort(vif(lm(return_1sw ~ cohort + delta_T_sst_minus_river + plankton_index +
              TriniteCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
              LabradorSea_TH + smolt_k + pfa_1sw, data = tri_1sw))) # Sans er1NAC_lag1
# Le modèle retenu est celui sans pfa_1sw. 


# -------------------------------------------------------- #
#                     ST. JEAN — 2SW                       #
# -------------------------------------------------------- #

# Test AVEC les variables er_WG_lag1 et er2NAC_lag2 non combinées en un seul indice (er_index)
sort(vif(lm(return_2sw ~ cohort + delta_T_sst_minus_river + plankton_index +
              StJeanCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
              LabradorSea_TH + smolt_k + er_WG_lag1 + er2NAC_lag2 + pfa_2sw, data = stj_2sw)))
# pfa_2sw a le VIF le plus élevé (23.61), suivit par cohort (21.08) et er2NAC_lag2 (15.61). 

# Test SANS les variables er_WG_lag1 et er2NAC_lag2, remplacées par er_index
sort(vif(lm(return_2sw ~ cohort + delta_T_sst_minus_river + plankton_index +
              StJeanCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
              LabradorSea_TH + smolt_k + er_index + pfa_2sw, data = stj_2sw)))
# pfa_2sw demeure le VIF le plus élevé (21.56), suivit par cohort (21.14). 

# Test sans pfa_2sw (modèle avec er_WG_lag1 + er2NAC_lag2 au lieu de er_index)
sort(vif(lm(return_2sw ~ cohort + delta_T_sst_minus_river + plankton_index +
              StJeanCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
              LabradorSea_TH + smolt_k + er_WG_lag1 + er2NAC_lag2, data = stj_2sw)))
# er2NAC_lag2 demeure élevé (14.17)

# Test sans pfa_2sw (modèle avec er_index au lieu de er_WG_lag1 + er2NAC_lag2)
sort(vif(lm(return_2sw ~ cohort + delta_T_sst_minus_river + plankton_index +
              StJeanCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
              LabradorSea_TH + smolt_k + er_index, data = stj_2sw)))
# Il semble que pfa_2sw doive être retiré. 

# Test avec pfa_2sw, mais sans aucune variable d'exploitation (dans ni er_index ni er_WG_lag1 + er2NAC_lag2)
sort(vif(lm(return_2sw ~ cohort + delta_T_sst_minus_river + plankton_index +
              StJeanCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
              LabradorSea_TH + smolt_k + pfa_2sw, data = stj_2sw)))
# Même en retirant complètement la pêche, pfa_2sw demeure élevé (VIF = 20.73). Il semble de plus en plus pertinent de la retirer malgré son importance biologique.

# Pour finir, test avec pfa, mais sans les variables qui concernent uniquement la première année en mer (qui ont été testées chez les 1SW; et notons que 1SW est corrélé avec 2SW)
sort(vif(lm(return_2sw ~ cohort + plankton_index +
              LabradorSea_TH + er_index + pfa_2sw, data = stj_2sw)))
# pfa_2sw (10.11) et cohort (10.57) restent légèrement > 10 

# Conclusion, er_index est conservé et pfa_2sw est retiré.


# -------------------------------------------------------- #
#                      TRINITE — 2SW                       #
# -------------------------------------------------------- #

# Test AVEC les variables er_WG_lag1 et er2NAC_lag2 non combinées en un seul indice (er_index)
sort(vif(lm(return_2sw ~ cohort + delta_T_sst_minus_river + plankton_index +
              TriniteCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
              LabradorSea_TH + smolt_k + er_WG_lag1 + er2NAC_lag2 + pfa_2sw, data = tri_2sw)))
# pfa_2sw a le VIF le plus élevé (21.67), suivit par er2NAC_lag2 (13.95).

# Test SANS les variables er_WG_lag1 et er2NAC_lag2, remplacées par er_index
sort(vif(lm(return_2sw ~ cohort + delta_T_sst_minus_river + plankton_index +
              TriniteCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
              LabradorSea_TH + smolt_k + er_index + pfa_2sw, data = tri_2sw)))
# pfa_2sw demeure le VIF le plus élevé (13.62)

# Test SANS pfa_2sw et on garde er_index
sort(vif(lm(return_2sw ~ cohort + delta_T_sst_minus_river + plankton_index +
              TriniteCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
              LabradorSea_TH + smolt_k + er_index, data = tri_2sw)))
# c'est le modèle qui sera utilisé


# |============================================================|
# |     RÉSUMÉ DES DÉCISIONS POST-VIF                          |
# |============================================================|

# STJ — 1SW
#   Aucune variable retirée. er1NAC_lag1 et pfa_1sw ont des VIF acceptables
#   ensemble (max VIF = 4.66). Un seul modèle est testé.

# TRI — 1SW
#   er1NAC_lag1 (VIF = 11.42) et pfa_1sw (VIF = 6.52) sont mutuellement
#   collinéaires. pfa_1sw est retirée pour des raisons biologiques et
#   méthodologiques : er1NAC_lag1 représente un forçage externe direct sur les
#   cohortes et constitue une hypothèse centrale du cadre d'analyse. Son retrait
#   ramène tous les VIF sous 3.5. Un seul modèle est testé.

# STJ — 2SW
#   pfa_2sw est structurellement collinéaire avec cohort (VIF > 10 même sans
#   aucune variable d'exploitation), ce qui indique une tendance temporelle
#   commune irrésoluble. pfa_2sw est retirée. er_WG_lag1 et er2NAC_lag2 sont
#   remplacés par er_index (PCA sur échelle log), ce qui résout la collinéarité
#   entre les deux taux d'exploitation et ramène tous les VIF sous 10.
#   Un seul modèle est testé.

# TRI — 2SW
#   pfa_2sw est structurellement collinéaire avec cohort (VIF = 10.25 même
#   après remplacement de er_WG_lag1 + er2NAC_lag2 par er_index), ce qui
#   indique une tendance temporelle commune irrésoluble. pfa_2sw est retirée.
#   er_WG_lag1 et er2NAC_lag2 sont remplacés par er_index (PCA sur échelle
#   log), ce qui résout la collinéarité entre les deux taux d'exploitation et
#   ramène tous les VIF sous 7. Un seul modèle est testé.


# |====================================|
# |     DEFINE FORMULAS FOR MODELS     |
# |====================================|

# STJ — 1SW
f_stj_1sw <- as.formula(
  "return_1sw ~ cohort + delta_T_sst_minus_river + plankton_index +
   StJeanCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
   LabradorSea_TH + smolt_k + er1NAC_lag1 + pfa_1sw"
)

# TRI — 1SW
f_tri_1sw <- as.formula(
  "return_1sw ~ cohort + delta_T_sst_minus_river + plankton_index +
   TriniteCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
   LabradorSea_TH + smolt_k + er1NAC_lag1"
)

# STJ — 2SW
f_stj_2sw <- as.formula(
  "return_2sw ~ cohort + delta_T_sst_minus_river + plankton_index +
   StJeanCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
   LabradorSea_TH + smolt_k + er_index"
)

# TRI — 2SW
f_tri_2sw <- as.formula(
  "return_2sw ~ cohort + delta_T_sst_minus_river + plankton_index +
   TriniteCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
   LabradorSea_TH + smolt_k + er_index"
)


# |==============================|
# |     EXPLORE DISTRIBUTION     |
# |==============================|

par(mfrow = c(1,1))

# -------------------------------------------------------- #
#                     ST. JEAN — 1SW                       #
# -------------------------------------------------------- #

summary(stj_1sw$return_1sw)

ggplot(stj_1sw, aes(x = return_1sw)) +
  geom_histogram(aes(y = after_stat(density)), bins = 10, fill = "grey80", color = "grey30") +
  geom_density(color = "black", linewidth = 1) +
  labs(title = "St. Jean - Distribution of return rate (1SW)", x = "Return rate", y = "Density") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

stj_1sw_beta <- gam(f_stj_1sw, family = betar(link = "logit"), data = stj_1sw)

stj_1sw_sim <- simulateResiduals(stj_1sw_beta)

plot(stj_1sw_sim)
testDispersion(stj_1sw_sim)
# Aucune déviation


# -------------------------------------------------------- #
#                      TRINITE — 1SW                       #
# -------------------------------------------------------- #

summary(tri_1sw$return_1sw)

ggplot(tri_1sw, aes(x = return_1sw)) +
  geom_histogram(aes(y = after_stat(density)), bins = 10, fill = "grey80", color = "grey30") +
  geom_density(color = "black", linewidth = 1) +
  labs(title = "Trinite - Distribution of return rate (1SW)", x = "Return rate", y = "Density") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

tri_1sw_beta <- gam(f_tri_1sw, family = betar(link = "logit"), data = tri_1sw)

tri_1sw_sim <- simulateResiduals(tri_1sw_beta)

plot(tri_1sw_sim)
testDispersion(tri_1sw_sim)
# Aucune déviation


# -------------------------------------------------------- #
#                     ST. JEAN — 2SW                       #
# -------------------------------------------------------- #

summary(stj_2sw$return_2sw)

ggplot(stj_2sw, aes(x = return_2sw)) +
  geom_histogram(aes(y = after_stat(density)), bins = 10, fill = "grey80", color = "grey30") +
  geom_density(color = "black", linewidth = 1) +
  labs(title = "St. Jean - Distribution of return rate (2SW)", x = "Return rate", y = "Density") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

stj_2sw_beta <- gam(f_stj_2sw, family = betar(link = "logit"), data = stj_2sw)

stj_2sw_sim <- simulateResiduals(stj_2sw_beta)

plot(stj_2sw_sim)
testDispersion(stj_2sw_sim)
# Aucune déviation


# -------------------------------------------------------- #
#                      TRINITE — 2SW                       #
# -------------------------------------------------------- #

summary(tri_2sw$return_2sw)

ggplot(tri_2sw, aes(x = return_2sw)) +
  geom_histogram(aes(y = after_stat(density)), bins = 10, fill = "grey80", color = "grey30") +
  geom_density(color = "black", linewidth = 1) +
  labs(title = "Trinite - Distribution of return rate (2SW)", x = "Return rate", y = "Density") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

tri_2sw_beta <- gam(f_tri_2sw, family = betar(link = "logit"), data = tri_2sw)

tri_2sw_sim <- simulateResiduals(tri_2sw_beta)

plot(tri_2sw_sim)
testDispersion(tri_2sw_sim)
# Aucune déviation


# |=======================================|
# |     TESTING FOR DATA INDEPENDENCE     |
# |=======================================|


# -------------------------------------------------------- #
#                     ST. JEAN — 1SW                       #
# -------------------------------------------------------- #

stj_1sw_ind <- betareg(f_stj_1sw, data = stj_1sw, link = "logit")

par(mfrow = c(1, 1))
acf(residuals(stj_1sw_ind))
lmtest::dwtest(stj_1sw_ind)
# Pas de corrélation temporelle importante


# -------------------------------------------------------- #
#                      TRINITE — 1SW                       #
# -------------------------------------------------------- #

tri_1sw_ind <- betareg(f_tri_1sw, data = tri_1sw, link = "logit")

par(mfrow = c(1, 1))
acf(residuals(tri_1sw_ind))
lmtest::dwtest(tri_1sw_ind)
# Pas de corrélation temporelle importante


# -------------------------------------------------------- #
#                     ST. JEAN — 2SW                       #
# -------------------------------------------------------- #

stj_2sw_ind <- betareg(f_stj_2sw, data = stj_2sw, link = "logit")

par(mfrow = c(1, 1))
acf(residuals(stj_2sw_ind))
lmtest::dwtest(stj_2sw_ind)
# Pas de corrélation temporelle importante


# -------------------------------------------------------- #
#                      TRINITE — 2SW                       #
# -------------------------------------------------------- #

tri_2sw_ind <- betareg(f_tri_2sw, data = tri_2sw, link = "logit")

par(mfrow = c(1, 1))
acf(residuals(tri_2sw_ind))
lmtest::dwtest(tri_2sw_ind)
# Légère corrélation temporelle, considérée comme acceptable


# |=========================================================|
# |     CHECK FOR OVERDISPERSION AND HETEROSCEDASTICITY     |
# |=========================================================|

# -------------------------------------------------------- #
#                     ST. JEAN — 1SW                       #
# -------------------------------------------------------- #

mean(stj_1sw$return_1sw)
var(stj_1sw$return_1sw)

plot(resid(stj_1sw_beta)); abline(h = 0, col = "red", lty = 2)
sum(resid(stj_1sw_beta, type = "pearson")^2) / stj_1sw_beta$df.resid


# -------------------------------------------------------- #
#                      TRINITE — 1SW                       #
# -------------------------------------------------------- #

mean(tri_1sw$return_1sw)
var(tri_1sw$return_1sw)

plot(resid(tri_1sw_beta)); abline(h = 0, col = "red", lty = 2)
sum(resid(tri_1sw_beta, type = "pearson")^2) / tri_1sw_beta$df.resid


# -------------------------------------------------------- #
#                     ST. JEAN — 2SW                       #
# -------------------------------------------------------- #

mean(stj_2sw$return_2sw)
var(stj_2sw$return_2sw)

plot(resid(stj_2sw_beta)); abline(h = 0, col = "red", lty = 2)
sum(resid(stj_2sw_beta, type = "pearson")^2) / stj_2sw_beta$df.resid


# -------------------------------------------------------- #
#                      TRINITE — 2SW                       #
# -------------------------------------------------------- #

mean(tri_2sw$return_2sw)
var(tri_2sw$return_2sw)

plot(resid(tri_2sw_beta)); abline(h = 0, col = "red", lty = 2)
sum(resid(tri_2sw_beta, type = "pearson")^2) / tri_2sw_beta$df.resid

# |================================|
# |     EXPORT FINAL DATA SETS     |
# |================================|

write.fst(stj_1sw, file.path(base_path, "04_Analysis/Sea/Final_DataSet_STJ_1SW.fst"))
write.fst(stj_2sw, file.path(base_path, "04_Analysis/Sea/Final_DataSet_STJ_2SW.fst"))
write.fst(tri_1sw, file.path(base_path, "04_Analysis/Sea/Final_DataSet_TRI_1SW.fst"))
write.fst(tri_2sw, file.path(base_path, "04_Analysis/Sea/Final_DataSet_TRI_2SW.fst"))


# |================================|
# |     FINAL MODELS               |
# |================================|

# -------------------------------------------------------- #
#                     ST. JEAN — 1SW                       #
# -------------------------------------------------------- #
# Variables retirées : aucune
# Un seul modèle testé.
# Formule finale : return_1sw ~ cohort + delta_T_sst_minus_river + plankton_index +
#                               StJeanCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
#                               LabradorSea_TH + smolt_k + er1NAC_lag1 + pfa_1sw


# -------------------------------------------------------- #
#                      TRINITE — 1SW                       #
# -------------------------------------------------------- #
# pfa_1sw est retirée et er1NAC_lag1 est conservée
#
# Formule finale : return_1sw ~ cohort + delta_T_sst_minus_river + plankton_index +
#                               TriniteCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
#                               LabradorSea_TH + smolt_k + er1NAC_lag1


# -------------------------------------------------------- #
#                     ST. JEAN — 2SW                       #
# -------------------------------------------------------- #
# retrait de pfa_2sw et utilisation de er_index 
#
# Formule finale : return_2sw ~ cohort + delta_T_sst_minus_river + plankton_index +
#                               StJeanCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
#                               LabradorSea_TH + smolt_k + er_index


# -------------------------------------------------------- #
#                      TRINITE — 2SW                       #
# -------------------------------------------------------- #
# er_index remplace er_WG_lag1 + er2NAC_lag2. Après ce remplacement,
# pfa_2sw passe sous le seuil (VIF = 9.36 < 10) et les deux variables
# coexistent dans un seul modèle.
#
# Un seul modèle testé.
# Formule finale : return_2sw ~ cohort + delta_T_sst_minus_river + plankton_index +
#                               TriniteCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
#                               LabradorSea_TH + smolt_k + er_index 





###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################

















