

# ══════════════════════════════════════════════════════════════════════════════
# HEADER
# ══════════════════════════════════════════════════════════════════════════════
#
# Assembles the four pre-selection datasets (St. Jean and Trinite, 1SW and 2SW)
# used in the marine survival models. For each river x sea-age combination,
# the response variable (return rate) is joined by cohort to the marine
# covariates: SST-river temperature difference at sea entry, food availability
# (Labrador Sea plankton, year t), annual thermal habitat, smolt Fulton's K,
# pre-fishery abundance (density), and exploitation rates. For 1SW, only
# er1NAC_lag1 is kept. For 2SW, er_WG_lag1 and er2NAC_lag2 are combined into a
# single exploitation index (er_index = PC1 of a PCA on their log values,
# oriented so that higher values mean higher fishing pressure).
#
# Inputs (R/03_clean_data/Sea/):
#   - DeltaT_SST_River_{stj,tri}_{1SW,2SW}.fst
#       SST minus river temperature difference at sea entry
#   - Plankton_LabradorSea_{stj,tri}_{1SW,2SW}.fst
#       Food availability index (Plankton_t)
#   - TH_annual_{stj,tri}_{1sw,2sw}.fst
#       Annual thermal habitat indices by polygon
#   - ExploitationRate_Sea.fst
#       Cohort-indexed exploitation rates (er1NAC_lag1, er_WG_lag1, er2NAC_lag2)
#   - Density_Lab.fst
#       Pre-fishery abundance (pfa_1sw, pfa_2sw)
#
# Inputs (R/00_raw_data/salmon/):
#   - Salmon_Monitoring_StJean_Trinite_2024.csv
#       Response variable (return_1sw_pct, return_2sw_pct) and smolt
#       Fulton's K (smolt_k); year = sea-entry year
#
# Outputs (R/04_analysis/Sea/Final_Data/):
#   - PreSelection_DataSet_1sw_stj.fst
#   - PreSelection_DataSet_2sw_stj.fst
#   - PreSelection_DataSet_1sw_tri.fst
#   - PreSelection_DataSet_2sw_tri.fst
#       One row per cohort; response + candidate covariates before VIF screening
# ══════════════════════════════════════════════════════════════════════════════


# ══════════════════════════════════════════════════════════════════════════════
# WORKING DIRECTORY
# ══════════════════════════════════════════════════════════════════════════════

# Sets the working directory to the repository root (the folder containing the
# .Rproj file), whatever the folder is named. 
find_repo_root <- function() {
  start <- tryCatch(dirname(rstudioapi::getSourceEditorContext()$path),
                    error = function(e) "")
  if (!nzchar(start)) {
    f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    start <- if (length(f)) dirname(normalizePath(f)) else getwd()
  }
  d <- normalizePath(start, winslash = "/")
  while (length(list.files(d, pattern = "\\.Rproj$")) == 0) {
    if (dirname(d) == d) stop("Repository root not found (no .Rproj file above this script).")
    d <- dirname(d)
  }
  d
}
setwd(find_repo_root())

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
# DeltaT_SST_River
#    Keep: river, cohort (= year), delta_T_sst_minus_river
#    Remove NAs in delta_T_sst_minus_river
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
# Food availability
#    Keep: all columns except sw_type
#    Remove NAs in plankton_index
# ─────────────────────────────────────────────

food_stj_1sw <- read_fst(file.path(base_path, "R/03_clean_data/Sea/Plankton_LabradorSea_stj_1SW.fst")) %>%
  dplyr::select(-sw_type) %>%
  filter(!is.na(Plankton_t))

food_stj_2sw <- read_fst(file.path(base_path, "R/03_clean_data/Sea/Plankton_LabradorSea_stj_2SW.fst")) %>%
  dplyr::select(-sw_type, -Plankton_tp1) %>%
  filter(!is.na(Plankton_t))

food_tri_1sw <- read_fst(file.path(base_path, "R/03_clean_data/Sea/Plankton_LabradorSea_tri_1SW.fst")) %>%
  dplyr::select(-sw_type) %>%
  filter(!is.na(Plankton_t))

food_tri_2sw <- read_fst(file.path(base_path, "R/03_clean_data/Sea/Plankton_LabradorSea_tri_2SW.fst")) %>%
  dplyr::select(-sw_type, -Plankton_tp1) %>%
  filter(!is.na(Plankton_t))


# ─────────────────────────────────────────────
# TH_annual
#    Keep: all columns + add river
#    Remove NAs in all columns except cohort
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
# Smolt Fulton's K
#    SOURCE: MELCCFP 2024 annual report (Table 5)
#    year = sea-entry year of the smolts
#    Remove NAs in smolt_k
# ─────────────────────────────────────────────

smolt_k <- read.csv(file.path(base_path,
                              "R/00_raw_data/salmon/Salmon_Monitoring_StJean_Trinite_2024.csv")) %>%
  dplyr::select(river, cohort = year, smolt_k) %>%
  filter(!is.na(smolt_k))


# ─────────────────────────────────────────────
# ExploitationRate_Sea
#    Keep: all columns
#    Remove NAs in er1NAC_lag1, er_WG_lag1 or er2NAC_lag2
#    (same dataset for stj and tri, joined by cohort)
# ─────────────────────────────────────────────

exploit <- read_fst(file.path(base_path, "R/03_Clean_Data/Sea/ExploitationRate_Sea.fst")) %>%
  filter(!is.na(er1NAC_lag1) & !is.na(er_WG_lag1) & !is.na(er2NAC_lag2))


# ─────────────────────────────────────────────
# Density_Lab
#    Keep: all columns
#    (same dataset for stj and tri, joined by cohort)
# ─────────────────────────────────────────────

density_lab <- read_fst(file.path(base_path, "R/03_Clean_Data/Sea/Density_Lab.fst"))


# ─────────────────────────────────────────────
# Return rates (response variable — no rows removed)
#    SOURCE: MELCCFP 2024 annual report (Tables 17-18)
#    year = sea-entry year of the smolts
# ─────────────────────────────────────────────

survival <- read.csv(file.path(base_path,
                               "R/00_raw_data/salmon/Salmon_Monitoring_StJean_Trinite_2024.csv")) %>%
  dplyr::select(river, cohort = year,
                return_1sw = return_1sw_pct,
                return_2sw = return_2sw_pct)


# ═════════════════════════════════════════════
# FINAL COMBINATION: 4 datasets
# full_join on cohort (and river when relevant)
# to keep all available years
#
# Combined datasets: survival (base), deltaT, food (Plankton_t), th,
#                    smolt_k, exploit, density.
# Removed from the pipeline: capelin, smolt_tl, smolt_weight.
# ═════════════════════════════════════════════

combine_data <- function(river_id, sw_label,
                         deltaT, food, th,
                         survival_data, exploit_data,
                         smolt_k_data, density_data) {
  
  # Survival (response variable): start here to keep all years
  base <- survival_data %>% filter(river == river_id)
  
  # Variables shared by both rivers -> filtered by river
  ex <- exploit_data                              # no river column in the source file
  sk <- smolt_k_data %>% filter(river == river_id)
  
  base %>%
    full_join(deltaT,       by = c("cohort", "river")) %>%
    full_join(food,         by = c("cohort", "river")) %>%
    full_join(th,           by = c("cohort", "river")) %>%
    full_join(sk,           by = c("cohort", "river")) %>%
    full_join(ex,           by = "cohort") %>%
    full_join(density_data, by = "cohort") %>%
    arrange(cohort)
}

data_stj_1sw <- combine_data(
  river_id = "stj", sw_label = "1sw",
  deltaT   = deltaT_stj_1sw,
  food     = food_stj_1sw,
  th       = th_stj_1sw,
  survival_data = survival %>% dplyr::select(cohort, river, return_1sw),
  exploit_data  = exploit,
  smolt_k_data  = smolt_k,
  density_data  = density_lab %>% dplyr::select(YEAR, pfa_1sw) %>% rename(cohort = YEAR)
) %>%
  dplyr::select(-er_WG_lag1, -er2NAC_lag2)

data_stj_2sw <- combine_data(
  river_id = "stj", sw_label = "2sw",
  deltaT   = deltaT_stj_2sw,
  food     = food_stj_2sw,
  th       = th_stj_2sw,
  survival_data = survival %>% dplyr::select(cohort, river, return_2sw),
  exploit_data  = exploit,
  smolt_k_data  = smolt_k,
  density_data  = density_lab %>% dplyr::select(YEAR, pfa_2sw) %>% dplyr::rename(cohort = YEAR)
)

data_tri_1sw <- combine_data(
  river_id = "tri", sw_label = "1sw",
  deltaT   = deltaT_tri_1sw,
  food     = food_tri_1sw,
  th       = th_tri_1sw,
  survival_data = survival %>% dplyr::select(cohort, river, return_1sw),
  exploit_data  = exploit,
  smolt_k_data  = smolt_k,
  density_data  = density_lab %>% dplyr::select(YEAR, pfa_1sw) %>% rename(cohort = YEAR)
) %>%
  dplyr::select(-er_WG_lag1, -er2NAC_lag2)

data_tri_2sw <- combine_data(
  river_id = "tri", sw_label = "2sw",
  deltaT   = deltaT_tri_2sw,
  food     = food_tri_2sw,
  th       = th_tri_2sw,
  survival_data = survival %>% dplyr::select(cohort, river, return_2sw),
  exploit_data  = exploit,
  smolt_k_data  = smolt_k,
  density_data  = density_lab %>% dplyr::select(YEAR, pfa_2sw) %>% rename(cohort = YEAR)
)

# ══════════════════════════════════════════════════════════════════════════════
# CREATION OF A SYNTHETIC INDEX FOR EXPLOITATION RATES (2SW)
# ══════════════════════════════════════════════════════════════════════════════
#
# er_WG_lag1 and er2NAC_lag2 measure the same phenomenon (fishing pressure)
# in two different but biologically linked geographic areas for 2SW salmon.
# A synthetic index derived from the first axis of a PCA is built to resolve
# their mutual collinearity while retaining the biological information of
# both variables.
#
# Steps:
#   1. Correlation between the two variables (raw and log)
#   2. PCA on log scale
#   3. Check of the biological interpretability of the score
#   4. Addition of er_index to the final datasets
#
# ══════════════════════════════════════════════════════════════════════════════


# ── 1. Correlation between er_WG_lag1 and er2NAC_lag2 ────────────────────────
#
# Necessary condition: both variables must be strongly and positively
# correlated for their shared variance to be biologically interpretable on a
# single axis.
#
# Interpretation rule:
#   r ≥ 0.70  → sufficient correlation, PCA index justified
#   r < 0.70  → the two variables capture dynamics that are too distinct;
#               a common index would erase important information
#
# Exploitation rates are also tested on log scale, because the raw
# distributions are strongly skewed. The log version is the one used for
# the PCA (see step 2).

# STJ
cor_stj_raw <- cor.test(data_stj_2sw$er_WG_lag1, data_stj_2sw$er2NAC_lag2, use = "complete.obs")
cor_stj_log <- cor.test(log(data_stj_2sw$er_WG_lag1), log(data_stj_2sw$er2NAC_lag2), use = "complete.obs")

# TRI
cor_tri_raw <- cor.test(data_tri_2sw$er_WG_lag1, data_tri_2sw$er2NAC_lag2, use = "complete.obs")
cor_tri_log <- cor.test(log(data_tri_2sw$er_WG_lag1), log(data_tri_2sw$er2NAC_lag2), use = "complete.obs")

# Results — raw and log correlations
sprintf("STJ — r raw = %.3f  (p = %.4f)", cor_stj_raw$estimate, cor_stj_raw$p.value)
sprintf("STJ — r log = %.3f  (p = %.4f)", cor_stj_log$estimate, cor_stj_log$p.value)
sprintf("TRI — r raw = %.3f  (p = %.4f)", cor_tri_raw$estimate, cor_tri_raw$p.value)
sprintf("TRI — r log = %.3f  (p = %.4f)", cor_tri_log$estimate, cor_tri_log$p.value)

# Side-by-side figures (raw vs log) by river
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
  color = "steelblue", title = "STJ 2SW — Raw scale",
  xlab = "er_WG_lag1", ylab = "er2NAC_lag2"
)
p_stj_log <- make_cor_plot(
  x = log(data_stj_2sw$er_WG_lag1), y = log(data_stj_2sw$er2NAC_lag2),
  r = cor_stj_log$estimate, p = cor_stj_log$p.value,
  color = "steelblue", title = "STJ 2SW — Log scale",
  xlab = "log(er_WG_lag1)", ylab = "log(er2NAC_lag2)"
)
p_tri_raw <- make_cor_plot(
  x = data_tri_2sw$er_WG_lag1, y = data_tri_2sw$er2NAC_lag2,
  r = cor_tri_raw$estimate, p = cor_tri_raw$p.value,
  color = "tomato", title = "TRI 2SW — Raw scale",
  xlab = "er_WG_lag1", ylab = "er2NAC_lag2"
)
p_tri_log <- make_cor_plot(
  x = log(data_tri_2sw$er_WG_lag1), y = log(data_tri_2sw$er2NAC_lag2),
  r = cor_tri_log$estimate, p = cor_tri_log$p.value,
  color = "tomato", title = "TRI 2SW — Log scale",
  xlab = "log(er_WG_lag1)", ylab = "log(er2NAC_lag2)"
)

ggarrange(p_stj_raw, p_stj_log, ncol = 2, nrow = 1)
ggarrange(p_tri_raw, p_tri_log, ncol = 2, nrow = 1)


# ── 2. PCA on er_WG_lag1 and er2NAC_lag2 (log scale) ─────────────────────────
#
# Exploitation rates are log-transformed before the PCA for two reasons.
# First, the effects of fishing mortality are multiplicative: an increase in
# the rate from 2% to 4% is not equivalent to an increase from 20% to 22%,
# even though the absolute difference is identical. The log scale therefore
# better reflects the biological nature of the phenomenon. Second, the raw
# distributions are strongly skewed, with a mass of values concentrated near
# zero and a few extreme observations. Without transformation, axis 1 would be
# dominated by these extreme values. The log transformation homogenizes the
# distribution and ensures that each observation contributes comparably to
# the covariance structure.
#
# The log transformation substantially improves the distribution, but does
# not make it perfectly symmetric. A few points remain away from the main
# linear relationship, which introduces residual uncertainty into the PCA
# score. This index must therefore be interpreted as an approximate gradient
# of overall fishing pressure, not as a precise measure.
#
# What we check:
#   a) % of variance explained by Dim 1
#      → Target: ≥ 70–80 %. Below 60 %, the score is too imprecise.
#   b) Loadings of both variables on Dim 1
#      → Both must point in the same direction (same sign) and have a
#        similar magnitude.
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
# Variance explained by Dim 1 and loadings — STJ
sprintf("STJ — Variance explained by Dim 1: %.1f %%", pca_stj$eig[1, 2])
round(pca_stj$var$coord[, 1, drop = FALSE], 3)
# Variance explained by Dim 1 and loadings — TRI
sprintf("TRI — Variance explained by Dim 1: %.1f %%", pca_tri$eig[1, 2])
round(pca_tri$var$coord[, 1, drop = FALSE], 3)
# Biplots
fviz_pca_var(pca_stj, col.var = "steelblue", repel = TRUE,
             title = "STJ 2SW — PCA log(er_WG_lag1) / log(er2NAC_lag2)")
fviz_pca_var(pca_tri, col.var = "tomato", repel = TRUE,
             title = "TRI 2SW — PCA log(er_WG_lag1) / log(er2NAC_lag2)")
# ── 3. Biological interpretability of the score ──────────────────────────────
#
# The PCA score (individual coordinate on Dim 1) is a centred and scaled
# variable: its mean is 0 and its standard deviation is 1. A high score must
# correspond to high fishing pressure for the index to be biologically
# interpretable in the model.
#
# We check this by computing the correlation of the score with each of the
# two original variables. If both correlations are positive and strong, the
# score is oriented in the right direction. If they are negative, the sign
# of the score is automatically reversed.

data_stj_2sw$er_index <- NA_real_
data_tri_2sw$er_index <- NA_real_
data_stj_2sw$er_index[as.integer(rownames(er_vars_log_stj))] <- pca_stj$ind$coord[, 1]
data_tri_2sw$er_index[as.integer(rownames(er_vars_log_tri))] <- pca_tri$ind$coord[, 1]
cor_stj_wg  <- cor(data_stj_2sw$er_index, data_stj_2sw$er_WG_lag1,  use = "complete.obs")
cor_stj_nac <- cor(data_stj_2sw$er_index, data_stj_2sw$er2NAC_lag2, use = "complete.obs")
cor_tri_wg  <- cor(data_tri_2sw$er_index, data_tri_2sw$er_WG_lag1,  use = "complete.obs")
cor_tri_nac <- cor(data_tri_2sw$er_index, data_tri_2sw$er2NAC_lag2, use = "complete.obs")

# Automatic orientation correction if needed
if (cor_stj_wg < 0) {
  data_stj_2sw$er_index <- data_stj_2sw$er_index * -1
  cor_stj_wg            <- cor_stj_wg  * -1
  cor_stj_nac           <- cor_stj_nac * -1
  # STJ — Score automatically reversed (er_index * -1)
} else {
  # STJ — Score correctly oriented, no change
}

if (cor_tri_wg < 0) {
  data_tri_2sw$er_index <- data_tri_2sw$er_index * -1
  cor_tri_wg            <- cor_tri_wg  * -1
  cor_tri_nac           <- cor_tri_nac * -1
  # TRI — Score automatically reversed (er_index * -1)
} else {
  # TRI — Score correctly oriented, no change
}

# Final correlations of the score with the original variables — STJ
sprintf("STJ — r(er_index, er_WG_lag1)  = %.3f", cor_stj_wg)
sprintf("STJ — r(er_index, er2NAC_lag2) = %.3f", cor_stj_nac)

# Final correlations of the score with the original variables — TRI
sprintf("TRI — r(er_index, er_WG_lag1)  = %.3f", cor_tri_wg)
sprintf("TRI — r(er_index, er2NAC_lag2) = %.3f", cor_tri_nac)




data_stj_2sw <- data_stj_2sw %>%
  dplyr::select(-er_WG_lag1, -er2NAC_lag2, -er1NAC_lag1)
data_tri_2sw <- data_tri_2sw %>%
  dplyr::select(-er_WG_lag1, -er2NAC_lag2, -er1NAC_lag1)



# |=========================|
# |     EXPORT DATASETS     |
# |=========================|

write_fst(data_stj_1sw, path = file.path(base_path, "R/04_analysis/Sea/Final_Data", "PreSelection_DataSet_1sw_stj.fst"), compress = 100)
write_fst(data_stj_2sw, path = file.path(base_path, "R/04_analysis/Sea/Final_Data", "PreSelection_DataSet_2sw_stj.fst"), compress = 100)
write_fst(data_tri_1sw, path = file.path(base_path, "R/04_analysis/Sea/Final_Data", "PreSelection_DataSet_1sw_tri.fst"), compress = 100)
write_fst(data_tri_2sw, path = file.path(base_path, "R/04_analysis/Sea/Final_Data", "PreSelection_DataSet_2sw_tri.fst"), compress = 100)



