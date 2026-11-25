


# ══════════════════════════════════════════════════════════════════════════════
# HEADER
# ══════════════════════════════════════════════════════════════════════════════
#
# Explores correlations among marine predictors, screens them by VIF, and
# checks residual diagnostics (DHARMa) of beta regression models of the 1SW
# and 2SW return rates for the Saint-Jean and Trinite rivers, then exports the
# final datasets.
#
# Inputs:
#   - R/04_analysis/Sea/Final_Data/PreSelection_DataSet_1sw_stj.fst (FST)
#   - R/04_analysis/Sea/Final_Data/PreSelection_DataSet_2sw_stj.fst (FST)
#   - R/04_analysis/Sea/Final_Data/PreSelection_DataSet_1sw_tri.fst (FST)
#   - R/04_analysis/Sea/Final_Data/PreSelection_DataSet_2sw_tri.fst (FST)
#
# Outputs:
#   - R/04_analysis/Sea/Final_Data/Final_DataSet_STJ_1SW.fst (FST)
#   - R/04_analysis/Sea/Final_Data/Final_DataSet_STJ_2SW.fst (FST)
#   - R/04_analysis/Sea/Final_Data/Final_DataSet_TRI_1SW.fst (FST)
#   - R/04_analysis/Sea/Final_Data/Final_DataSet_TRI_2SW.fst (FST)
#
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

# ══════════════════════════════════════════════════════════════════════════════
# LOAD LIBRARIES
# ══════════════════════════════════════════════════════════════════════════════

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
library(MuMIn)


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

# ── plot_cor ──────────────────────────────────────────────────────────────────
# Purpose: Computes pairwise Pearson correlations among all numeric columns
#   (excluding cohort) and returns a colour-coded upper-triangle heatmap.
# Parameters:
#   data    data frame with numeric columns
#   title_  character; plot title

plot_cor <- function(data, title_ = ""){
  
  cor_mat <- cor(
    data %>% dplyr::select(where(is.numeric), -cohort),
    use = "pairwise.complete.obs"
  )
  
  cor_df <- as.data.frame(cor_mat) %>%
    rownames_to_column("var1") %>%
    pivot_longer(-var1, names_to = "var2", values_to = "Cor")
  
  cor_df <- cor_df %>%
    mutate(
      var1 = factor(var1, levels = colnames(cor_mat)),
      var2 = factor(var2, levels = colnames(cor_mat))
    ) %>%
    filter(as.numeric(var1) < as.numeric(var2))
  
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


# ── plot_pca ──────────────────────────────────────────────────────────────────
# Purpose: Runs a PCA on all numeric columns (excluding cohort and constant
#   columns) and returns a biplot with individual scores and variable arrows,
#   highlighting one target variable.
# Parameters:
#   data           data frame with numeric columns
#   highlight_var  character; name of the variable to highlight in red
#   arrow_scale    numeric; scaling factor applied to variable arrows
#   title_         character; plot title

plot_pca <- function(data, highlight_var = "surv_egg_smolt", arrow_scale = 3, title_ = "") {
  
  data <- data %>%
    dplyr::select(where(is.numeric), -cohort)
  
  keep <- sapply(data, function(x) length(unique(x)) != 1)
  data <- data[, keep]
  
  pca_res <- PCA(data, scale.unit = TRUE, graph = FALSE)
  
  ind_df <- as.data.frame(pca_res$ind$coord)
  ind_df$ind <- rownames(ind_df)
  
  var_df <- as.data.frame(pca_res$var$coord)
  var_df$var <- rownames(var_df)
  
  var_df$col <- ifelse(var_df$var == highlight_var, "red", "black")
  
  var_df_scaled <- var_df %>%
    mutate(
      Dim.1 = Dim.1 * arrow_scale,
      Dim.2 = Dim.2 * arrow_scale
    )
  
  ggplot() +
    geom_point(data = ind_df, aes(x = Dim.1, y = Dim.2), color = "grey30") +
    geom_segment(data = var_df_scaled,
                 aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2, color = col),
                 arrow = arrow(length = unit(0.3, "cm")), size = 0.7) +
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

# ══════════════════════════════════════════════════════════════════════════════
# IMPORT DATA
# ══════════════════════════════════════════════════════════════════════════════

stj_1sw <- read_fst(file.path(base_path, "R/04_analysis/Sea/Final_Data", "PreSelection_DataSet_1sw_stj.fst")) %>% na.omit()
stj_2sw <- read_fst(file.path(base_path, "R/04_analysis/Sea/Final_Data", "PreSelection_DataSet_2sw_stj.fst")) %>% na.omit()
tri_1sw <- read_fst(file.path(base_path, "R/04_analysis/Sea/Final_Data", "PreSelection_DataSet_1sw_tri.fst")) %>% na.omit()
tri_2sw <- read_fst(file.path(base_path, "R/04_analysis/Sea/Final_Data", "PreSelection_DataSet_2sw_tri.fst")) %>% na.omit()


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


# ══════════════════════════════════════════════════════════════════════════════
# CREATE GROUPS OF VARIABLES
# ══════════════════════════════════════════════════════════════════════════════

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


# ══════════════════════════════════════════════════════════════════════════════
# EXPLORE CORRELATIONS
# ══════════════════════════════════════════════════════════════════════════════

# Temperature - STJ
plot_pca(stj_1sw_temp, arrow_scale = 5, highlight_var = "return_1sw", title_ = "Temperature St. Jean")
plot_cor(data = stj_1sw_temp, title_ = "Temperature St. Jean")

plot_pca(stj_2sw_temp, arrow_scale = 5, highlight_var = "return_2sw", title_ = "Temperature St. Jean")
plot_cor(data = stj_2sw_temp, title_ = "Temperature St. Jean")

# Temperature - TRI
plot_pca(tri_1sw_temp, arrow_scale = 5, highlight_var = "return_1sw", title_ = "Temperature Trinite")
plot_cor(data = tri_1sw_temp, title_ = "Temperature Trinite")

plot_pca(tri_2sw_temp, arrow_scale = 5, highlight_var = "return_2sw", title_ = "Temperature Trinite")
plot_cor(data = tri_2sw_temp, title_ = "Temperature Trinite")

# No variable is removed a priori on the basis of correlations alone.



# ══════════════════════════════════════════════════════════════════════════════
# EXPLORE ALL VARIABLES CORRELATIONS
# ══════════════════════════════════════════════════════════════════════════════

# ── St. Jean — 1SW ────────────────────────────────────────────────────────────
plot_pca(stj_1sw, arrow_scale = 5, highlight_var = "return_1sw", title_ = "All variables - St. Jean (1sw)")
plot_cor(data = stj_1sw, title_ = "All variables - St. Jean (1sw)")


# ── St. Jean — 2SW ────────────────────────────────────────────────────────────
plot_pca(stj_2sw, arrow_scale = 5, highlight_var = "return_2sw", title_ = "All variables - St. Jean (2sw)")
plot_cor(data = stj_2sw, title_ = "All variables - St. Jean (2sw)")


# ── Trinite — 1SW ─────────────────────────────────────────────────────────────
plot_pca(tri_1sw, arrow_scale = 5, highlight_var = "return_1sw", title_ = "All variables - Trinite (1sw)")
plot_cor(data = tri_1sw, title_ = "All variables - Trinite (1sw)")

# ── Trinite — 2SW ─────────────────────────────────────────────────────────────
plot_pca(tri_2sw, arrow_scale = 5, highlight_var = "return_2sw", title_ = "All variables - Trinite (2sw)")
plot_cor(data = tri_2sw, title_ = "All variables - Trinite (2sw)")


# ══════════════════════════════════════════════════════════════════════════════
# SELECT VARIABLES BASED ON VIF: TEMPERATURE
# ══════════════════════════════════════════════════════════════════════════════

sort(vif(lm(return_1sw ~ cohort + delta_T_sst_minus_river + StJeanCorridor_TH  + BelleIsle_TH + LabradorEntrance_TH + LabradorSea_TH, data = stj_1sw)))
sort(vif(lm(return_2sw ~ cohort + delta_T_sst_minus_river + StJeanCorridor_TH  + BelleIsle_TH + LabradorEntrance_TH + LabradorSea_TH, data = stj_2sw)))
sort(vif(lm(return_1sw ~ cohort + delta_T_sst_minus_river + TriniteCorridor_TH + BelleIsle_TH + LabradorEntrance_TH + LabradorSea_TH, data = tri_1sw)))
sort(vif(lm(return_2sw ~ cohort + delta_T_sst_minus_river + TriniteCorridor_TH + BelleIsle_TH + LabradorEntrance_TH + LabradorSea_TH, data = tri_2sw)))
# No variable has a VIF > 10; none is removed.



# ══════════════════════════════════════════════════════════════════════════════
# SELECT VARIABLES BASED ON VIF: ALL VARIABLES
# ══════════════════════════════════════════════════════════════════════════════

# ── St. Jean — 1SW ────────────────────────────────────────────────────────────
sort(vif(lm(return_1sw ~ cohort + delta_T_sst_minus_river + Plankton_t +
              StJeanCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
              LabradorSea_TH + smolt_k + er1NAC_lag1 + pfa_1sw, data = stj_1sw)))
# pfa_1sw is removed

sort(vif(lm(return_1sw ~ cohort + delta_T_sst_minus_river + Plankton_t +
              StJeanCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
              LabradorSea_TH + smolt_k + er1NAC_lag1, data = stj_1sw)))


# ── Trinite — 1SW ─────────────────────────────────────────────────────────────
sort(vif(lm(return_1sw ~ cohort + delta_T_sst_minus_river + Plankton_t +
              TriniteCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
              LabradorSea_TH + smolt_k + er1NAC_lag1 + pfa_1sw, data = tri_1sw)))
# Although the usual rule is to remove the variable with the highest VIF,
# pfa_1sw was excluded. Both alternatives are compared below.


sort(vif(lm(return_1sw ~ cohort + delta_T_sst_minus_river + Plankton_t +
              TriniteCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
              LabradorSea_TH + smolt_k + er1NAC_lag1, data = tri_1sw))) # Without pfa_1sw
sort(vif(lm(return_1sw ~ cohort + delta_T_sst_minus_river + Plankton_t +
              TriniteCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
              LabradorSea_TH + smolt_k + pfa_1sw, data = tri_1sw))) # Without er1NAC_lag1
# Retained model: without pfa_1sw.


# ── St. Jean — 2SW ────────────────────────────────────────────────────────────
# Exploitation is represented by er_index (PCA of er_WG_lag1 and er2NAC_lag2;
# see 01_Combine_Covariates_Sea.R)
sort(vif(lm(return_2sw ~ cohort + delta_T_sst_minus_river + Plankton_t +
              StJeanCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
              LabradorSea_TH + smolt_k + er_index + pfa_2sw, data = stj_2sw)))
# pfa_2sw has the highest VIF, followed by cohort.

# Without pfa_2sw
sort(vif(lm(return_2sw ~ cohort + delta_T_sst_minus_river + Plankton_t +
              StJeanCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
              LabradorSea_TH + smolt_k + er_index, data = stj_2sw)))
# pfa_2sw is removed.


# ── Trinite — 2SW ─────────────────────────────────────────────────────────────
sort(vif(lm(return_2sw ~ cohort + delta_T_sst_minus_river + Plankton_t +
              TriniteCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
              LabradorSea_TH + smolt_k + er_index + pfa_2sw, data = tri_2sw)))
# pfa_2sw has the highest VIF

# Without pfa_2sw
sort(vif(lm(return_2sw ~ cohort + delta_T_sst_minus_river + Plankton_t +
              TriniteCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
              LabradorSea_TH + smolt_k + er_index, data = tri_2sw)))

# Without er_index
sort(vif(lm(return_2sw ~ cohort + delta_T_sst_minus_river + Plankton_t +
              TriniteCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
              LabradorSea_TH + smolt_k, data = tri_2sw)))
# Retained model


# ══════════════════════════════════════════════════════════════════════════════
# DEFINE FORMULAS FOR MODELS
# ══════════════════════════════════════════════════════════════════════════════

# STJ — 1SW
f_stj_1sw <- as.formula(
  "return_1sw ~ cohort + delta_T_sst_minus_river + Plankton_t +
   StJeanCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
   LabradorSea_TH + smolt_k + er1NAC_lag1"
)

# TRI — 1SW
f_tri_1sw <- as.formula(
  "return_1sw ~ cohort + delta_T_sst_minus_river + Plankton_t +
   TriniteCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
   LabradorSea_TH + smolt_k + er1NAC_lag1"
)

# STJ — 2SW
f_stj_2sw <- as.formula(
  "return_2sw ~ cohort + delta_T_sst_minus_river + Plankton_t +
   StJeanCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
   LabradorSea_TH + smolt_k + er_index"
)

# TRI — 2SW
f_tri_2sw <- as.formula(
  "return_2sw ~ cohort + delta_T_sst_minus_river + Plankton_t +
   TriniteCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
   LabradorSea_TH + smolt_k"
)


# ══════════════════════════════════════════════════════════════════════════════
# EXPLORE DISTRIBUTION
# ══════════════════════════════════════════════════════════════════════════════

par(mfrow = c(1,1))

# ── St. Jean — 1SW ────────────────────────────────────────────────────────────
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
# No deviation


# ── Trinite — 1SW ─────────────────────────────────────────────────────────────
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
# No deviation


# ── St. Jean — 2SW ────────────────────────────────────────────────────────────
summary(stj_2sw$return_2sw)

ggplot(stj_2sw, aes(x = return_2sw)) +
  geom_histogram(aes(y = after_stat(density)), bins = 10, fill = "grey80", color = "grey30") +
  geom_density(color = "black", linewidth = 1) +
  labs(title = "St. Jean - Distribution of return rate (2SW)", x = "Return rate", y = "Density") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())


# ── Scale predictors (cohort kept in years) ───────────────────────────────────
stj_2sw_scale_params <- scale(stj_2sw %>% dplyr::select(-cohort, -return_2sw, -river))
stj_2sw_center       <- attr(stj_2sw_scale_params, "scaled:center")
stj_2sw_scale        <- attr(stj_2sw_scale_params, "scaled:scale")

scaled_stj_2sw <- cbind(
  data.frame(
    cohort     = stj_2sw$cohort,
    return_2sw = stj_2sw$return_2sw
  ),
  as.data.frame(stj_2sw_scale_params)
)

f_stj_2sw <- as.formula(
  "return_2sw ~ cohort + delta_T_sst_minus_river + Plankton_t +
   StJeanCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
   LabradorSea_TH + smolt_k + er_index"
)


# ── Residual diagnostics of the global model ──────────────────────────────────
stj_2sw_beta <- gam(f_stj_2sw, family = betar(link = "logit"), data = scaled_stj_2sw)
stj_2sw_beta$converged

set.seed(1)
stj_2sw_sim <- simulateResiduals(stj_2sw_beta, n = 1000)

plot(stj_2sw_sim)
testDispersion(stj_2sw_sim)

# KS (P = 0.95), dispersion (P = 0.52) and outlier (P = 1) tests are
# fine: the beta distribution is appropriate. However, residuals vs. predicted
# (rank-transformed) show a wave-shaped pattern and the combined quantile test
# is significant. Three possible causes are checked:
#   (a) temporal autocorrelation;
#   (b) a missing nonlinear effect of one predictor;
#   (c) a pattern specific to the global model.

# ── Check if the pattern is sensitive to the predictor scale ──────────────────
# plot() tests quantiles against RANK-transformed predictions; check both.

testQuantiles(stj_2sw_sim, predictor = rank(stj_2sw_sim$fittedPredictedResponse),
              plot = TRUE)
plotResiduals(stj_2sw_sim, rank = FALSE)

# Significant on the rank scale (P = 0.018), but not on raw
# predictions. The pattern is real (not simulation noise: it persists with
# 1000 simulations) but weak, and depends on how predictions are spread.


# ── Check for temporal autocorrelation ────────────────────────────────────────
testTemporalAutocorrelation(stj_2sw_sim, time = stj_2sw_beta$model$cohort,
                            plot = TRUE)

# Durbin-Watson DW = 1.90, P = 0.78.
# -> No temporal autocorrelation (cause a ruled out).


# ── Check for nonlinear effect of a predictor ─────────────────────────────────
# Residuals against each predictor: a single predictor carrying the
# curvature would indicate a missing nonlinear term.

mf_stj_2sw <- stj_2sw_beta$model
par(mfrow = c(3, 3))
for (v in names(mf_stj_2sw)[-1]) {
  plotResiduals(stj_2sw_sim, form = mf_stj_2sw[[v]], xlab = v, main = v)
}
par(mfrow = c(1, 1))

# No significant problem for any predictor (slight, non-significant
# curvature for Plankton_t, StJeanCorridor_TH and LabradorSea_TH).
# -> No single predictor explains the pattern (cause b ruled out).


# The pattern seems specific to the global model (cause c). If so, the global
# model should be heavily penalized by AICc, and the retained models should
# be free of it. This is checked after model selection (below).


# ── Model selection ───────────────────────────────────────────────────────────
stj_2sw_glm <- betareg(
  formula   = f_stj_2sw,
  data      = scaled_stj_2sw,
  link      = "logit",
  na.action = na.fail
)
stj_2sw_all_models <- dredge(stj_2sw_glm, rank = "AICc")


# ── Check for weight of the global model in the selection ─────────────────────
all_tab <- as.data.frame(stj_2sw_all_models)
all_tab[all_tab$df == max(all_tab$df), c("df", "AICc", "delta", "weight")]

# delta = 16.6, weight = 1.8e-05: the global model is excluded from
# the top set and contributes nothing to model-averaged estimates.


# ── Check for residuals of all models in the top set (delta AICc < 2) ─────────
check_top <- function(top_models, data, nsim = 1000) {
  res <- lapply(seq_along(top_models), function(i) {
    f   <- formula(top_models[[i]])
    fit <- gam(f, family = betar(link = "logit"), data = data)
    set.seed(1)
    s   <- simulateResiduals(fit, n = nsim)
    data.frame(
      model          = i,
      formula        = paste(deparse(f[[3]]), collapse = ""),
      converged      = fit$converged,
      KS             = testUniformity(s, plot = FALSE)$p.value,
      dispersion     = testDispersion(s, plot = FALSE)$p.value,
      quantiles_rank = testQuantiles(s, predictor = rank(s$fittedPredictedResponse),
                                     plot = FALSE)$p.value,
      DW             = testTemporalAutocorrelation(s, time = data$cohort,
                                                   plot = FALSE)$p.value
    )
  })
  out <- do.call(rbind, res)
  out[, 4:7] <- round(out[, 4:7], 3)
  out
}

top_stj_2sw <- get.models(stj_2sw_all_models, subset = delta < 2)
check_top(top_stj_2sw, data = scaled_stj_2sw)

# Main residual plot of the best model
best_stj_2sw <- gam(formula(top_stj_2sw[[1]]), family = betar(link = "logit"),
                    data = scaled_stj_2sw)
set.seed(1)
plot(simulateResiduals(best_stj_2sw, n = 1000))

# All top-set models converge and show no problem (all P > 0.05),
# including the best model (quantiles on rank scale P = 0.85), except model 7
# (cohort + er_index + LabradorEntrance_TH + Plankton_t + smolt_k; quantiles
# on rank scale P = 0.001). Models with er_index alone (model 5, P = 0.77) or
# LabradorEntrance_TH alone (model 3, P = 0.40) show no deviation.
# -> The pattern arises when both predictors are in the same model, as in
#    the global model. Checked below.


# ── Check top-set model showing the deviation (model 7) ───────────────────────
# Are the two predictors related?
cor.test(stj_2sw$er_index, stj_2sw$LabradorEntrance_TH) # not strongly related

# Weight of model 7 within the top set
top_tab <- as.data.frame(subset(stj_2sw_all_models, delta < 2))
top_tab$weight_top <- top_tab$weight / sum(top_tab$weight)
top_tab[, c("df", "AICc", "delta", "weight_top")]

# Sensitivity: model-averaged estimates with vs. without model 7
avg_all <- model.avg(top_stj_2sw)
avg_no7 <- model.avg(top_stj_2sw[-7])

round(cbind(full_all = coef(avg_all, full = TRUE),
            full_no7 = coef(avg_no7, full = TRUE)[names(coef(avg_all, full = TRUE))]), 4)

summary(avg_all)$coefmat.full
summary(avg_no7)$coefmat.full

# The two predictors are moderately correlated (r = -0.40,
# P = 0.028). Model 7 has a weight of 0.09 within the top set. Excluding it
# does not change model-averaged estimates or their significance: cohort and
# Plankton_t remain significant (P < 0.001 and P = 0.005; estimates change
# < 2 %), smolt_k remains non-significant in the full average (P = 0.086 vs
# 0.095), and er_index remains far from significance (P = 0.71 vs 0.79).


# ── Final conclusion ──────────────────────────────────────────────────────────
# The quantile deviation appears only in models including both er_index and
# LabradorEntrance_TH (global model and top-set model 7, w = 0.09). It does
# not affect inference: the global model is excluded by AICc, all other
# top-set models show no problem, and excluding model 7 leaves
# model-averaged estimates and their significance unchanged. No change to the
# analysis is required.

# ── Trinite — 2SW ─────────────────────────────────────────────────────────────
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

# ── Full residual diagnostics: global model + top set ─────────────────────────
# Scales predictors (cohort kept in years), checks the global model, runs
# dredge, and checks every model of the top set (delta AICc < 2).
run_diagnostics <- function(data, formula, label, nsim = 1000) {
  pred <- setdiff(all.vars(formula)[-1], "cohort")
  d    <- data
  d[pred] <- lapply(d[pred], function(x) as.numeric(scale(x)))
  
  # Global model
  g <- gam(formula, family = betar(link = "logit"), data = d)
  set.seed(1)
  sg <- simulateResiduals(g, n = nsim)
  plot(sg)
  global_diag <- data.frame(
    converged      = g$converged,
    KS             = testUniformity(sg, plot = FALSE)$p.value,
    dispersion     = testDispersion(sg, plot = FALSE)$p.value,
    quantiles_rank = testQuantiles(sg, predictor = rank(sg$fittedPredictedResponse),
                                   plot = FALSE)$p.value,
    DW             = testTemporalAutocorrelation(sg, time = g$model$cohort,
                                                 plot = FALSE)$p.value
  )
  
  # Model selection
  b   <- do.call("betareg", list(formula = formula, data = d,
                                 link = "logit", na.action = "na.fail"))
  dd  <- dredge(b, rank = "AICc")
  tab <- as.data.frame(dd)
  global_sel <- tab[tab$df == max(tab$df), c("df", "AICc", "delta", "weight")]
  
  # Top set
  top <- get.models(dd, subset = delta < 2)
  
  cat("\n==========", label, "==========\n")
  cat("\n-- Global model: residual tests --\n"); print(round(global_diag[, -1], 3))
  cat("\n-- Global model: selection --\n");      print(global_sel)
  cat("\n-- Top set (delta AICc < 2) --\n");     print(check_top(top, data = d))
  
  invisible(list(dredge = dd, top = top, data_scaled = d))
}

diag_tri_2sw <- run_diagnostics(tri_2sw, f_tri_2sw, "TRI 2SW")

d_tri  <- diag_tri_2sw$data_scaled
top_tri <- diag_tri_2sw$top

# Does smolt_k carry the residual pattern of model 1?
fit1 <- gam(formula(top_tri[[1]]), family = betar(link = "logit"), data = d_tri)
set.seed(1)
sim1 <- simulateResiduals(fit1, n = 1000)
plot(sim1)
plotResiduals(sim1, form = d_tri$smolt_k, xlab = "smolt_k (not in model 1)")

# Weights within the top set
top_tab_tri <- as.data.frame(subset(diag_tri_2sw$dredge, delta < 2))
top_tab_tri$weight_top <- top_tab_tri$weight / sum(top_tab_tri$weight)
top_tab_tri[, c("df", "AICc", "delta", "weight_top")]

# Sensitivity: model-averaged estimates with vs. without model 1
avg_tri_all <- model.avg(top_tri)
avg_tri_no1 <- model.avg(top_tri[-1])
summary(avg_tri_all)$coefmat.full
summary(avg_tri_no1)$coefmat.full

# Wave in the lower quantile of the global model (quantiles on rank
# scale P = 0.005); global model excluded by AICc (delta = 13.4,
# w = 1.4e-04). In the top set, the best model (BelleIsle_TH + Plankton_t,
# w = 0.44) shows the same deviation (P = 0.001), whereas the two models
# adding smolt_k show none (P = 0.73 and 0.96). Residuals of the best model
# show a significant negative trend against smolt_k: the deviation reflects
# the omitted effect of smolt condition. No temporal autocorrelation
# (DW P = 0.13-0.15).
#
# Excluding the best model does not change conclusions: BelleIsle_TH
# (P = 0.014 vs 0.017) and Plankton_t (P < 0.001) remain significant with
# similar estimates; smolt_k remains non-significant in the full average
# (P = 0.42 vs 0.12). No change to the analysis is required.


# ══════════════════════════════════════════════════════════════════════════════
# TESTING FOR DATA INDEPENDENCE
# ══════════════════════════════════════════════════════════════════════════════


# ── St. Jean — 1SW ────────────────────────────────────────────────────────────
stj_1sw_ind <- betareg(f_stj_1sw, data = stj_1sw, link = "logit")

par(mfrow = c(1, 1))
acf(residuals(stj_1sw_ind))
lmtest::dwtest(stj_1sw_ind)
# No substantial temporal autocorrelation


# ── Trinite — 1SW ─────────────────────────────────────────────────────────────
tri_1sw_ind <- betareg(f_tri_1sw, data = tri_1sw, link = "logit")

par(mfrow = c(1, 1))
acf(residuals(tri_1sw_ind))
lmtest::dwtest(tri_1sw_ind)
# No substantial temporal autocorrelation


# ── St. Jean — 2SW ────────────────────────────────────────────────────────────
stj_2sw_ind <- betareg(f_stj_2sw, data = stj_2sw, link = "logit")

par(mfrow = c(1, 1))
acf(residuals(stj_2sw_ind))
lmtest::dwtest(stj_2sw_ind)
# No substantial temporal autocorrelation


# ── Trinite — 2SW ─────────────────────────────────────────────────────────────
tri_2sw_ind <- betareg(f_tri_2sw, data = tri_2sw, link = "logit")

par(mfrow = c(1, 1))
acf(residuals(tri_2sw_ind))
lmtest::dwtest(tri_2sw_ind)
# Slight temporal autocorrelation, considered acceptable


# ══════════════════════════════════════════════════════════════════════════════
# CHECK FOR OVERDISPERSION AND HETEROSCEDASTICITY
# ══════════════════════════════════════════════════════════════════════════════

# ── St. Jean — 1SW ────────────────────────────────────────────────────────────
mean(stj_1sw$return_1sw)
var(stj_1sw$return_1sw)

plot(resid(stj_1sw_beta)); abline(h = 0, col = "red", lty = 2)
sum(resid(stj_1sw_beta, type = "pearson")^2) / stj_1sw_beta$df.resid


# ── Trinite — 1SW ─────────────────────────────────────────────────────────────
mean(tri_1sw$return_1sw)
var(tri_1sw$return_1sw)

plot(resid(tri_1sw_beta)); abline(h = 0, col = "red", lty = 2)
sum(resid(tri_1sw_beta, type = "pearson")^2) / tri_1sw_beta$df.resid


# ── St. Jean — 2SW ────────────────────────────────────────────────────────────
mean(stj_2sw$return_2sw)
var(stj_2sw$return_2sw)

plot(resid(stj_2sw_beta)); abline(h = 0, col = "red", lty = 2)
sum(resid(stj_2sw_beta, type = "pearson")^2) / stj_2sw_beta$df.resid


# ── Trinite — 2SW ─────────────────────────────────────────────────────────────
mean(tri_2sw$return_2sw)
var(tri_2sw$return_2sw)

plot(resid(tri_2sw_beta)); abline(h = 0, col = "red", lty = 2)
sum(resid(tri_2sw_beta, type = "pearson")^2) / tri_2sw_beta$df.resid

# ══════════════════════════════════════════════════════════════════════════════
# EXPORT FINAL DATA SETS
# ══════════════════════════════════════════════════════════════════════════════

write.fst(stj_1sw, file.path(base_path, "R/04_analysis/Sea/Final_Data/Final_DataSet_STJ_1SW.fst"))
write.fst(stj_2sw, file.path(base_path, "R/04_analysis/Sea/Final_Data/Final_DataSet_STJ_2SW.fst"))
write.fst(tri_1sw, file.path(base_path, "R/04_analysis/Sea/Final_Data/Final_DataSet_TRI_1SW.fst"))
write.fst(tri_2sw, file.path(base_path, "R/04_analysis/Sea/Final_Data/Final_DataSet_TRI_2SW.fst"))

