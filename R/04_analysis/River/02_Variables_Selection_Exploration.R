

# ══════════════════════════════════════════════════════════════════════════════
# HEADER
# ══════════════════════════════════════════════════════════════════════════════
#
# Selects predictor variables via VIF screening across temperature, discharge,
# and FDD groups, then fits beta regression models for egg-to-smolt survival
# for the Saint-Jean and Trinite rivers and exports the final cleaned datasets.
#
# Inputs:
#   - R/04_analysis/River/Final_Data/PreVariablesSelection_River_StJean.fst (FST)
#   - R/04_analysis/River/Final_Data/PreVariablesSelection_River_Trinite.fst (FST)
#
# Outputs:
#   - R/04_analysis/River/Final_Data/Final_DataSet_STJ.fst (FST)
#   - R/04_analysis/River/Final_Data/Final_DataSet_TRI.fst (FST)
#
# ══════════════════════════════════════════════════════════════════════════════


# ══════════════════════════════════════════════════════════════════════════════
# LIBRARIES
# ══════════════════════════════════════════════════════════════════════════════

library(dplyr)
library(tidyr)
library(mgcv)
library(fst)
library(car)
library(FactoMineR)
library(ggrepel)
library(DHARMa)
library(ggplot2)
library(ggpubr)
library(lmtest)
library(tibble)
library(betareg)


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
# DATA IMPORT
# ══════════════════════════════════════════════════════════════════════════════

# ── Pre-selection dataset — Saint-Jean ────────────────────────────────────────
# SOURCE: Merged predictor and response variables
stj <- read_fst(
  file.path(base_path,
            "R/04_analysis/River/Final_Data/PreVariablesSelection_River_StJean.fst"))


# ── Pre-selection dataset — Trinite ───────────────────────────────────────────
# SOURCE: Merged predictor and response variables
tri <- read_fst(
  file.path(base_path,
            "R/04_analysis/River/Final_Data/PreVariablesSelection_River_Trinite.fst"))


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

# ── plot_pca ──────────────────────────────────────────────────────────────────
# Purpose: Runs a PCA on all numeric columns and returns a biplot with
#   individual scores and variable arrows, highlighting one target variable.
# Parameters:
#   data           data frame with numeric columns
#   highlight_var  character; name of the variable to highlight in red
#   arrow_scale    numeric; scaling factor applied to variable arrows
#   title_         character; plot title

plot_pca <- function(data, highlight_var = "surv_egg_smolt",
                     arrow_scale = 3, title_ = "") {
  
  data    <- data %>% dplyr::select(where(is.numeric))
  pca_res <- PCA(data, scale.unit = TRUE, graph = FALSE)
  
  ind_df      <- as.data.frame(pca_res$ind$coord)
  ind_df$ind  <- rownames(ind_df)
  
  var_df      <- as.data.frame(pca_res$var$coord)
  var_df$var  <- rownames(var_df)
  var_df$col  <- ifelse(var_df$var == highlight_var, "red", "black")
  
  var_df_scaled <- var_df %>%
    mutate(Dim.1 = Dim.1 * arrow_scale,
           Dim.2 = Dim.2 * arrow_scale)
  
  ggplot() +
    geom_point(data = ind_df,
               aes(x = Dim.1, y = Dim.2), color = "grey30") +
    geom_segment(data = var_df_scaled,
                 aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2,
                     color = col),
                 arrow = arrow(length = unit(0.3, "cm")), size = 0.7) +
    geom_text_repel(data = var_df_scaled,
                    aes(x = Dim.1, y = Dim.2, label = var, color = col),
                    size = 4, box.padding = 0.8, point.padding = 0.8,
                    max.overlaps = Inf, segment.color = "grey50") +
    scale_color_identity() +
    labs(x = paste0("Dim 1 (", round(pca_res$eig[1, 2], 1), "%)"),
         y = paste0("Dim 2 (", round(pca_res$eig[2, 2], 1), "%)"),
         title = title_) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    coord_fixed()
}


# ── plot_cor ──────────────────────────────────────────────────────────────────
# Purpose: Computes pairwise Pearson correlations among all numeric columns
#   (excluding cohort) and returns a colour-coded upper-triangle heatmap.
# Parameters:
#   data    data frame with numeric columns
#   title_  character; plot title

plot_cor <- function(data, title_ = "") {
  
  cor_mat <- cor(
    data %>% dplyr::select(where(is.numeric), -cohort),
    use = "pairwise.complete.obs"
  )
  
  cor_df <- as.data.frame(cor_mat) %>%
    rownames_to_column("var1") %>%
    pivot_longer(-var1, names_to = "var2", values_to = "Cor") %>%
    mutate(
      var1 = factor(var1, levels = colnames(cor_mat)),
      var2 = factor(var2, levels = colnames(cor_mat))
    ) %>%
    filter(as.numeric(var1) < as.numeric(var2))
  
  ggplot(cor_df, aes(var1, var2, fill = Cor)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(Cor, 2)), size = 3, color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-1, 1)) +
    scale_y_discrete(limits = rev) +
    labs(title = title_) +
    theme_minimal() +
    theme(axis.text.x      = element_text(angle = 45, hjust = 1),
          axis.title        = element_blank(),
          panel.grid        = element_blank(),
          plot.title        = element_text(hjust = 0.5)) +
    coord_fixed()
}



# ══════════════════════════════════════════════════════════════════════════════
# MAIN SCRIPT
# ══════════════════════════════════════════════════════════════════════════════

# ── Rescale survival and remove unwanted columns ──────────────────────────────

# Survival was exported as a percentage; rescale to [0, 1] for beta regression
if (any(stj$surv_egg_smolt > 1)) {
  stj <- stj %>% mutate(surv_egg_smolt = surv_egg_smolt / 100)
}
if (any(tri$surv_egg_smolt > 1)) {
  tri <- tri %>% mutate(surv_egg_smolt = surv_egg_smolt / 100)
}

# Columns excluded a priori before VIF screening

# 1 - PGTI is 100% correlated with PGTI by construction
cols_opposite_pgti <- c("em_sum_inv_pgti", "co_year_mean_sum_inv_ptgi")

# Retain only co_fdd_per_oc_ja and co_fdd_bimon_ja_fe from the FDD family;
# all other FDD metrics are removed to avoid redundancy before VIF screening
cols_fdd <- c(
  "co_fdd_per_oc_no", "co_fdd_per_oc_de", "co_fdd_per_oc_fe",
  "co_fdd_per_oc_ma", "co_fdd_per_oc_ap",
  "co_fdd_mon_oc",    "co_fdd_mon_no",    "co_fdd_mon_de",
  "co_fdd_mon_ja",    "co_fdd_mon_fe",    "co_fdd_mon_ma",
  "co_fdd_mon_ap",    "co_fdd_bimon_no_de", "co_fdd_bimon_ma_ap"
)

cols_to_remove <- c(cols_var_disc, cols_opposite_pgti, cols_fdd)

clean_dataframe <- function(df) {
  df %>%
    dplyr::filter(!is.na(surv_egg_smolt)) %>%
    dplyr::select(where(~ !all(is.na(.)))) %>%
    dplyr::select(where(~ !all(. == 0))) %>%
    dplyr::select(-any_of(cols_to_remove)) %>%
    na.omit()
}

stj <- clean_dataframe(stj)
tri <- clean_dataframe(tri)


# ── VIF screening — temperature variables ────────────────────────────────────

# Temperature variables were divided into two biological periods:
#   1. Emergence — short annual window; simple summary statistics (mean,
#      median) are adequate. Cumulative indices such as PGTI or CDD are
#      considered less informative over such a short period.
#      em_mean_temp and em_median_temp were therefore prioritised over
#      em_sum_ptgi and em_dd_u_16_temp.
#   2. Growth — long season (approx. April to October); cumulative variables
#      (co_year_mean_sum_ptgi, co_gr_dd_u_16_temp, co_gr_dd_o_20_temp) are
#      considered more biologically relevant than mean or median temperatures,
#      because temperature varies substantially within the season and simple
#      summary statistics may obscure meaningful thermal signals.
#
# No variables were removed solely on the VIF. Variables representing a
# similar ecological mechanism were treated as interchangeable descriptors,
# and the VIF was used only to retain a single representative metric.

# Full temperature model
sort(vif(lm(surv_egg_smolt ~ cohort +
              co_gr_mean_temp + co_gr_median_temp + co_gr_var_temp +
              em_mean_temp + em_median_temp + em_var_temp +
              co_gr_dd_u_16_temp + co_gr_dd_o_20_temp +
              em_sum_ptgi + co_year_mean_sum_ptgi,
            data = stj)))
sort(vif(lm(surv_egg_smolt ~ cohort +
              co_gr_mean_temp + co_gr_median_temp + co_gr_var_temp +
              em_mean_temp + em_median_temp + em_var_temp +
              co_gr_dd_u_16_temp + co_gr_dd_o_20_temp +
              em_sum_ptgi + co_year_mean_sum_ptgi,
            data = tri)))
# em_sum_ptgi is removed since it has the highest VIF

# Without em_sum_ptgi
sort(vif(lm(surv_egg_smolt ~ cohort +
              co_gr_mean_temp + co_gr_median_temp + co_gr_var_temp +
              em_mean_temp + em_median_temp + em_var_temp +
              co_gr_dd_u_16_temp + co_gr_dd_o_20_temp +
              co_year_mean_sum_ptgi,
            data = stj)))
sort(vif(lm(surv_egg_smolt ~ cohort +
              co_gr_mean_temp + co_gr_median_temp + co_gr_var_temp +
              em_mean_temp + em_median_temp + em_var_temp +
              co_gr_dd_u_16_temp + co_gr_dd_o_20_temp +
              co_year_mean_sum_ptgi,
            data = tri)))
# co_gr_mean_temp showed the highest VIF 

# Without co_gr_mean_temp
sort(vif(lm(surv_egg_smolt ~ cohort +
              co_gr_median_temp + co_gr_var_temp +
              em_mean_temp + em_median_temp + em_var_temp +
              co_gr_dd_u_16_temp + co_gr_dd_o_20_temp +
              co_year_mean_sum_ptgi,
            data = stj)))
sort(vif(lm(surv_egg_smolt ~ cohort +
              co_gr_median_temp + co_gr_var_temp +
              em_mean_temp + em_median_temp + em_var_temp +
              co_gr_dd_u_16_temp + co_gr_dd_o_20_temp +
              co_year_mean_sum_ptgi,
            data = tri)))
# co_gr_var_temp had the highest VIF in STJ
# em_mean_temp had the highest VIF in TRI

# Without co_gr_var_temp in STJ and without em_mean_temp in TRI
sort(vif(lm(surv_egg_smolt ~ cohort +
              co_gr_median_temp +
              em_mean_temp + em_median_temp + em_var_temp +
              co_gr_dd_u_16_temp + co_gr_dd_o_20_temp +
              co_year_mean_sum_ptgi,
            data = stj)))
sort(vif(lm(surv_egg_smolt ~ cohort +
              co_gr_median_temp + co_gr_var_temp +
              em_median_temp + em_var_temp +
              co_gr_dd_u_16_temp + co_gr_dd_o_20_temp +
              co_year_mean_sum_ptgi,
            data = tri)))
# Remove em_mean_temp in STJ
# TRI stays like that

# Without em_mean_temp in STJ, no change in TRI
sort(vif(lm(surv_egg_smolt ~ cohort +
              co_gr_median_temp +
              em_median_temp + em_var_temp +
              co_gr_dd_u_16_temp + co_gr_dd_o_20_temp +
              co_year_mean_sum_ptgi,
            data = stj)))
sort(vif(lm(surv_egg_smolt ~ cohort +
              co_gr_median_temp + co_gr_var_temp +
              em_median_temp + em_var_temp +
              co_gr_dd_u_16_temp + co_gr_dd_o_20_temp +
              co_year_mean_sum_ptgi,
            data = tri)))

# Correlation check — temperature
stj_temp <- stj %>%
  dplyr::select(surv_egg_smolt, cohort, co_gr_median_temp, em_median_temp,
                em_var_temp, co_gr_dd_u_16_temp, co_gr_dd_o_20_temp,
                  co_year_mean_sum_ptgi)
ggarrange(
  plot_pca(data = stj_temp, highlight_var = "surv_egg_smolt",
           title_ = "St. Jean R. - Temperature"),
  plot_cor(data = stj_temp, title_ = "St. Jean R. - Temperature"),
  nrow = 1
)

tri_temp <- tri %>%
  dplyr::select(surv_egg_smolt, cohort, co_gr_median_temp, co_gr_var_temp, 
                em_median_temp, em_var_temp, co_gr_dd_u_16_temp, co_gr_dd_o_20_temp,
                co_year_mean_sum_ptgi)
ggarrange(
  plot_pca(data = tri_temp, highlight_var = "surv_egg_smolt",
           title_ = "Trinite R. - Temperature"),
  plot_cor(data = tri_temp, title_ = "Trinite R. - Temperature"),
  nrow = 1
)


# ── VIF screening — discharge variables ──────────────────────────────────────

# Discharge variables were also divided into emergence and growth periods.
# Emergence data for Trinite River are not available due to missing data.

# Full discharge model
sort(vif(lm(surv_egg_smolt ~ cohort +
              co_gr_mean_disc + co_gr_median_disc +
              co_gr_mean_75qt_disc + co_30d_mean_disch +
              em_mean_disc + em_median_disc + em_quant_75qt_disc,
            data = stj)))
sort(vif(lm(surv_egg_smolt ~ cohort +
              co_gr_mean_disc + co_gr_median_disc +
              co_gr_mean_75qt_disc + co_30d_mean_disch,
            data = tri)))
# em_mean_disc removed from St. Jean (highest VIF).
# Trinite unchanged (em_mean_disc not available).

# Without em_mean_disc
sort(vif(lm(surv_egg_smolt ~ cohort +
              co_gr_mean_disc + co_gr_median_disc +
              co_gr_mean_75qt_disc + co_30d_mean_disch +
              em_median_disc + em_quant_75qt_disc,
            data = stj)))
sort(vif(lm(surv_egg_smolt ~ cohort +
              co_gr_mean_disc + co_gr_median_disc +
              co_gr_mean_75qt_disc + co_30d_mean_disch,
            data = tri)))
# co_gr_mean_75qt_disc above 10 in St. Jean; no VIF > 10 in Trinite.
# Removed from both rivers for consistency.

# Without co_gr_mean_75qt_disc
sort(vif(lm(surv_egg_smolt ~ cohort +
              co_gr_mean_disc + co_gr_median_disc +
              co_30d_mean_disch + em_median_disc + em_quant_75qt_disc,
            data = stj)))
sort(vif(lm(surv_egg_smolt ~ cohort +
              co_gr_mean_disc + co_gr_median_disc + co_30d_mean_disch,
            data = tri)))


# Correlation check — discharge
stj_disc <- stj %>%
  dplyr::select(surv_egg_smolt, cohort,
                co_gr_mean_disc, co_gr_median_disc, co_30d_mean_disch,
                em_median_disc, em_quant_75qt_disc)
ggarrange(
  plot_pca(data = stj_disc, highlight_var = "surv_egg_smolt",
           title_ = "St. Jean R. - Discharge"),
  plot_cor(data = stj_disc, title_ = "St. Jean R. - Discharge"),
  nrow = 1
)

tri_disc <- tri %>%
  dplyr::select(surv_egg_smolt, cohort,
                co_gr_mean_disc, co_gr_median_disc, co_30d_mean_disch)
ggarrange(
  plot_pca(data = tri_disc, highlight_var = "surv_egg_smolt",
           title_ = "Trinite R. - Discharge"),
  plot_cor(data = tri_disc, title_ = "Trinite R. - Discharge"),
  nrow = 1
)


# ── VIF screening — FDD variables ────────────────────────────────────────────

sort(vif(lm(surv_egg_smolt ~ cohort +
              co_fdd_per_oc_ja + co_fdd_bimon_ja_fe,
            data = stj)))
sort(vif(lm(surv_egg_smolt ~ cohort +
              co_fdd_per_oc_ja + co_fdd_bimon_ja_fe,
            data = tri)))
# co_fdd_per_oc_ja is removed


# ── Build clean datasets for combined VIF screening ──────────────────────────

stj_clean <- stj %>%
  dplyr::select(
    cohort, surv_egg_smolt,
    co_gr_median_temp, em_median_temp, em_var_temp,
    co_gr_dd_u_16_temp, co_gr_dd_o_20_temp, co_year_mean_sum_ptgi,
    co_gr_mean_disc, co_gr_median_disc, co_30d_mean_disch,
    em_median_disc, em_quant_75qt_disc,
    co_fdd_bimon_ja_fe,
    co_mean_forest
  )

tri_clean <- tri %>%
  dplyr::select(
    cohort, surv_egg_smolt,
    co_gr_median_temp, co_gr_var_temp, co_year_mean_sum_ptgi,
    em_median_temp, em_var_temp, co_gr_dd_u_16_temp, co_gr_dd_o_20_temp,
    co_gr_mean_disc, co_gr_median_disc, co_30d_mean_disch,
    co_fdd_bimon_ja_fe,
    co_mean_forest
  )

# ── VIF screening — all remaining variables combined ─────────────────────────

# Full combined model
sort(vif(lm(surv_egg_smolt ~ cohort +
              co_gr_median_temp + em_median_temp + em_var_temp + 
              co_gr_dd_u_16_temp + co_gr_dd_o_20_temp + co_year_mean_sum_ptgi + 
              co_gr_mean_disc + co_gr_median_disc + co_30d_mean_disch + 
              em_median_disc + em_quant_75qt_disc + 
              co_fdd_bimon_ja_fe + 
              co_mean_forest,
            data = stj_clean)))
sort(vif(lm(surv_egg_smolt ~ cohort +
              co_gr_median_temp + co_gr_var_temp + co_year_mean_sum_ptgi + 
              em_median_temp + em_var_temp + co_gr_dd_u_16_temp + co_gr_dd_o_20_temp + 
              co_gr_mean_disc + co_gr_median_disc + co_30d_mean_disch + 
              co_fdd_bimon_ja_fe + 
              co_mean_forest,
              data = tri_clean)))
# co_gr_dd_u_16_temp is removed in STJ
# co_gr_dd_o_20_temp is removed in TRI


sort(vif(lm(surv_egg_smolt ~ cohort +
              co_gr_median_temp + em_median_temp + em_var_temp + 
              co_gr_dd_o_20_temp + co_year_mean_sum_ptgi + 
              co_gr_mean_disc + co_gr_median_disc + co_30d_mean_disch + 
              em_median_disc + em_quant_75qt_disc + 
              co_fdd_bimon_ja_fe + 
              co_mean_forest,
            data = stj_clean)))
sort(vif(lm(surv_egg_smolt ~ cohort +
              co_gr_median_temp + co_gr_var_temp + co_year_mean_sum_ptgi + 
              em_median_temp + em_var_temp + co_gr_dd_u_16_temp + 
              co_gr_mean_disc + co_gr_median_disc + co_30d_mean_disch + 
              co_fdd_bimon_ja_fe + 
              co_mean_forest,
            data = tri_clean)))
# em_median_disc removed from STJ
# Nothing is removed from TRI for now because no em_median_disc in TRI

sort(vif(lm(surv_egg_smolt ~ cohort +
              co_gr_median_temp + em_median_temp + em_var_temp + 
              co_gr_dd_o_20_temp + co_year_mean_sum_ptgi + 
              co_gr_mean_disc + co_gr_median_disc + co_30d_mean_disch + 
              em_quant_75qt_disc + 
              co_fdd_bimon_ja_fe + 
              co_mean_forest,
            data = stj_clean)))
sort(vif(lm(surv_egg_smolt ~ cohort +
              co_gr_median_temp + co_gr_var_temp + co_year_mean_sum_ptgi + 
              em_median_temp + em_var_temp + co_gr_dd_u_16_temp + 
              co_gr_mean_disc + co_gr_median_disc + co_30d_mean_disch + 
              co_fdd_bimon_ja_fe + 
              co_mean_forest,
            data = tri_clean)))
# co_gr_median_disc removed in both

sort(vif(lm(surv_egg_smolt ~ cohort +
              co_gr_median_temp + em_median_temp + em_var_temp + 
              co_gr_dd_o_20_temp + co_year_mean_sum_ptgi + 
              co_gr_mean_disc + co_30d_mean_disch + 
              em_quant_75qt_disc + 
              co_fdd_bimon_ja_fe + 
              co_mean_forest,
            data = stj_clean)))
sort(vif(lm(surv_egg_smolt ~ cohort +
              co_gr_median_temp + co_gr_var_temp + co_year_mean_sum_ptgi + 
              em_median_temp + em_var_temp + co_gr_dd_u_16_temp + 
              co_gr_mean_disc + co_30d_mean_disch + 
              co_fdd_bimon_ja_fe + 
              co_mean_forest,
            data = tri_clean)))
# co_30d_mean_disch is removed from STJ
# TRI stays like that

sort(vif(lm(surv_egg_smolt ~ cohort +
              co_gr_median_temp + em_median_temp + em_var_temp + 
              co_gr_dd_o_20_temp + co_year_mean_sum_ptgi + 
              co_gr_mean_disc + 
              em_quant_75qt_disc + 
              co_fdd_bimon_ja_fe + 
              co_mean_forest,
            data = stj_clean)))
# co_fdd_bimon_ja_fe removed from STJ

sort(vif(lm(surv_egg_smolt ~ cohort +
              co_gr_median_temp + em_median_temp + em_var_temp + 
              co_gr_dd_o_20_temp + co_year_mean_sum_ptgi + 
              co_gr_mean_disc + 
              em_quant_75qt_disc + 
              co_mean_forest,
            data = stj_clean)))


# Correlation check — all remaining variables
stj_clean_plot <- stj_clean %>%
  dplyr::select(surv_egg_smolt, cohort,
                co_gr_median_temp, em_median_temp, em_var_temp, 
                  co_gr_dd_o_20_temp, co_year_mean_sum_ptgi, 
                  co_gr_mean_disc, em_quant_75qt_disc, 
                  co_mean_forest)
ggarrange(
  plot_pca(data = stj_clean_plot, highlight_var = "surv_egg_smolt",
           title_ = "St. Jean R. - Remaining variables"),
  plot_cor(data = stj_clean_plot,
           title_ = "St. Jean R. - Remaining variables"),
  nrow = 1
)

tri_clean_plot <- tri_clean %>%
  dplyr::select(surv_egg_smolt, cohort,
                co_gr_median_temp, co_gr_var_temp, co_year_mean_sum_ptgi, 
                  em_median_temp, em_var_temp, co_gr_dd_u_16_temp, 
                  co_gr_mean_disc, co_30d_mean_disch, 
                  co_fdd_bimon_ja_fe, co_mean_forest)
ggarrange(
  plot_pca(data = tri_clean_plot, highlight_var = "surv_egg_smolt",
           title_ = "Trinite R. - Remaining variables"),
  plot_cor(data = tri_clean_plot,
           title_ = "Trinite R. - Remaining variables"),
  nrow = 1
)


# ── Test for temporal autocorrelation ─────────────────────────────────────────

par(mfrow = c(1, 1))

stj_ind <- betareg(
  surv_egg_smolt ~ cohort +
    co_gr_median_temp + em_median_temp + em_var_temp + 
    co_gr_dd_o_20_temp + co_year_mean_sum_ptgi + 
    co_gr_mean_disc + 
    em_quant_75qt_disc + 
    co_mean_forest,
  data = stj_clean, 
  link = "logit"
)
acf(residuals(stj_ind))
lmtest::dwtest(stj_ind)

tri_ind <- betareg(
  surv_egg_smolt ~ cohort +
    co_gr_median_temp + co_gr_var_temp + co_year_mean_sum_ptgi + 
    em_median_temp + em_var_temp + co_gr_dd_u_16_temp + 
    co_gr_mean_disc + co_30d_mean_disch + 
    co_fdd_bimon_ja_fe + 
    co_mean_forest,
  data = tri_clean, 
  link = "logit"
)
acf(residuals(tri_ind))
lmtest::dwtest(tri_ind)

# No evidence of temporal autocorrelation in residuals. Observations were
# therefore treated as independent.


# ── Explore response variable distribution ────────────────────────────────────

# Saint-Jean
summary(stj_clean$surv_egg_smolt)  # No 0, no 1

ggplot(stj_clean, aes(x = surv_egg_smolt)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 10, fill = "grey80", color = "grey30") +
  geom_density(color = "black", linewidth = 1) +
  labs(title = "St. Jean - Distribution of Egg-to-Smolt Survival",
       x = "Egg-to-Smolt survival", y = "Density") +
  theme_minimal(base_size = 13) +
  theme(plot.title        = element_text(face = "bold"),
        panel.grid.minor  = element_blank())

stj_beta <- gam(
  surv_egg_smolt ~ cohort +
    co_gr_median_temp + em_median_temp + em_var_temp + 
    co_gr_dd_o_20_temp + co_year_mean_sum_ptgi + 
    co_gr_mean_disc + 
    em_quant_75qt_disc + 
    co_mean_forest,
  family = betar(link = "logit"), data = stj_clean
)

stj_sim <- simulateResiduals(stj_beta)
plot(stj_sim)
testDispersion(stj_sim)

# Trinite
summary(tri_clean$surv_egg_smolt)  # No 0, no 1

ggplot(tri_clean, aes(x = surv_egg_smolt)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 10, fill = "grey80", color = "grey30") +
  geom_density(color = "black", linewidth = 1) +
  labs(title = "Trinite - Distribution of Egg-to-Smolt Survival",
       x = "Egg-to-Smolt survival", y = "Density") +
  theme_minimal(base_size = 13) +
  theme(plot.title        = element_text(face = "bold"),
        panel.grid.minor  = element_blank())

tri_beta <- gam(
  surv_egg_smolt ~ cohort +
    co_gr_median_temp + co_gr_var_temp + co_year_mean_sum_ptgi + 
    em_median_temp + em_var_temp + co_gr_dd_u_16_temp + 
    co_gr_mean_disc + co_30d_mean_disch + 
    co_fdd_bimon_ja_fe + 
    co_mean_forest,
  family = betar(link = "logit"), data = tri_clean
)

tri_sim <- simulateResiduals(tri_beta)
plot(tri_sim)


# ── Check for overdispersion and heteroscedasticity ───────────────────────────

# Saint-Jean
mean(stj_clean$surv_egg_smolt)
var(stj_clean$surv_egg_smolt)
testDispersion(stj_beta)
plot(resid(stj_beta))
# Ratio ~1 = no overdispersion; >1.2-1.5 = suspicious; <0.8 = underdispersion
sum(resid(stj_beta, type = "pearson")^2) / stj_beta$df.resid

# Trinite
mean(tri_clean$surv_egg_smolt)
var(tri_clean$surv_egg_smolt)
testDispersion(tri_beta)
plot(resid(tri_beta))
sum(resid(tri_beta, type = "pearson")^2) / tri_beta$df.resid


# ── Export final datasets ─────────────────────────────────────────────────────

stj_clean_final <- stj_clean %>%
  dplyr::select(surv_egg_smolt, cohort,
                co_gr_median_temp, em_median_temp, em_var_temp, 
                co_gr_dd_o_20_temp, co_year_mean_sum_ptgi, 
                co_gr_mean_disc, em_quant_75qt_disc, 
                co_mean_forest)

tri_clean_final <- tri_clean %>%
  dplyr::select(surv_egg_smolt, cohort,
                co_gr_median_temp, co_gr_var_temp, co_year_mean_sum_ptgi, 
                em_median_temp, em_var_temp, co_gr_dd_u_16_temp, 
                co_gr_mean_disc, co_30d_mean_disch, 
                co_fdd_bimon_ja_fe, co_mean_forest)

write.fst(stj_clean_final,
          file.path(base_path, "R/04_analysis/River/Final_Data/Final_DataSet_STJ.fst"))
write.fst(tri_clean_final,
          file.path(base_path, "R/04_analysis/River/Final_Data/Final_DataSet_TRI.fst"))



