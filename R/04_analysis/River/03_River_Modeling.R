

# ══════════════════════════════════════════════════════════════════════════════
# HEADER
# ══════════════════════════════════════════════════════════════════════════════
#
# Fits beta regression models of egg-to-smolt survival for the Saint-Jean
# and Trinite rivers. For each river: linearity check (GAM with smooths),
# model selection by AICc (MuMIn::dredge), model averaging over the top set
# (delta AICc < 2), pseudo-R2, leave-one-out cross-validation, and
# model-averaged marginal effects with percentile bootstrap confidence
# intervals.
#
# Inputs:
#   - R/04_analysis/River/Final_Data/Final_DataSet_STJ.fst (FST)
#   - R/04_analysis/River/Final_Data/Final_DataSet_TRI.fst (FST)
#
# ══════════════════════════════════════════════════════════════════════════════


# ══════════════════════════════════════════════════════════════════════════════
# LIBRARIES
# ══════════════════════════════════════════════════════════════════════════════

library(dplyr)
library(mgcv)
library(fst)
library(ggplot2)
library(MuMIn)
library(betareg)
library(ggrepel)


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

stj <- read_fst(file.path(base_path, "R/04_analysis/River/Final_Data/Final_DataSet_STJ.fst"))
tri <- read_fst(file.path(base_path, "R/04_analysis/River/Final_Data/Final_DataSet_TRI.fst"))


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

# ── loo_cv ────────────────────────────────────────────────────────────────────
# Purpose: Leave-one-out cross-validation of a betareg model; refits the
#   model without each observation in turn and returns its out-of-sample
#   prediction.
# Parameters:
#   model  fitted betareg model
#   data   data used to fit the model

loo_cv <- function(model, data) {
  
  n     <- nrow(data)
  preds <- numeric(n)
  
  for (i in 1:n) {
    
    mod_i <- tryCatch(
      betareg(
        formula = formula(model),
        data    = data[-i, ],
        link    = "logit"
      ),
      error = function(e) NULL
    )
    
    if (!is.null(mod_i)) {
      preds[i] <- predict(mod_i, newdata = data[i, ], type = "response")
    } else {
      preds[i] <- NA
    }
  }
  
  return(preds)
}


# ── build_marginal_plot_data_avg ──────────────────────────────────────────────
# Purpose: Computes the model-averaged marginal effect of one predictor
#   across the top-set models. The focal predictor varies over its observed
#   range while all other predictors are held at their mean; predictions of
#   all top-set models are averaged with AICc weights. Confidence intervals
#   are percentile bootstrap intervals (predict.betareg() provides no
#   standard errors, so delta-method intervals cannot be computed).
# Parameters:
#   top_models    list of betareg models (top set, delta AICc < 2)
#   var           character; focal predictor
#   data_scaled   scaled data used to fit the models
#   data_orig     data in original units (observed points)
#   center        scaling centers (attr "scaled:center")
#   scale_params  scaling SDs (attr "scaled:scale")
#   response_var  character; response variable
#   n_points      integer; number of points in the prediction grid
#   n_boot        integer; number of bootstrap iterations

build_marginal_plot_data_avg <- function(top_models, var, data_scaled, data_orig,
                                         center, scale_params,
                                         response_var = "surv_egg_smolt",
                                         n_points = 10, n_boot = 20) {
  
  # AICc weights
  aicc_vals <- sapply(top_models, AICc)
  delta     <- aicc_vals - min(aicc_vals)
  weights   <- exp(-0.5 * delta) / sum(exp(-0.5 * delta))
  
  # Prediction grid
  var_seq <- seq(min(data_scaled[[var]], na.rm = TRUE),
                 max(data_scaled[[var]], na.rm = TRUE),
                 length.out = n_points)
  
  nd <- data_scaled[rep(1, n_points),
                    !names(data_scaled) %in% response_var,
                    drop = FALSE]
  for (col in names(nd)) nd[[col]] <- mean(data_scaled[[col]], na.rm = TRUE)
  nd[[var]] <- var_seq
  
  # Predictions of each model — ALL models are used for the weighted mean
  preds_by_model <- lapply(top_models, function(m) {
    as.numeric(predict(m, newdata = nd, type = "response"))
  })
  
  # Weighted mean (uses all models, including those without the variable)
  fit_avg <- Reduce("+", mapply(function(p, w) p * w,
                                preds_by_model, weights,
                                SIMPLIFY = FALSE))
  
  # Bootstrap CI
  # At each iteration: resampling, refit of all models,
  # recomputation of AICc weights, new weighted mean.
  # The uncertainty reflects both coefficient uncertainty AND
  # model selection uncertainty.
  boot_mat <- matrix(NA, nrow = n_points, ncol = n_boot)
  for (b in seq_len(n_boot)) {
    idx    <- sample(nrow(data_scaled), replace = TRUE)
    mods_b <- lapply(top_models, function(m) {
      tryCatch(
        betareg(formula(m), data = data_scaled[idx, ], link = "logit"),
        error = function(e) NULL
      )
    })
    valid  <- !sapply(mods_b, is.null)
    mods_b <- mods_b[valid]
    if (length(mods_b) == 0) next
    
    aicc_b  <- sapply(mods_b, AICc)
    delta_b <- aicc_b - min(aicc_b)
    wts_b   <- exp(-0.5 * delta_b) / sum(exp(-0.5 * delta_b))
    
    preds_b <- lapply(mods_b, function(m) {
      as.numeric(predict(m, newdata = nd, type = "response"))
    })
    boot_mat[, b] <- Reduce("+", mapply(function(p, w) p * w,
                                        preds_b, wts_b,
                                        SIMPLIFY = FALSE))
  }
  lower <- apply(boot_mat, 1, quantile, 0.025, na.rm = TRUE)
  upper <- apply(boot_mat, 1, quantile, 0.975, na.rm = TRUE)
  
  x_orig <- if (var == "cohort") var_seq else var_seq * scale_params[var] + center[var]
  
  # Identify which models contain the focal variable
  # (to plot only the relevant individual curves)
  contains_var <- sapply(top_models, function(m) {
    var %in% attr(terms(m), "term.labels")
  })
  
  list(
    var      = var,
    pred_avg = data.frame(x     = x_orig,
                          fit   = as.numeric(fit_avg),
                          lower = as.numeric(lower),
                          upper = as.numeric(upper)),
    # Only the models that contain the focal variable
    pred_models = lapply(which(contains_var), function(i) {
      vars_i <- paste(attr(terms(top_models[[i]]), "term.labels"), collapse = " + ")
      data.frame(
        x     = x_orig,
        fit   = preds_by_model[[i]],
        model = paste0("M", i, " (w=", round(weights[i], 2), "): ", vars_i)
      )
    }),
    obs_df = data.frame(
      x      = data_orig[[var]],
      y      = data_orig[[response_var]],
      cohort = data_orig$cohort
    )
  )
}


# ── build_marginal_ggplot_avg ─────────────────────────────────────────────────
# Purpose: Plots the output of build_marginal_plot_data_avg: model-averaged
#   curve with bootstrap CI (blue), individual top-set models containing the
#   focal predictor (grey), and observed values coloured and labelled by
#   cohort.

build_marginal_ggplot_avg <- function(plot_data, base_size = 14, point_size = 4,
                                      label_size = 3, line_width = 1.5,
                                      ribbon_alpha = 0.2,
                                      y_label = "Egg-to-smolt survival") {
  d          <- plot_data
  cohort_mid <- median(d$obs_df$cohort, na.rm = TRUE)
  
  # Individual model curves (grey, only those containing the variable)
  pred_models_df <- bind_rows(d$pred_models)
  
  p <- ggplot(d$pred_avg, aes(x = x, y = fit))
  
  # Add individual curves if at least one model contains the variable
  if (nrow(pred_models_df) > 0) {
    p <- p +
      geom_line(data = pred_models_df,
                aes(x = x, y = fit, group = model, linetype = model),
                color = "grey50", linewidth = 0.8)
  }
  
  p +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = ribbon_alpha, fill = "steelblue") +
    geom_line(linewidth = line_width, color = "steelblue4") +
    geom_point(data = d$obs_df, aes(x = x, y = y, fill = cohort),
               inherit.aes = FALSE,
               size = point_size, color = "black", shape = 21, stroke = 0.5) +
    scale_fill_gradient2(low = "#00897B", mid = "#BDBDBD", high = "#8E24AA",
                         midpoint = cohort_mid, name = "Cohort") +
    geom_text_repel(data = d$obs_df, aes(x = x, y = y, label = cohort),
                    inherit.aes = FALSE,
                    size = label_size, max.overlaps = 20) +
    labs(x = d$var, y = y_label) +
    theme_classic(base_size = base_size)
}


# ══════════════════════════════════════════════════════════════════════════════
# SCALE DATA
# ══════════════════════════════════════════════════════════════════════════════

# Save scaling parameters before creating scaled datasets.
# These must be saved to back-transform scaled values to original units for plotting.

stj_scale_params <- scale(stj %>% dplyr::select(-cohort, -surv_egg_smolt))
stj_center       <- attr(stj_scale_params, "scaled:center")
stj_scale        <- attr(stj_scale_params, "scaled:scale")

scaled_stj <- cbind(
  data.frame(
    cohort         = stj$cohort,
    surv_egg_smolt = stj$surv_egg_smolt
  ),
  as.data.frame(stj_scale_params)
)

tri_scale_params <- scale(tri %>% dplyr::select(-cohort, -surv_egg_smolt))
tri_center       <- attr(tri_scale_params, "scaled:center")
tri_scale        <- attr(tri_scale_params, "scaled:scale")

scaled_tri <- cbind(
  data.frame(
    cohort         = tri$cohort,
    surv_egg_smolt = tri$surv_egg_smolt
  ),
  as.data.frame(tri_scale_params)
)


# ══════════════════════════════════════════════════════════════════════════════
# ST. JEAN RIVER (STJ)
# ══════════════════════════════════════════════════════════════════════════════

# ── STJ — Linearity check (GAM with smooths) ──────────────────────────────────

# Testing the linearity of all variables in a GAM using smoothness.
# If edf ≈ 1 -> relationship is linear -> smooth is unnecessary.
# If edf > 1 -> non-linearity detected -> plot the smooth to evaluate.

stj_gam_smooth <- gam(
  surv_egg_smolt ~ s(cohort) +
    s(co_gr_median_temp) + s(em_median_temp) + s(em_var_temp) + 
    s(co_gr_dd_o_20_temp) + s(co_year_mean_sum_ptgi) + 
    s(co_gr_mean_disc) + 
    s(em_quant_75qt_disc) + 
    s(co_mean_forest),
  family = betar(link = "logit"),
  data = scaled_stj
)
summary(stj_gam_smooth)


# ── STJ — All available variables ─────────────────────────────────────────────

stj_vars <- c(
  "cohort",
  "co_gr_median_temp", "em_median_temp", "em_var_temp", 
  "co_gr_dd_o_20_temp", "co_year_mean_sum_ptgi", 
  "co_gr_mean_disc", 
  "em_quant_75qt_disc", 
  "co_mean_forest"
)
stj_formula <- as.formula(
  paste("surv_egg_smolt ~", paste(stj_vars, collapse = " + "))
)


# ── STJ — Fit beta regression ─────────────────────────────────────────────────

stj_glm <- betareg(
  formula   = stj_formula,
  data      = scaled_stj,
  link      = "logit",
  na.action = na.fail
)
summary(stj_glm)


# ── STJ — Dredge and model averaging ──────────────────────────────────────────

stj_all_models <- dredge(stj_glm, rank = "AICc")
head(stj_all_models, 10)

stj_top_models <- get.models(stj_all_models, subset = delta < 2)
stj_avg_model  <- model.avg(stj_top_models, full = TRUE)

summary(stj_avg_model)
stj_avg_model

# Relative importance of each variable
sw(stj_avg_model)

confint(stj_avg_model, full = TRUE) # Confidence intervals on coefficient estimates (not fitted values).
confint(stj_avg_model, full = FALSE) # In order to get confint for conditionnal average


# ── STJ — Pseudo-R2 ───────────────────────────────────────────────────────────

post_avg_stj <- bind_rows(
  lapply(stj_top_models, function(x) {
    
    summ_x <- summary(x)
    
    return(
      data.frame(
        model         = paste(as.character(x$formula)[-1], collapse = " ~ "),
        pseudo_r2     = summ_x$pseudo.r.squared
      )
    )
  })
)
post_avg_stj

mean(post_avg_stj$pseudo_r2)
range(post_avg_stj$pseudo_r2)


# ── STJ — Leave-one-out cross-validation (LOO-CV) ─────────────────────────────
# Since there are only 30 observations, LOO-CV is used to test for overfitting.

stj_loo_results <- lapply(seq_along(stj_top_models), function(i) {
  
  mod       <- stj_top_models[[i]]
  loo_preds <- loo_cv(mod, scaled_stj)
  observed  <- scaled_stj$surv_egg_smolt
  
  insample_preds <- fitted(mod)
  
  # RMSE: root mean squared error — lower is better
  # If LOO RMSE >> in-sample RMSE, the model overfits
  rmse_insample <- sqrt(mean((observed - insample_preds)^2, na.rm = TRUE))
  rmse_loo      <- sqrt(mean((observed - loo_preds)^2,      na.rm = TRUE))
  
  ss_res_in  <- sum((observed - insample_preds)^2, na.rm = TRUE)
  ss_res_loo <- sum((observed - loo_preds)^2,      na.rm = TRUE)
  ss_tot     <- sum((observed - mean(observed))^2)
  r2_insample <- 1 - ss_res_in  / ss_tot
  r2_loo      <- 1 - ss_res_loo / ss_tot
  
  vars <- paste(attr(terms(mod), "term.labels"), collapse = " + ")
  
  data.frame(
    model         = paste0("Model ", i, ": ", vars),
    rmse_insample = round(rmse_insample, 6),
    rmse_loo      = round(rmse_loo,      6),
    rmse_ratio    = round(rmse_loo / rmse_insample, 3),
    r2_insample   = round(r2_insample, 3),
    r2_loo        = round(r2_loo,      3)
  )
})

stj_loo_summary <- bind_rows(stj_loo_results)
stj_loo_summary

mean(stj_loo_summary$rmse_ratio, na.rm = TRUE)
range(stj_loo_summary$rmse_ratio, na.rm = TRUE)

mean(stj_loo_summary$r2_loo)
range(stj_loo_summary$r2_loo)

# Interpretation:
# rmse_ratio close to 1       -> no overfitting, model generalizes well
# rmse_ratio > 1.5            -> moderate overfitting
# rmse_ratio > 2              -> serious overfitting
# r2_loo close to r2_insample -> model generalizes well
# r2_loo much lower           -> model is overfit to training data


# ── STJ — Model-averaged marginal effects ─────────────────────────────────────

# All variables present in at least one top-set model
all_vars <- unique(unlist(lapply(stj_top_models, function(m) {
  attr(terms(m), "term.labels")
})))

set.seed(42)
avg_plot_data <- lapply(all_vars, function(v) {
  build_marginal_plot_data_avg(
    stj_top_models, v,
    scaled_stj, stj, stj_center, stj_scale,
    n_points = 10,
    n_boot   = 20
  )
})

avg_plots <- lapply(avg_plot_data, build_marginal_ggplot_avg)

# Blue curve = weighted mean
# Grey curves = individual models
avg_plots[[1]] +
  labs(title = paste("Variable :", all_vars[1]))

avg_plots[[2]] +
  labs(title = paste("Variable :", all_vars[2]))

avg_plots[[3]] +
  labs(title = paste("Variable :", all_vars[3]))

avg_plots[[4]] +
  labs(title = paste("Variable :", all_vars[4]))

avg_plots[[5]] +
  labs(title = paste("Variable :", all_vars[5]))

avg_plots[[6]] +
  labs(title = paste("Variable :", all_vars[6]))


# ══════════════════════════════════════════════════════════════════════════════
# TRINITE RIVER (TRI)
# ══════════════════════════════════════════════════════════════════════════════

# ── TRI — Linearity check (GAM with smooths) ──────────────────────────────────

tri_gam_smooth <- gam(
  surv_egg_smolt ~ s(cohort) +
    s(co_gr_median_temp) + s(co_gr_var_temp) + s(co_year_mean_sum_ptgi) + 
    s(em_median_temp) + s(em_var_temp) + s(co_gr_dd_u_16_temp) + 
    s(co_gr_mean_disc) + s(co_30d_mean_disch) + 
    s(co_fdd_bimon_ja_fe) + 
    s(co_mean_forest),
  family    = betar(link = "logit"),
  data      = scaled_tri,
  na.action = na.fail
)
summary(tri_gam_smooth)

par(mfrow = c(2, 2))
plot(tri_gam_smooth, select = 1,  shade = TRUE) 
plot(tri_gam_smooth, select = 5,  shade = TRUE) 
plot(tri_gam_smooth, select = 10,  shade = TRUE)
plot(tri_gam_smooth, select = 11, shade = TRUE) 
par(mfrow = c(1, 1))


# ── TRI — All available variables ─────────────────────────────────────────────

tri_vars <- c(
  "cohort",
  "co_gr_median_temp", "co_gr_var_temp", "co_year_mean_sum_ptgi",  
  "em_median_temp", "em_var_temp", "co_gr_dd_u_16_temp",  
  "co_gr_mean_disc", "co_30d_mean_disch", 
  "co_fdd_bimon_ja_fe", 
  "co_mean_forest"
)
tri_formula <- as.formula(
  paste("surv_egg_smolt ~", paste(tri_vars, collapse = " + "))
)


# ── TRI — Fit beta regression ─────────────────────────────────────────────────

tri_glm <- betareg(
  formula   = tri_formula,
  data      = scaled_tri,
  link      = "logit",
  na.action = na.fail
)
summary(tri_glm)


# ── TRI — Dredge and model averaging ──────────────────────────────────────────

tri_all_models <- dredge(tri_glm, rank = "AICc")
head(tri_all_models, 20)

tri_top_models <- get.models(tri_all_models, subset = delta < 2)
tri_avg_model  <- model.avg(tri_top_models, full = TRUE)

summary(tri_avg_model)
tri_avg_model

sw(tri_avg_model)

confint(tri_avg_model, full = TRUE)

length(tri_top_models)


# ── TRI — Pseudo-R2 ───────────────────────────────────────────────────────────

post_avg_tri <- bind_rows(
  lapply(tri_top_models, function(x){
    
    summ_x <- summary(x)
    
    return(
      data.frame(
        model         = paste(as.character(x$formula)[-1], collapse = " + "),
        pseudo_r2     = summ_x$pseudo.r.squared
      )
    )
  })
)
post_avg_tri

mean(post_avg_tri$pseudo_r2, na.rm = TRUE)
range(post_avg_tri$pseudo_r2, na.rm = TRUE)


# ── TRI — Leave-one-out cross-validation (LOO-CV) ─────────────────────────────

tri_loo_results <- lapply(seq_along(tri_top_models), function(i) {
  
  mod       <- tri_top_models[[i]]
  loo_preds <- loo_cv(mod, scaled_tri)
  observed  <- scaled_tri$surv_egg_smolt
  
  insample_preds <- fitted(mod)
  
  rmse_insample <- sqrt(mean((observed - insample_preds)^2, na.rm = TRUE))
  rmse_loo      <- sqrt(mean((observed - loo_preds)^2,      na.rm = TRUE))
  
  ss_res_in  <- sum((observed - insample_preds)^2, na.rm = TRUE)
  ss_res_loo <- sum((observed - loo_preds)^2,      na.rm = TRUE)
  ss_tot     <- sum((observed - mean(observed))^2)
  r2_insample <- 1 - ss_res_in  / ss_tot
  r2_loo      <- 1 - ss_res_loo / ss_tot
  
  vars <- paste(attr(terms(mod), "term.labels"), collapse = " + ")
  
  data.frame(
    model         = paste0("Model ", i, ": ", vars),
    rmse_insample = round(rmse_insample, 6),
    rmse_loo      = round(rmse_loo,      6),
    rmse_ratio    = round(rmse_loo / rmse_insample, 3),
    r2_insample   = round(r2_insample, 3),
    r2_loo        = round(r2_loo,      3)
  )
})

tri_loo_summary <- bind_rows(tri_loo_results)
tri_loo_summary

mean(tri_loo_summary$rmse_ratio,  na.rm = TRUE)
range(tri_loo_summary$rmse_ratio, na.rm = TRUE)
mean(tri_loo_summary$r2_loo,      na.rm = TRUE)
range(tri_loo_summary$r2_loo,     na.rm = TRUE)


# ── TRI — Model-averaged marginal effects ─────────────────────────────────────

all_vars_tri <- unique(unlist(lapply(tri_top_models, function(m) {
  attr(terms(m), "term.labels")
})))

set.seed(42)
avg_plot_data_tri <- lapply(all_vars_tri, function(v) {
  build_marginal_plot_data_avg(
    tri_top_models, v,
    scaled_tri, tri, tri_center, tri_scale,
    n_points = 10,
    n_boot   = 20
  )
})

avg_plots_tri <- lapply(avg_plot_data_tri, build_marginal_ggplot_avg)

# Blue curve = weighted mean 
# Grey curves = individual models
avg_plots_tri[[1]] +
  labs(title = paste("Variable :", all_vars_tri[1]))

avg_plots_tri[[2]] +
  labs(title = paste("Variable :", all_vars_tri[2]))

avg_plots_tri[[3]] +
  labs(title = paste("Variable :", all_vars_tri[3]))

avg_plots_tri[[4]] +
  labs(title = paste("Variable :", all_vars_tri[4]))

avg_plots_tri[[5]] +
  labs(title = paste("Variable :", all_vars_tri[5]))

avg_plots_tri[[6]] +
  labs(title = paste("Variable :", all_vars_tri[6]))

avg_plots_tri[[7]] +
  labs(title = paste("Variable :", all_vars_tri[7]))


