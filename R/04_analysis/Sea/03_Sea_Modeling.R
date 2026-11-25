
# ══════════════════════════════════════════════════════════════════════════════
# HEADER
# ══════════════════════════════════════════════════════════════════════════════
#
# Fits beta regression models of smolt-to-adult return rates (1SW and 2SW)
# for the St. Jean and Trinite rivers. For each river and sea age:
# linearity check (GAM with smooths), model selection by AICc
# (MuMIn::dredge), model averaging over the top set (delta AICc < 2),
# pseudo-R2, leave-one-out cross-validation, and model-averaged marginal
# effects with percentile bootstrap confidence intervals.
#
# Inputs:
#   - R/04_analysis/Sea/Final_Data/Final_DataSet_STJ_1SW.fst (FST)
#   - R/04_analysis/Sea/Final_Data/Final_DataSet_STJ_2SW.fst (FST)
#   - R/04_analysis/Sea/Final_Data/Final_DataSet_TRI_1SW.fst (FST)
#   - R/04_analysis/Sea/Final_Data/Final_DataSet_TRI_2SW.fst (FST)
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

stj_1sw <- read_fst(file.path(base_path, "R/04_analysis/Sea/Final_Data/Final_DataSet_STJ_1SW.fst")) %>% dplyr::select(-river)
stj_2sw <- read_fst(file.path(base_path, "R/04_analysis/Sea/Final_Data/Final_DataSet_STJ_2SW.fst")) %>% dplyr::select(-river)
tri_1sw <- read_fst(file.path(base_path, "R/04_analysis/Sea/Final_Data/Final_DataSet_TRI_1SW.fst")) %>% dplyr::select(-river)
tri_2sw <- read_fst(file.path(base_path, "R/04_analysis/Sea/Final_Data/Final_DataSet_TRI_2SW.fst")) %>% dplyr::select(-river)


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

# ── loo_cv ────────────────────────────────────────────────────────────────────
# Purpose: Leave-one-out cross-validation of the competing betareg models
#   (no variable-dispersion model here); refits the model without each
#   observation in turn and returns its out-of-sample prediction.
# Parameters:
#   model  fitted betareg model
#   data   data used to fit the model

loo_cv <- function(model, data) {
  
  n     <- nrow(data)
  preds <- numeric(n)
  
  for (i in 1:n) {
    
    mod_i <- tryCatch(
      betareg(
        formula  = formula(model),
        data     = data[-i, ],
        link     = "logit"
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
                                         response_var,
                                         n_points = 100, n_boot = 500) {
  
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
                                      y_label = "Return rate") {
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

# ── St. Jean — 1SW ────────────────────────────────────────────────────────────

stj_1sw_scale_params <- scale(stj_1sw %>% dplyr::select(-cohort, -return_1sw))
stj_1sw_center       <- attr(stj_1sw_scale_params, "scaled:center")
stj_1sw_scale        <- attr(stj_1sw_scale_params, "scaled:scale")

scaled_stj_1sw <- cbind(
  data.frame(
    cohort     = stj_1sw$cohort,
    return_1sw = stj_1sw$return_1sw
  ),
  as.data.frame(stj_1sw_scale_params)
)

# ── St. Jean — 2SW ────────────────────────────────────────────────────────────

stj_2sw_scale_params <- scale(stj_2sw %>% dplyr::select(-cohort, -return_2sw))
stj_2sw_center       <- attr(stj_2sw_scale_params, "scaled:center")
stj_2sw_scale        <- attr(stj_2sw_scale_params, "scaled:scale")

scaled_stj_2sw <- cbind(
  data.frame(
    cohort     = stj_2sw$cohort,
    return_2sw = stj_2sw$return_2sw
  ),
  as.data.frame(stj_2sw_scale_params)
)

# ── Trinite — 1SW ─────────────────────────────────────────────────────────────

tri_1sw_scale_params <- scale(tri_1sw %>% dplyr::select(-cohort, -return_1sw))
tri_1sw_center       <- attr(tri_1sw_scale_params, "scaled:center")
tri_1sw_scale        <- attr(tri_1sw_scale_params, "scaled:scale")

scaled_tri_1sw <- cbind(
  data.frame(
    cohort     = tri_1sw$cohort,
    return_1sw = tri_1sw$return_1sw
  ),
  as.data.frame(tri_1sw_scale_params)
)

# ── Trinite — 2SW ─────────────────────────────────────────────────────────────

tri_2sw_scale_params <- scale(tri_2sw %>% dplyr::select(-cohort, -return_2sw))
tri_2sw_center       <- attr(tri_2sw_scale_params, "scaled:center")
tri_2sw_scale        <- attr(tri_2sw_scale_params, "scaled:scale")

scaled_tri_2sw <- cbind(
  data.frame(
    cohort     = tri_2sw$cohort,
    return_2sw = tri_2sw$return_2sw
  ),
  as.data.frame(tri_2sw_scale_params)
)


# ══════════════════════════════════════════════════════════════════════════════
# ST. JEAN RIVER — 1SW
# ══════════════════════════════════════════════════════════════════════════════

# ── STJ 1SW — Linearity check (GAM with smooths) ──────────────────────────────

stj_1sw_gam_smooth <- gam(
  return_1sw ~ s(cohort) +
    s(delta_T_sst_minus_river) +
    s(StJeanCorridor_TH) + s(BelleIsle_TH) + s(LabradorEntrance_TH) + s(LabradorSea_TH) +
    s(Plankton_t) +
    s(smolt_k) +
    s(er1NAC_lag1),
  family = betar(link = "logit"),
  data = scaled_stj_1sw
)
summary(stj_1sw_gam_smooth)

par(mfrow = c(2, 2))
plot(stj_1sw_gam_smooth, select = 2, shade = TRUE)
plot(stj_1sw_gam_smooth, select = 3, shade = TRUE)
plot(stj_1sw_gam_smooth, select = 4, shade = TRUE)
par(mfrow = c(1, 1))

# They are basically all linear


# ── STJ 1SW — Fit beta regression ─────────────────────────────────────────────

f_stj_1sw <- as.formula(
  "return_1sw ~ cohort + delta_T_sst_minus_river + Plankton_t +
   StJeanCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
   LabradorSea_TH + smolt_k + er1NAC_lag1"
)

stj_1sw_glm <- betareg(
  formula   = f_stj_1sw,
  data      = scaled_stj_1sw,
  link      = "logit",
  na.action = na.fail
)
summary(stj_1sw_glm)


# ── STJ 1SW — Dredge and model averaging ──────────────────────────────────────

stj_1sw_all_models <- dredge(stj_1sw_glm, rank = "AICc")
head(stj_1sw_all_models)

stj_1sw_top_models <- get.models(stj_1sw_all_models, subset = delta < 2)
stj_1sw_avg_model  <- model.avg(stj_1sw_top_models, full = TRUE)

summary(stj_1sw_avg_model)
stj_1sw_avg_model

length(stj_1sw_top_models)

# Relative importance of each variable (sum of AICc weights)
sw(stj_1sw_avg_model)
confint(stj_1sw_avg_model, full = TRUE)


# ── STJ 1SW — Pseudo-R2 ───────────────────────────────────────────────────────

post_avg_stj_1sw <- bind_rows(
  lapply(stj_1sw_top_models, function(x) {
    
    summ_x  <- summary(x)
    
    data.frame(
      model         = paste(as.character(x$formula)[-1], collapse = " ~ "),
      pseudo_R2    = summ_x$pseudo.r.squared
    )
  })
)

post_avg_stj_1sw
mean(post_avg_stj_1sw$pseudo_R2)
range(post_avg_stj_1sw$pseudo_R2)


# ── STJ 1SW — Leave-one-out cross-validation (LOO-CV) ─────────────────────────

stj_1sw_loo_results <- lapply(seq_along(stj_1sw_top_models), function(i) {
  
  mod            <- stj_1sw_top_models[[i]]
  loo_preds      <- loo_cv(mod, scaled_stj_1sw)
  observed       <- scaled_stj_1sw$return_1sw
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

stj_1sw_loo_summary <- bind_rows(stj_1sw_loo_results)
stj_1sw_loo_summary

mean(stj_1sw_loo_summary$rmse_ratio)
range(stj_1sw_loo_summary$rmse_ratio)

mean(stj_1sw_loo_summary$r2_loo)
range(stj_1sw_loo_summary$r2_loo)


# ── STJ 1SW — Model-averaged marginal effects ─────────────────────────────────

# All variables present in at least one top-set model
all_vars <- unique(unlist(lapply(stj_1sw_top_models, function(m) {
  attr(terms(m), "term.labels")
})))

set.seed(42)
avg_plot_data <- lapply(all_vars, function(v) {
  build_marginal_plot_data_avg(
    stj_1sw_top_models, v,
    scaled_stj_1sw, stj_1sw, stj_1sw_center, stj_1sw_scale,
    response_var = "return_1sw",
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
# TRINITE RIVER — 1SW
# ══════════════════════════════════════════════════════════════════════════════

# ── TRI 1SW — Linearity check (GAM with smooths) ──────────────────────────────

tri_1sw_gam_smooth_er <- gam(
  return_1sw ~ s(cohort) +
    s(delta_T_sst_minus_river) +
    s(TriniteCorridor_TH) + s(BelleIsle_TH) + s(LabradorEntrance_TH) + s(LabradorSea_TH) +
    s(Plankton_t) +
    s(smolt_k) +
    s(er1NAC_lag1),
  family = betar(link = "logit"),
  data = scaled_tri_1sw
)
summary(tri_1sw_gam_smooth_er)

par(mfrow = c(1, 2))
plot(tri_1sw_gam_smooth_er, select = 2, shade = TRUE)
plot(tri_1sw_gam_smooth_er, select = 8, shade = TRUE)
par(mfrow = c(1, 1))

# linear


# ── TRI 1SW — Fit beta regression ─────────────────────────────────────────────

f_tri_1sw_er <- as.formula(
  "return_1sw ~ cohort + delta_T_sst_minus_river + Plankton_t +
   TriniteCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
   LabradorSea_TH + smolt_k + er1NAC_lag1"
)


tri_1sw_glm_er <- betareg(
  formula   = f_tri_1sw_er,
  data      = scaled_tri_1sw,
  link      = "logit",
  na.action = na.fail
)
summary(tri_1sw_glm_er)


# ── TRI 1SW — Dredge and model averaging ──────────────────────────────────────

tri_1sw_all_models_er <- dredge(tri_1sw_glm_er, rank = "AICc")
head(tri_1sw_all_models_er, 15)

tri_1sw_top_models_er <- get.models(tri_1sw_all_models_er, subset = delta < 2)
tri_1sw_avg_model_er  <- model.avg(tri_1sw_top_models_er)

summary(tri_1sw_avg_model_er)
length(tri_1sw_top_models_er)

sw(tri_1sw_avg_model_er)
confint(tri_1sw_avg_model_er, full = TRUE)
confint(tri_1sw_avg_model_er)            # conditional average


# ── TRI 1SW — Pseudo-R2 ───────────────────────────────────────────────────────

post_avg_tri_1sw_er <- bind_rows(
  lapply(tri_1sw_top_models_er, function(x) {
    data.frame(
      model         = paste(attr(terms(x), "term.labels"), collapse = " + "),
      pseudo_R2     = summary(x)$pseudo.r.squared
    )
  })
)

post_avg_tri_1sw_er
mean(post_avg_tri_1sw_er$pseudo_R2)
range(post_avg_tri_1sw_er$pseudo_R2)


# ── TRI 1SW — Leave-one-out cross-validation (LOO-CV) ─────────────────────────

tri_1sw_loo_results_er <- lapply(seq_along(tri_1sw_top_models_er), function(i) {
  
  mod            <- tri_1sw_top_models_er[[i]]
  loo_preds      <- loo_cv(mod, scaled_tri_1sw)          # TRI data
  observed       <- scaled_tri_1sw$return_1sw
  insample_preds <- fitted(mod)
  
  ss_res_in  <- sum((observed - insample_preds)^2, na.rm = TRUE)
  ss_res_loo <- sum((observed - loo_preds)^2,      na.rm = TRUE)
  ss_tot     <- sum((observed - mean(observed))^2)
  
  rmse_insample <- sqrt(mean((observed - insample_preds)^2, na.rm = TRUE))
  rmse_loo      <- sqrt(mean((observed - loo_preds)^2,      na.rm = TRUE))
  
  data.frame(
    model         = paste0("Model ", i, ": ",
                           paste(attr(terms(mod), "term.labels"), collapse = " + ")),
    rmse_insample = round(rmse_insample, 6),
    rmse_loo      = round(rmse_loo,      6),
    rmse_ratio    = round(rmse_loo / rmse_insample, 3),
    r2_insample   = round(1 - ss_res_in  / ss_tot, 3),
    r2_loo        = round(1 - ss_res_loo / ss_tot, 3)
  )
})

tri_1sw_loo_summary_er <- bind_rows(tri_1sw_loo_results_er)
tri_1sw_loo_summary_er

mean(tri_1sw_loo_summary_er$rmse_ratio)
range(tri_1sw_loo_summary_er$rmse_ratio)
mean(tri_1sw_loo_summary_er$r2_loo)
range(tri_1sw_loo_summary_er$r2_loo)


# ── TRI 1SW — Model-averaged marginal effects ─────────────────────────────────

all_vars_tri_1sw_er <- unique(unlist(lapply(tri_1sw_top_models_er, function(m) {
  attr(terms(m), "term.labels")
})))

set.seed(42)
avg_plot_data_tri_1sw_er <- lapply(all_vars_tri_1sw_er, function(v) {
  build_marginal_plot_data_avg(
    tri_1sw_top_models_er, v,
    scaled_tri_1sw, tri_1sw, tri_1sw_center, tri_1sw_scale,
    response_var = "return_1sw",
    n_points = 10,
    n_boot   = 20
  )
})

avg_plots_tri_1sw_er <- lapply(
  avg_plot_data_tri_1sw_er,
  function(x) build_marginal_ggplot_avg(x, y_label = "Return rate (1SW)")
)


avg_plots_tri_1sw_er[[1]] +
  labs(title = paste("Variable :", all_vars_tri_1sw_er[1]))

avg_plots_tri_1sw_er[[2]] +
  labs(title = paste("Variable :", all_vars_tri_1sw_er[2]))

avg_plots_tri_1sw_er[[3]] +
  labs(title = paste("Variable :", all_vars_tri_1sw_er[3]))

avg_plots_tri_1sw_er[[4]] +
  labs(title = paste("Variable :", all_vars_tri_1sw_er[4]))

avg_plots_tri_1sw_er[[5]] +
  labs(title = paste("Variable :", all_vars_tri_1sw_er[5]))

avg_plots_tri_1sw_er[[6]] +
  labs(title = paste("Variable :", all_vars_tri_1sw_er[6]))

avg_plots_tri_1sw_er[[7]] +
  labs(title = paste("Variable :", all_vars_tri_1sw_er[7]))


# ══════════════════════════════════════════════════════════════════════════════
# ST. JEAN RIVER — 2SW
# ══════════════════════════════════════════════════════════════════════════════

# ── STJ 2SW — Linearity check (GAM with smooths) ──────────────────────────────

stj_2sw_gam_smooth_er <- gam(
  return_2sw ~ s(cohort) +
    s(delta_T_sst_minus_river) +
    s(StJeanCorridor_TH) + s(BelleIsle_TH) + s(LabradorEntrance_TH) + s(LabradorSea_TH) +
    s(Plankton_t) +
    s(smolt_k) +
    s(er_index),
  family = betar(link = "logit"),
  data = scaled_stj_2sw
)
summary(stj_2sw_gam_smooth_er)

plot(stj_2sw_gam_smooth_er, select = 3, shade = TRUE)
# linear


# ── STJ 2SW — Fit beta regression ─────────────────────────────────────────────

f_stj_2sw_er <- as.formula(
  "return_2sw ~ cohort + delta_T_sst_minus_river + Plankton_t +
   StJeanCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
   LabradorSea_TH + smolt_k + er_index"
)

stj_2sw_glm_er <- betareg(
  formula   = f_stj_2sw_er,
  data      = scaled_stj_2sw,
  link      = "logit",
  na.action = na.fail
)
summary(stj_2sw_glm_er)


# ── STJ 2SW — Dredge and model averaging ──────────────────────────────────────

stj_2sw_all_models_er <- dredge(stj_2sw_glm_er, rank = "AICc")
head(stj_2sw_all_models_er, 10)

stj_2sw_top_models_er <- get.models(stj_2sw_all_models_er, subset = delta < 2)
stj_2sw_avg_model_er  <- model.avg(stj_2sw_top_models_er, full = TRUE)

summary(stj_2sw_avg_model_er)
stj_2sw_avg_model_er

length(stj_2sw_top_models_er)

# Relative importance of each variable (sum of AICc weights)
sw(stj_2sw_avg_model_er)
confint(stj_2sw_avg_model_er, full = TRUE)


# ── STJ 2SW — Pseudo-R2 ───────────────────────────────────────────────────────

stj_2sw_r2_dev_exp <- bind_rows(
  lapply(stj_2sw_top_models_er, function(x) {
    
    summ_x  <- summary(x)
    
    data.frame(
      model         = paste(as.character(x$formula)[-1], collapse = " ~ "),
      pseudo_R2     = summ_x$pseudo.r.squared
    )
  })
)

mean(stj_2sw_r2_dev_exp$pseudo_R2)
range(stj_2sw_r2_dev_exp$pseudo_R2)


# ── STJ 2SW — Leave-one-out cross-validation (LOO-CV) ─────────────────────────
# With only ~30 observations, LOO-CV is used to
# assess the actual predictive ability of the models (generalisation).
# 1 - SS_res/SS_tot is used for both in-sample AND LOO R2 so
# that the two metrics are directly comparable.

stj_2sw_loo_results_er <- lapply(seq_along(stj_2sw_top_models_er), function(i) {
  
  mod            <- stj_2sw_top_models_er[[i]]
  loo_preds      <- loo_cv(mod, scaled_stj_2sw)
  observed       <- scaled_stj_2sw$return_2sw
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

stj_2sw_loo_summary_er <- bind_rows(stj_2sw_loo_results_er)
stj_2sw_loo_summary_er

mean(stj_2sw_loo_summary_er$rmse_ratio)
range(stj_2sw_loo_summary_er$rmse_ratio)

mean(stj_2sw_loo_summary_er$r2_loo)
range(stj_2sw_loo_summary_er$r2_loo)


stj_2sw_top_models_er[2]
stj_2sw_top_models_er[7]

head(stj_2sw_all_models_er, 10)


# ── STJ 2SW — Model-averaged marginal effects ─────────────────────────────────

all_vars_stj_2sw_er <- unique(unlist(lapply(stj_2sw_top_models_er, function(m) {
  attr(terms(m), "term.labels")
})))

set.seed(42)
avg_plot_data_stj_2sw_er <- lapply(all_vars_stj_2sw_er, function(v) {
  build_marginal_plot_data_avg(
    stj_2sw_top_models_er, v,
    scaled_stj_2sw, stj_2sw, stj_2sw_center, stj_2sw_scale,
    response_var = "return_2sw",
    n_points = 10,
    n_boot   = 20
  )
})

avg_plots_stj_2sw_er <- lapply(
  avg_plot_data_stj_2sw_er,
  function(x) build_marginal_ggplot_avg(x, y_label = "Return rate (2SW)")
)

# Blue curve = weighted mean
# Grey curves = individual models
avg_plots_stj_2sw_er[[1]] +
  labs(title = paste("Variable :", all_vars_stj_2sw_er[1]))

avg_plots_stj_2sw_er[[2]] +
  labs(title = paste("Variable :", all_vars_stj_2sw_er[2]))

avg_plots_stj_2sw_er[[3]] +
  labs(title = paste("Variable :", all_vars_stj_2sw_er[3]))

avg_plots_stj_2sw_er[[4]] +
  labs(title = paste("Variable :", all_vars_stj_2sw_er[4]))

avg_plots_stj_2sw_er[[5]] +
  labs(title = paste("Variable :", all_vars_stj_2sw_er[5]))

avg_plots_stj_2sw_er[[6]] +
  labs(title = paste("Variable :", all_vars_stj_2sw_er[6]))

avg_plots_stj_2sw_er[[7]] +
  labs(title = paste("Variable :", all_vars_stj_2sw_er[7]))


# ══════════════════════════════════════════════════════════════════════════════
# TRINITE RIVER — 2SW
# ══════════════════════════════════════════════════════════════════════════════

# ── TRI 2SW — Linearity check (GAM with smooths) ──────────────────────────────

tri_2sw_gam_smooth <- gam(
  return_2sw ~ s(cohort) +
    s(delta_T_sst_minus_river) +
    s(TriniteCorridor_TH) + s(BelleIsle_TH) + s(LabradorEntrance_TH) + s(LabradorSea_TH) +
    s(Plankton_t) +
    s(smolt_k) +
    s(er_index),
  family = betar(link = "logit"),
  data = scaled_tri_2sw
)
summary(tri_2sw_gam_smooth)

par(mfrow = c(2,2))
plot(tri_2sw_gam_smooth, select = 1, shade = TRUE)
plot(tri_2sw_gam_smooth, select = 2, shade = TRUE)
plot(tri_2sw_gam_smooth, select = 9, shade = TRUE)
par(mfrow = c(1,1))

# Basicaly all linear


# ── TRI 2SW — Fit beta regression ─────────────────────────────────────────────

f_tri_2sw <- as.formula(
  "return_2sw ~ cohort + delta_T_sst_minus_river + Plankton_t +
   TriniteCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
   LabradorSea_TH + smolt_k"
)

tri_2sw_glm <- betareg(
  formula   = f_tri_2sw,
  data      = scaled_tri_2sw,
  link      = "logit",
  na.action = na.fail
)
summary(tri_2sw_glm)


# ── TRI 2SW — Dredge and model averaging ──────────────────────────────────────

tri_2sw_all_models <- dredge(tri_2sw_glm, rank = "AICc")
head(tri_2sw_all_models, 10)

tri_2sw_top_models <- get.models(tri_2sw_all_models, subset = delta < 2)
tri_2sw_avg_model  <- model.avg(tri_2sw_top_models, full = TRUE)

summary(tri_2sw_avg_model)
tri_2sw_avg_model

length(tri_2sw_top_models)

# Relative importance of each variable (sum of AICc weights)
sw(tri_2sw_avg_model)
confint(tri_2sw_avg_model, full = TRUE)


# ── TRI 2SW — Pseudo-R2 ───────────────────────────────────────────────────────

tri_2sw_r2_dev_exp <- bind_rows(
  lapply(tri_2sw_top_models, function(x) {
    
    summ_x  <- summary(x)
    
    data.frame(
      model         = paste(as.character(x$formula)[-1], collapse = " ~ "),
      pseudo_R2     = summ_x$pseudo.r.squared
    )
  })
)

tri_2sw_r2_dev_exp
mean(tri_2sw_r2_dev_exp$pseudo_R2)
range(tri_2sw_r2_dev_exp$pseudo_R2)


# ── TRI 2SW — Leave-one-out cross-validation (LOO-CV) ─────────────────────────

tri_2sw_loo_results <- lapply(seq_along(tri_2sw_top_models), function(i) {
  
  mod            <- tri_2sw_top_models[[i]]
  loo_preds      <- loo_cv(mod, scaled_tri_2sw)
  observed       <- scaled_tri_2sw$return_2sw
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

tri_2sw_loo_summary <- bind_rows(tri_2sw_loo_results)
tri_2sw_loo_summary

mean(tri_2sw_loo_summary$rmse_ratio[-nrow(tri_2sw_loo_summary)])
range(tri_2sw_loo_summary$rmse_ratio[-nrow(tri_2sw_loo_summary)])

mean(tri_2sw_loo_summary$r2_loo[-nrow(tri_2sw_loo_summary)])
range(tri_2sw_loo_summary$r2_loo[-nrow(tri_2sw_loo_summary)])


# ── TRI 2SW — Model-averaged marginal effects ─────────────────────────────────

all_vars_tri_2sw <- unique(unlist(lapply(tri_2sw_top_models, function(m) {
  attr(terms(m), "term.labels")
})))

set.seed(42)
avg_plot_data_tri_2sw <- lapply(all_vars_tri_2sw, function(v) {
  build_marginal_plot_data_avg(
    tri_2sw_top_models, v,
    scaled_tri_2sw, tri_2sw, tri_2sw_center, tri_2sw_scale,
    response_var = "return_2sw",
    n_points = 10,
    n_boot   = 20
  )
})

avg_plots_tri_2sw <- lapply(
  avg_plot_data_tri_2sw,
  function(x) build_marginal_ggplot_avg(x, y_label = "Return rate (2SW)")
)

# Blue curve = weighted mean
# Grey curves = individual models
avg_plots_tri_2sw[[1]] +
  labs(title = paste("Variable :", all_vars_tri_2sw[1]))

avg_plots_tri_2sw[[2]] +
  labs(title = paste("Variable :", all_vars_tri_2sw[2]))

avg_plots_tri_2sw[[3]] +
  labs(title = paste("Variable :", all_vars_tri_2sw[3]))

avg_plots_tri_2sw[[4]] +
  labs(title = paste("Variable :", all_vars_tri_2sw[4]))




