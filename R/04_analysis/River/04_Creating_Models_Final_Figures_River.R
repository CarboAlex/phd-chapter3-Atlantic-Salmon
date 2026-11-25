# =============================================================================
#
#   Fits the final St. Jean global model (same covariates as 03_River_Modeling.R),
#   runs dredge, averages the top models (delta AICc < 2), computes
#   model-averaged marginal effects with bootstrap 95% CIs.
#
#   Trinité: no figure (null model in the top set, no covariate significant
#   in the model average; see 03_River_Modeling.R).
#
#   Input : R/04_analysis/River/Final_Data/Final_DataSet_STJ.fst
#   Output: R/04_analysis/River/Combined_MarginalEffects_River.png
#
# =============================================================================


# |========================|
# |     LOAD LIBRARIES     |
# |========================|

library(dplyr)
library(fst)
library(ggplot2)
library(MuMIn)
library(betareg)
library(cowplot)
library(here)


# |=================|
# |      PATHS      |
# |=================|

base_path   <- here::here()
output_path <- file.path(base_path, "R/04_analysis/River")


# |=====================|
# |     IMPORT DATA     |
# |=====================|

stj <- read_fst(file.path(base_path, "R/04_analysis/River/Final_Data/Final_DataSet_STJ.fst"))


# |====================|
# |     SCALE DATA     |
# |====================|

stj_scale_params <- scale(stj %>% dplyr::select(-cohort, -surv_egg_smolt))
stj_center       <- attr(stj_scale_params, "scaled:center")
stj_scale        <- attr(stj_scale_params, "scaled:scale")

scaled_stj <- cbind(
  data.frame(cohort = stj$cohort, surv_egg_smolt = stj$surv_egg_smolt),
  as.data.frame(stj_scale_params)
)


# |========================|
# |     GLOBAL MODEL       |
# |========================|
# Same covariates as the final St. Jean global model in 03_River_Modeling.R

stj_vars <- c(
  "cohort",
  "co_gr_median_temp", "em_median_temp", "em_var_temp",
  "co_gr_dd_o_20_temp", "co_year_mean_sum_ptgi",
  "co_gr_mean_disc",
  "em_quant_75qt_disc",
  "co_mean_forest"
)
stj_formula <- as.formula(paste("surv_egg_smolt ~", paste(stj_vars, collapse = " + ")))

stj_glm <- betareg(
  formula   = stj_formula,
  data      = scaled_stj,
  link      = "logit",
  na.action = na.fail
)


# |=========================|
# |     DREDGE + TOP M.     |
# |=========================|

stj_all_models <- dredge(stj_glm, rank = "AICc")
stj_top_models <- get.models(stj_all_models, subset = delta < 2)
stj_avg_model  <- model.avg(stj_top_models)

length(stj_top_models)
summary(stj_avg_model)


# =============================================================================
#   FUNCTION — MODEL-AVERAGED MARGINAL EFFECTS (+ BOOTSTRAP CI)
# =============================================================================

build_marginal_plot_data_avg <- function(top_models, var, data_scaled, data_orig,
                                         center, scale_params,
                                         response_var = "surv_egg_smolt",
                                         n_points = 100, n_boot = 500) {

  aicc_vals <- sapply(top_models, AICc)
  delta     <- aicc_vals - min(aicc_vals)
  weights   <- exp(-0.5 * delta) / sum(exp(-0.5 * delta))

  var_seq <- seq(min(data_scaled[[var]], na.rm = TRUE),
                 max(data_scaled[[var]], na.rm = TRUE),
                 length.out = n_points)

  nd <- data_scaled[rep(1, n_points), !names(data_scaled) %in% response_var, drop = FALSE]
  for (col in names(nd)) nd[[col]] <- mean(data_scaled[[col]], na.rm = TRUE)
  nd[[var]] <- var_seq

  preds_by_model <- lapply(top_models, function(m) as.numeric(predict(m, newdata = nd, type = "response")))
  fit_avg <- Reduce("+", mapply(function(p, w) p * w, preds_by_model, weights, SIMPLIFY = FALSE))

  # Bootstrap: resample cohorts, refit every top model, recompute AICc weights
  boot_mat <- matrix(NA, nrow = n_points, ncol = n_boot)
  for (b in seq_len(n_boot)) {
    idx    <- sample(nrow(data_scaled), replace = TRUE)
    mods_b <- lapply(top_models, function(m) {
      tryCatch(betareg(formula(m), data = data_scaled[idx, ], link = "logit"),
               error = function(e) NULL)
    })
    mods_b <- mods_b[!sapply(mods_b, is.null)]
    if (length(mods_b) == 0) next

    aicc_b  <- sapply(mods_b, AICc)
    wts_b   <- exp(-0.5 * (aicc_b - min(aicc_b)))
    wts_b   <- wts_b / sum(wts_b)
    preds_b <- lapply(mods_b, function(m) as.numeric(predict(m, newdata = nd, type = "response")))
    boot_mat[, b] <- Reduce("+", mapply(function(p, w) p * w, preds_b, wts_b, SIMPLIFY = FALSE))
  }
  lower <- apply(boot_mat, 1, quantile, 0.025, na.rm = TRUE)
  upper <- apply(boot_mat, 1, quantile, 0.975, na.rm = TRUE)

  x_orig <- if (var == "cohort") var_seq else var_seq * scale_params[var] + center[var]

  contains_var <- sapply(top_models, function(m) var %in% attr(terms(m), "term.labels"))

  list(
    var      = var,
    pred_avg = data.frame(x = x_orig, fit = as.numeric(fit_avg),
                          lower = as.numeric(lower), upper = as.numeric(upper)),
    pred_models = lapply(which(contains_var), function(i) {
      data.frame(x = x_orig, fit = preds_by_model[[i]],
                 model = paste0("M", i, " (w=", round(weights[i], 2), ")"))
    }),
    obs_df = data.frame(x = data_orig[[var]], y = data_orig[[response_var]],
                        cohort = data_orig$cohort)
  )
}


# =============================================================================
#   FUNCTION — ONE PANEL
# =============================================================================

build_marginal_ggplot_avg <- function(plot_data,
                                      x_label      = NULL,
                                      y_label      = "Egg-to-smolt survival (%)",
                                      conditional  = FALSE,
                                      base_size    = 13,
                                      text_scale   = 1.4,
                                      point_size   = 3.5,
                                      line_width   = 1.2,
                                      ribbon_alpha = 0.20) {

  bs             <- base_size * text_scale
  d              <- plot_data
  cohort_mid     <- median(d$obs_df$cohort, na.rm = TRUE)
  pred_models_df <- dplyr::bind_rows(d$pred_models)

  p <- ggplot(d$pred_avg, aes(x = x, y = fit))

  if (nrow(pred_models_df) > 0) {
    p <- p + geom_line(data = pred_models_df, aes(x = x, y = fit, group = model),
                       color = "grey55", linewidth = 0.7)
  }

  p <- p +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = ribbon_alpha, fill = "steelblue") +
    geom_line(linewidth = line_width, color = "steelblue4") +
    geom_point(data = d$obs_df, aes(x = x, y = y, fill = cohort), inherit.aes = FALSE,
               size = point_size, color = "black", shape = 21, stroke = 0.45) +
    scale_fill_gradient2(low = "#00897B", mid = "#BDBDBD", high = "#8E24AA",
                         midpoint = cohort_mid, name = "Cohort") +
    scale_y_continuous(labels = function(x) paste0(round(x * 100, 1))) +
    labs(x = if (!is.null(x_label)) x_label else d$var, y = y_label) +
    theme_classic(base_size = bs) +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.45 * text_scale, "cm"),
      legend.text     = element_text(size = 9  * text_scale),
      legend.title    = element_text(size = 10 * text_scale),
      axis.title      = element_text(size = bs)
    )

  # Asterisk = significant in the conditional average only (top-left corner)
  if (conditional) {
    p <- p + annotate("text", x = -Inf, y = Inf, label = "*",
                      hjust = -0.6, vjust = 1.1, size = 9, fontface = "bold")
  }
  p
}


# =============================================================================
#   FIGURE DATA — ST. JEAN
# =============================================================================

fig_vars <- c("cohort", "co_gr_median_temp", "co_year_mean_sum_ptgi")

N_POINTS <- 100
N_BOOT   <- 1000

set.seed(42)
avg_plot_data_stj <- lapply(fig_vars, function(v) {
  build_marginal_plot_data_avg(
    stj_top_models, v,
    scaled_stj, stj, stj_center, stj_scale,
    response_var = "surv_egg_smolt",
    n_points = N_POINTS, n_boot = N_BOOT
  )
})
names(avg_plot_data_stj) <- fig_vars


# =============================================================================
#   FIGURE 
# =============================================================================

export_stj_figure <- function(plot_data_list,
                              col_specs,
                              filename,
                              base_size    = 13,
                              text_scale   = 1.4,
                              point_size   = 3,
                              line_width   = 1.1,
                              ribbon_alpha = 0.20,
                              dpi          = 300,
                              width_mm     = 260,
                              height_mm    = 230) {

  make_panel <- function(spec) {
    build_marginal_ggplot_avg(
      plot_data_list[[spec$var]],
      x_label      = spec$x_label,
      conditional  = spec$conditional,
      base_size    = base_size,
      text_scale   = text_scale,
      point_size   = point_size,
      line_width   = line_width,
      ribbon_alpha = ribbon_alpha
    )
  }

  p_cohort <- make_panel(col_specs$cohort)            + theme(legend.position = "none")
  p_temp   <- make_panel(col_specs$co_gr_median_temp) + theme(legend.position = "none")
  p_pgti   <- make_panel(col_specs$co_year_mean_sum_ptgi) + theme(legend.position = "none")

  # Cohort colour legend (vertical colourbar, recent cohorts on top)
  p_leg <- make_panel(col_specs$cohort) +
    guides(fill = guide_colorbar(
      title          = "Cohort",
      barheight      = unit(6 * text_scale, "cm"),
      barwidth       = unit(0.5 * text_scale, "cm"),
      direction      = "vertical",
      title.position = "top",
      reverse        = TRUE
    )) +
    theme(legend.position = "right",
          legend.text     = element_text(size = 9  * text_scale),
          legend.title    = element_text(size = 10 * text_scale, face = "bold"))
  cohort_legend <- cowplot::get_legend(p_leg)

  # Align the three panels so their plotting areas line up
  aligned <- cowplot::align_plots(p_cohort, p_temp, p_pgti, align = "hv", axis = "tblr")

  final_figure <- cowplot::plot_grid(
    aligned[[1]], cohort_legend,
    aligned[[2]], aligned[[3]],
    ncol = 2, nrow = 2
  )

  ggplot2::ggsave(filename = filename, plot = final_figure, device = "png", dpi = dpi,
                  width = width_mm, height = height_mm, units = "mm", bg = "white")
  message("Exported : ", basename(filename))
  invisible(final_figure)
}


# ── Panel specifications ──────────────────────────────────────────────────────

col_specs_stj <- list(
  cohort = list(
    var         = "cohort",
    x_label     = "Cohort",
    conditional = FALSE
  ),
  co_gr_median_temp = list(
    var         = "co_gr_median_temp",
    x_label     = "River median temperature\n(growth period, °C)",
    conditional = TRUE
  ),
  co_year_mean_sum_ptgi = list(
    var         = "co_year_mean_sum_ptgi",
    x_label     = "Potential growth thermal index\n(growth period)",
    conditional = TRUE
  )
)


# ── Export ────────────────────────────────────────────────────────────────────

export_stj_figure(
  plot_data_list = avg_plot_data_stj,
  col_specs      = col_specs_stj,
  filename       = file.path(output_path, "Combined_MarginalEffects_River.png")
)


