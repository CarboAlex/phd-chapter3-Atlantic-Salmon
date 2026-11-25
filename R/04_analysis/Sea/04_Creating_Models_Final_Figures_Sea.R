# =============================================================================
#
#   Fits the four final global models (same formulas as the marine modelling
#   script), runs dredge, averages the top models (delta AICc < 2), computes
#   model-averaged marginal effects with bootstrap 95% CIs, and exports one
#   combined figure.
#
#   Input : R/04_analysis/Sea/Final_Data/Final_DataSet_{STJ,TRI}_{1SW,2SW}.fst
#   Output: R/04_analysis/Creating_Models_Final_Figures/
#           Combined_MarginalEffects_Sea.png
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
library(patchwork)
library(cowplot)
library(here)


# |=================|
# |      PATHS      |
# |=================|

base_path   <- here::here()
data_path   <- file.path(base_path, "R/04_analysis/Sea/Final_Data")
output_path <- file.path(base_path, "R/04_analysis/Sea")

N_POINTS <- 100
N_BOOT   <- 1000


# |=====================|
# |     IMPORT DATA     |
# |=====================|

stj_1sw <- read_fst(file.path(data_path, "Final_DataSet_STJ_1SW.fst")) %>% dplyr::select(-river)
stj_2sw <- read_fst(file.path(data_path, "Final_DataSet_STJ_2SW.fst")) %>% dplyr::select(-river)
tri_1sw <- read_fst(file.path(data_path, "Final_DataSet_TRI_1SW.fst")) %>% dplyr::select(-river)
tri_2sw <- read_fst(file.path(data_path, "Final_DataSet_TRI_2SW.fst")) %>% dplyr::select(-river)


# |====================|
# |     SCALE DATA     |
# |====================|
# Cohort kept in years; all other covariates standardized (mean 0, SD 1).

scale_dataset <- function(data, response) {
  sp <- scale(data %>% dplyr::select(-cohort, -dplyr::all_of(response)))
  scaled <- cbind(
    data.frame(cohort = data$cohort, setNames(data.frame(data[[response]]), response)),
    as.data.frame(sp)
  )
  list(scaled = scaled,
       center = attr(sp, "scaled:center"),
       scale  = attr(sp, "scaled:scale"))
}

sc_stj_1sw <- scale_dataset(stj_1sw, "return_1sw")
sc_stj_2sw <- scale_dataset(stj_2sw, "return_2sw")
sc_tri_1sw <- scale_dataset(tri_1sw, "return_1sw")
sc_tri_2sw <- scale_dataset(tri_2sw, "return_2sw")


# |========================|
# |     GLOBAL MODELS      |
# |========================|
# Same formulas as the final marine modelling script

f_stj_1sw <- as.formula(
  "return_1sw ~ cohort + delta_T_sst_minus_river + Plankton_t +
   StJeanCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
   LabradorSea_TH + smolt_k + er1NAC_lag1"
)
f_tri_1sw <- as.formula(
  "return_1sw ~ cohort + delta_T_sst_minus_river + Plankton_t +
   TriniteCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
   LabradorSea_TH + smolt_k + er1NAC_lag1"
)
f_stj_2sw <- as.formula(
  "return_2sw ~ cohort + delta_T_sst_minus_river + Plankton_t +
   StJeanCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
   LabradorSea_TH + smolt_k + er_index"
)
f_tri_2sw <- as.formula(
  "return_2sw ~ cohort + delta_T_sst_minus_river + Plankton_t +
   TriniteCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
   LabradorSea_TH + smolt_k"
)

# Fitted at top level (not inside a function): get.models() re-evaluates the
# model call, so the formula and data must be found in the global environment.
scaled_stj_1sw <- sc_stj_1sw$scaled
scaled_stj_2sw <- sc_stj_2sw$scaled
scaled_tri_1sw <- sc_tri_1sw$scaled
scaled_tri_2sw <- sc_tri_2sw$scaled

stj_1sw_glm <- betareg(f_stj_1sw, data = scaled_stj_1sw, link = "logit", na.action = na.fail)
stj_2sw_glm <- betareg(f_stj_2sw, data = scaled_stj_2sw, link = "logit", na.action = na.fail)
tri_1sw_glm <- betareg(f_tri_1sw, data = scaled_tri_1sw, link = "logit", na.action = na.fail)
tri_2sw_glm <- betareg(f_tri_2sw, data = scaled_tri_2sw, link = "logit", na.action = na.fail)


# |=================================|
# |     DREDGE + TOP MODELS         |
# |=================================|

stj_1sw_top_models <- get.models(dredge(stj_1sw_glm, rank = "AICc"), subset = delta < 2)
stj_2sw_top_models <- get.models(dredge(stj_2sw_glm, rank = "AICc"), subset = delta < 2)
tri_1sw_top_models <- get.models(dredge(tri_1sw_glm, rank = "AICc"), subset = delta < 2)
tri_2sw_top_models <- get.models(dredge(tri_2sw_glm, rank = "AICc"), subset = delta < 2)

sapply(list(STJ_1SW = stj_1sw_top_models, STJ_2SW = stj_2sw_top_models,
            TRI_1SW = tri_1sw_top_models, TRI_2SW = tri_2sw_top_models), length)


# =============================================================================
#   FUNCTION — MODEL-AVERAGED MARGINAL EFFECTS (+ BOOTSTRAP CI)
# =============================================================================

build_marginal_plot_data_avg <- function(top_models, var, data_scaled, data_orig,
                                         center, scale_params, response_var,
                                         n_points = 100, n_boot = 500) {

  aicc_vals <- sapply(top_models, AICc)
  weights   <- exp(-0.5 * (aicc_vals - min(aicc_vals)))
  weights   <- weights / sum(weights)

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

compute_all <- function(top_models, vars, sc, data_orig, response_var) {
  set.seed(42)
  out <- lapply(vars, function(v) {
    build_marginal_plot_data_avg(top_models, v, sc$scaled, data_orig, sc$center, sc$scale,
                                 response_var = response_var,
                                 n_points = N_POINTS, n_boot = N_BOOT)
  })
  setNames(out, vars)
}


# =============================================================================
#   PANEL FUNCTIONS
# =============================================================================

# Common colour scale for cohorts across all rivers (one legend for the figure)
all_cohorts   <- c(stj_1sw$cohort, stj_2sw$cohort, tri_1sw$cohort, tri_2sw$cohort)
COHORT_LIMITS <- range(all_cohorts, na.rm = TRUE)
COHORT_MID    <- median(all_cohorts, na.rm = TRUE)

cohort_fill_scale <- function() {
  scale_fill_gradient2(low = "#00897B", mid = "#BDBDBD", high = "#8E24AA",
                       midpoint = COHORT_MID, limits = COHORT_LIMITS, name = "Cohort")
}

build_marginal_ggplot_avg <- function(plot_data, x_label = NULL, y_label = "Return rate (%)",
                                      base_size = 11, point_size = 3, line_width = 1.1,
                                      ribbon_alpha = 0.20) {
  d              <- plot_data
  pred_models_df <- dplyr::bind_rows(d$pred_models)

  p <- ggplot(d$pred_avg, aes(x = x, y = fit))
  if (nrow(pred_models_df) > 0) {
    p <- p + geom_line(data = pred_models_df, aes(x = x, y = fit, group = model),
                       color = "grey55", linewidth = 0.7)
  }
  p +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = ribbon_alpha, fill = "steelblue") +
    geom_line(linewidth = line_width, color = "steelblue4") +
    geom_point(data = d$obs_df, aes(x = x, y = y, fill = cohort), inherit.aes = FALSE,
               size = point_size, color = "black", shape = 21, stroke = 0.45) +
    cohort_fill_scale() +
    scale_y_continuous(labels = function(x) paste0(round(x * 100, 1))) +
    labs(x = if (!is.null(x_label)) x_label else d$var, y = y_label) +
    theme_classic(base_size = base_size) +
    theme(legend.position = "none", axis.title = element_text(size = base_size))
}

# n.r. panel: grey scatter + lm line + r, p (kept for completeness)
make_corr_panel_nr <- function(x, y, x_label = "", y_label = "", base_size = 11) {
  df <- data.frame(x = x, y = y)
  df <- df[complete.cases(df), ]
  ct <- cor.test(df$x, df$y)
  p_label <- if (ct$p.value < 0.001) "p < 0.001" else paste0("p = ", round(ct$p.value, 3))
  ggplot(df, aes(x = x, y = y)) +
    geom_smooth(method = "lm", se = TRUE, color = "grey50", fill = "grey80",
                linewidth = 0.8, alpha = 0.35) +
    geom_point(color = "grey55", fill = "grey70", shape = 21, size = 2.5, stroke = 0.4) +
    annotate("text", x = -Inf, y = Inf,
             label = paste0("n.r.\nr = ", round(ct$estimate, 2), "\n", p_label),
             hjust = -0.10, vjust = 1.3, size = (base_size - 3) / ggplot2::.pt,
             color = "grey45", fontface = "italic", lineheight = 1.3) +
    scale_y_continuous(labels = function(x) paste0(round(x * 100, 1))) +
    labs(x = x_label, y = y_label) +
    theme_classic(base_size = base_size) +
    theme(axis.line  = element_line(color = "grey65"),
          axis.ticks = element_line(color = "grey65"),
          axis.text  = element_text(color = "grey55"),
          axis.title = element_text(color = "grey50", size = base_size - 2))
}

add_status_mark <- function(p, status, base_size = 11) {
  if (status == "ns") {
    p <- p + annotate("text", x = -Inf, y = Inf, label = "n.s.", hjust = -0.25, vjust = 1.5,
                      size = (base_size + 3) / ggplot2::.pt, fontface = "italic")
  }
  if (status == "cond_only") {
    p <- p + annotate("text", x = -Inf, y = Inf, label = "*", hjust = -0.5, vjust = 1.1,
                      size = 9, fontface = "bold")
  }
  p
}


# =============================================================================
#   FIGURE BUILDERS
# =============================================================================

# One river block: grey category header row + one row of panels per sea age
build_river_block <- function(row_specs, col_specs, x_labels,
                              base_size = 11, strip_size = 11) {

  row_panels <- lapply(row_specs, function(r) {
    plots_in_row <- lapply(col_specs, function(col) {
      v      <- col$var
      status <- r$var_status[[v]]
      x_lab  <- x_labels[[v]]
      if (status == "nr") {
        return(make_corr_panel_nr(r$data_orig[[v]], r$data_orig[[r$response_col]],
                                  x_label = x_lab, y_label = r$y_label, base_size = base_size))
      }
      p <- build_marginal_ggplot_avg(r$plot_data[[v]], x_label = x_lab, y_label = r$y_label,
                                     base_size = base_size)
      add_status_mark(p, status, base_size)
    })
    title_grob <- cowplot::ggdraw() +
      cowplot::draw_label(r$label, x = 0, hjust = 0, size = strip_size, fontface = "bold")
    cowplot::plot_grid(title_grob,
                       patchwork::wrap_plots(plots_in_row, nrow = 1),
                       ncol = 1, rel_heights = c(0.12, 1))
  })

  cat_panels <- lapply(col_specs, function(col) {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = col$category,
               size = (base_size - 2) / ggplot2::.pt) +
      theme_void() +
      theme(plot.background = element_rect(fill = "grey85", color = "grey50", linewidth = 0.4))
  })

  cowplot::plot_grid(
    patchwork::wrap_plots(cat_panels, nrow = 1),
    cowplot::plot_grid(plotlist = row_panels, ncol = 1),
    ncol = 1, rel_heights = c(0.10, 1)
  )
}

export_combined_figure <- function(stj_row_specs, stj_col_specs,
                                   tri_row_specs, tri_col_specs,
                                   x_labels, filename,
                                   base_size = 11, dpi = 300,
                                   width_mm = 310, height_mm = 300) {

  stj_block <- build_river_block(stj_row_specs, stj_col_specs, x_labels, base_size)
  tri_block <- build_river_block(tri_row_specs, tri_col_specs, x_labels, base_size)

  # Common cohort legend (vertical colourbar, recent cohorts on top)
  p_leg <- build_marginal_ggplot_avg(stj_row_specs[[1]]$plot_data[[1]], base_size = base_size) +
    guides(fill = guide_colorbar(title = "Cohort", barheight = unit(8, "cm"),
                                 barwidth = unit(0.45, "cm"), direction = "vertical",
                                 title.position = "top", reverse = TRUE)) +
    theme(legend.position = "right",
          legend.text     = element_text(size = base_size - 2),
          legend.title    = element_text(size = base_size - 1, face = "bold"))
  cohort_legend <- cowplot::get_legend(p_leg)

  final_figure <- cowplot::plot_grid(
    cowplot::plot_grid(stj_block, tri_block, ncol = 1),
    cohort_legend,
    nrow = 1, rel_widths = c(1, 0.07)
  )

  ggplot2::ggsave(filename = filename, plot = final_figure, device = "png", dpi = dpi,
                  width = width_mm, height = height_mm, units = "mm", bg = "white")
  message("Exported : ", basename(filename))
  invisible(final_figure)
}


# =============================================================================
#   PANELS: VARIABLES, ORDER AND STATUS
# =============================================================================

x_labels <- list(
  cohort         = "Cohort",
  BelleIsle_TH   = "Thermal habitat\n(Strait of Belle Isle)",
  LabradorSea_TH = "Thermal habitat\n(Labrador Sea)",
  smolt_k        = "Smolt condition (K)",
  Plankton_t     = "Food availability (g)\n(Labrador Sea)"
)

# ── St. Jean: cohort | SOBI TH | smolt condition | food ──────────────────────
stj_col_specs <- list(
  list(var = "cohort",       category = "Temporal trend"),
  list(var = "BelleIsle_TH", category = "Marine thermal habitat"),
  list(var = "smolt_k",      category = "Smolt condition"),
  list(var = "Plankton_t",   category = "Food availability")
)
stj_vars <- sapply(stj_col_specs, `[[`, "var")

# ── Trinité: SOBI TH | Labrador Sea TH | smolt condition | food ──────────────
tri_col_specs <- list(
  list(var = "BelleIsle_TH",   category = "Marine thermal habitat"),
  list(var = "LabradorSea_TH", category = "Marine thermal habitat"),
  list(var = "smolt_k",        category = "Smolt condition"),
  list(var = "Plankton_t",     category = "Food availability")
)
tri_vars <- sapply(tri_col_specs, `[[`, "var")


# =============================================================================
#   MARGINAL EFFECTS
# =============================================================================

avg_plot_data_stj_1sw <- compute_all(stj_1sw_top_models, stj_vars, sc_stj_1sw, stj_1sw, "return_1sw")
avg_plot_data_stj_2sw <- compute_all(stj_2sw_top_models, stj_vars, sc_stj_2sw, stj_2sw, "return_2sw")
avg_plot_data_tri_1sw <- compute_all(tri_1sw_top_models, tri_vars, sc_tri_1sw, tri_1sw, "return_1sw")
avg_plot_data_tri_2sw <- compute_all(tri_2sw_top_models, tri_vars, sc_tri_2sw, tri_2sw, "return_2sw")


# =============================================================================
#   FIGURE
# =============================================================================

stj_row_specs <- list(
  list(label = "St. Jean 1SW", plot_data = avg_plot_data_stj_1sw, data_orig = stj_1sw,
       response_col = "return_1sw", y_label = "1SW return rate (%)",
       var_status = list(cohort = "sig", BelleIsle_TH = "sig", smolt_k = "sig", Plankton_t = "sig")),
  list(label = "St. Jean 2SW", plot_data = avg_plot_data_stj_2sw, data_orig = stj_2sw,
       response_col = "return_2sw", y_label = "2SW return rate (%)",
       var_status = list(cohort = "sig", BelleIsle_TH = "ns", smolt_k = "cond_only", Plankton_t = "sig"))
)

tri_row_specs <- list(
  list(label = "Trinité 1SW", plot_data = avg_plot_data_tri_1sw, data_orig = tri_1sw,
       response_col = "return_1sw", y_label = "1SW return rate (%)",
       var_status = list(BelleIsle_TH = "ns", LabradorSea_TH = "sig", smolt_k = "cond_only",
                         Plankton_t = "cond_only")),
  list(label = "Trinité 2SW", plot_data = avg_plot_data_tri_2sw, data_orig = tri_2sw,
       response_col = "return_2sw", y_label = "2SW return rate (%)",
       var_status = list(BelleIsle_TH = "sig", LabradorSea_TH = "ns", smolt_k = "ns", Plankton_t = "sig"))
)

export_combined_figure(
  stj_row_specs = stj_row_specs, stj_col_specs = stj_col_specs,
  tri_row_specs = tri_row_specs, tri_col_specs = tri_col_specs,
  x_labels      = x_labels,
  filename      = file.path(output_path, "Combined_MarginalEffects_Sea.png"),
  width_mm      = 310,
  height_mm     = 300
)
