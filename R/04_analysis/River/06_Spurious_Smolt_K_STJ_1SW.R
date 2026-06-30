

# =============================================================================
#
#   STJ 1SW — RELATION POTENTIELLEMENT SPURIEUSE DE smolt_k
#
#   Contexte : la relation négative smolt_k ~ return_1sw observée dans M1
#   est biologiquement inattendue (un meilleur coefficient de Fulton devrait
#   améliorer la survie, non la diminuer). On teste ici si cette relation
#   reflète un effet causal ou une covariation temporelle commune.
#
# =============================================================================


# |=================|
# |    BASE PATH    |
# |=================|

base_path <- "C:/Users/carbo/Documents/Alex/Ecole/Doctorat/Chapitre_3"


# |========================|
# |     LOAD LIBRARIES     |
# |========================|

library(dplyr)
library(fst)
library(ggplot2)
library(patchwork)
library(betareg)
library(MuMIn)
library(ggrepel)
library(ppcor)


# |=====================|
# |     IMPORT DATA     |
# |=====================|

stj_1sw <- read_fst(file.path(base_path, "04_Analysis/Sea/Final_DataSet_STJ_1SW.fst")) %>%
  dplyr::select(-river)


# |====================|
# |     SCALE DATA     |
# |====================|

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


# |======================|
# |     MODÈLES M1-VD    |
# |======================|

# M1-VD : modèle retenu (beta regression, phi ~ cohort)
stj_1sw_m1_vd <- betareg(
  formula = return_1sw ~ BelleIsle_TH + cohort + plankton_index + smolt_k | cohort,
  data    = scaled_stj_1sw,
  link    = "logit"
)

# Top models du dredge (nécessaires pour build_marginal_plot_data_avg)
stj_1sw_m1 <- betareg(
  formula   = return_1sw ~ BelleIsle_TH + cohort + plankton_index + smolt_k,
  data      = scaled_stj_1sw,
  link      = "logit",
  na.action = na.fail
)
stj_1sw_all_models <- dredge(stj_1sw_m1, rank = "AICc")
stj_1sw_top_models <- get.models(stj_1sw_all_models, subset = delta < 2)
stj_1sw_avg_model  <- model.avg(stj_1sw_top_models, full = TRUE)


# |==============================|
# |     FONCTIONS MARGINALES     |
# |==============================|

build_marginal_plot_data <- function(model, var, data_scaled, data_orig,
                                     center, scale_params,
                                     response_var = "return_1sw",
                                     n_points = 10, n_boot = 20) {
  
  var_seq <- seq(min(data_scaled[[var]], na.rm = TRUE),
                 max(data_scaled[[var]], na.rm = TRUE),
                 length.out = n_points)
  
  nd <- data_scaled[rep(1, n_points),
                    !names(data_scaled) %in% response_var, drop = FALSE]
  for (col in names(nd)) nd[[col]] <- mean(data_scaled[[col]], na.rm = TRUE)
  nd[[var]] <- var_seq
  
  fit <- as.numeric(predict(model, newdata = nd, type = "response"))
  
  boot_mat <- matrix(NA, nrow = n_points, ncol = n_boot)
  for (b in seq_len(n_boot)) {
    idx   <- sample(nrow(data_scaled), replace = TRUE)
    mod_b <- tryCatch(
      betareg(formula(model), data = data_scaled[idx, ], link = "logit"),
      error = function(e) NULL
    )
    if (!is.null(mod_b))
      boot_mat[, b] <- predict(mod_b, newdata = nd, type = "response")
  }
  lower <- apply(boot_mat, 1, quantile, 0.025, na.rm = TRUE)
  upper <- apply(boot_mat, 1, quantile, 0.975, na.rm = TRUE)
  
  x_orig <- if (var == "cohort") var_seq else
    var_seq * scale_params[var] + center[var]
  
  list(
    var     = var,
    pred_df = data.frame(x = x_orig, fit = fit,
                         lower = as.numeric(lower),
                         upper = as.numeric(upper)),
    obs_df  = data.frame(x      = data_orig[[var]],
                         y      = data_orig[[response_var]],
                         cohort = data_orig$cohort)
  )
}

build_marginal_ggplot <- function(plot_data, base_size = 13, point_size = 3.5,
                                  label_size = 3, line_width = 1.4,
                                  ribbon_alpha = 0.20,
                                  y_label = "1SW return rate") {
  d          <- plot_data
  cohort_mid <- median(d$obs_df$cohort, na.rm = TRUE)
  
  ggplot(d$pred_df, aes(x, fit)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = ribbon_alpha, fill = "steelblue") +
    geom_line(linewidth = line_width, colour = "steelblue4") +
    geom_point(data = d$obs_df, aes(x, y, fill = cohort),
               inherit.aes = FALSE,
               size = point_size, colour = "black", shape = 21, stroke = 0.5) +
    scale_fill_gradient2(low = "#00897B", mid = "#BDBDBD", high = "#8E24AA",
                         midpoint = cohort_mid, name = "Cohort") +
    geom_text_repel(data = d$obs_df, aes(x, y, label = cohort),
                    inherit.aes = FALSE,
                    size = label_size, max.overlaps = 20) +
    labs(x = d$var, y = y_label) +
    theme_classic(base_size = base_size)
}

build_marginal_plot_data_avg <- function(top_models, var, data_scaled, data_orig,
                                         center, scale_params,
                                         response_var = "return_1sw",
                                         n_points = 10, n_boot = 20) {
  
  aicc_vals <- sapply(top_models, AICc)
  delta_w   <- aicc_vals - min(aicc_vals)
  weights   <- exp(-0.5 * delta_w) / sum(exp(-0.5 * delta_w))
  
  var_seq <- seq(min(data_scaled[[var]], na.rm = TRUE),
                 max(data_scaled[[var]], na.rm = TRUE),
                 length.out = n_points)
  
  nd <- data_scaled[rep(1, n_points),
                    !names(data_scaled) %in% response_var, drop = FALSE]
  for (col in names(nd)) nd[[col]] <- mean(data_scaled[[col]], na.rm = TRUE)
  nd[[var]] <- var_seq
  
  preds_by_model <- lapply(top_models, function(m)
    as.numeric(predict(m, newdata = nd, type = "response")))
  
  fit_avg <- Reduce("+", mapply(function(p, w) p * w,
                                preds_by_model, weights, SIMPLIFY = FALSE))
  
  boot_mat <- matrix(NA, nrow = n_points, ncol = n_boot)
  for (b in seq_len(n_boot)) {
    idx    <- sample(nrow(data_scaled), replace = TRUE)
    mods_b <- lapply(top_models, function(m)
      tryCatch(betareg(formula(m), data = data_scaled[idx, ], link = "logit"),
               error = function(e) NULL))
    valid  <- !sapply(mods_b, is.null)
    mods_b <- mods_b[valid]
    if (!length(mods_b)) next
    
    aicc_b <- sapply(mods_b, AICc)
    wts_b  <- exp(-0.5 * (aicc_b - min(aicc_b)))
    wts_b  <- wts_b / sum(wts_b)
    preds_b <- lapply(mods_b, function(m)
      as.numeric(predict(m, newdata = nd, type = "response")))
    boot_mat[, b] <- Reduce("+", mapply(function(p, w) p * w,
                                        preds_b, wts_b, SIMPLIFY = FALSE))
  }
  lower <- apply(boot_mat, 1, quantile, 0.025, na.rm = TRUE)
  upper <- apply(boot_mat, 1, quantile, 0.975, na.rm = TRUE)
  
  x_orig       <- if (var == "cohort") var_seq else
    var_seq * scale_params[var] + center[var]
  contains_var <- sapply(top_models, function(m) var %in% attr(terms(m), "term.labels"))
  
  list(
    var      = var,
    pred_avg = data.frame(x = x_orig, fit = as.numeric(fit_avg),
                          lower = as.numeric(lower),
                          upper = as.numeric(upper)),
    pred_models = lapply(which(contains_var), function(i) {
      data.frame(
        x     = x_orig,
        fit   = preds_by_model[[i]],
        model = paste0("M", i, " (w=", round(weights[i], 2), "): ",
                       paste(attr(terms(top_models[[i]]), "term.labels"),
                             collapse = " + "))
      )
    }),
    obs_df = data.frame(x      = data_orig[[var]],
                        y      = data_orig[[response_var]],
                        cohort = data_orig$cohort)
  )
}

build_marginal_ggplot_avg <- function(plot_data, base_size = 13, point_size = 3.5,
                                      label_size = 3, line_width = 1.4,
                                      ribbon_alpha = 0.20,
                                      y_label = "1SW return rate") {
  d             <- plot_data
  cohort_mid    <- median(d$obs_df$cohort, na.rm = TRUE)
  pred_models_df <- dplyr::bind_rows(d$pred_models)
  
  p <- ggplot(d$pred_avg, aes(x, fit))
  if (nrow(pred_models_df) > 0)
    p <- p + geom_line(data = pred_models_df,
                       aes(x, fit, group = model, linetype = model),
                       colour = "grey50", linewidth = 0.8)
  
  p +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = ribbon_alpha, fill = "steelblue") +
    geom_line(linewidth = line_width, colour = "steelblue4") +
    geom_point(data = d$obs_df, aes(x, y, fill = cohort),
               inherit.aes = FALSE,
               size = point_size, colour = "black", shape = 21, stroke = 0.5) +
    scale_fill_gradient2(low = "#00897B", mid = "#BDBDBD", high = "#8E24AA",
                         midpoint = cohort_mid, name = "Cohort") +
    geom_text_repel(data = d$obs_df, aes(x, y, label = cohort),
                    inherit.aes = FALSE,
                    size = label_size, max.overlaps = 20) +
    labs(x = d$var, y = y_label) +
    theme_classic(base_size = base_size)
}


# =============================================================================
#
#   TEST 1 — TENDANCE TEMPORELLE DE smolt_k
#
#   Si smolt_k a décliné dans le temps au même rythme que return_1sw,
#   la relation négative dans le modèle pourrait être spurieuse —
#   les deux variables déclinent ensemble sans lien causal direct.
#
# =============================================================================

cor.test(stj_1sw$cohort, stj_1sw$smolt_k)

# Visualisation de la tendance
ggplot(stj_1sw, aes(x = cohort, y = smolt_k)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, colour = "steelblue4") +
  geom_text_repel(aes(label = cohort), size = 3) +
  labs(
    title    = "Tendance temporelle de smolt_k — STJ 1SW",
    subtitle = paste0(
      "r = ", round(cor(stj_1sw$cohort, stj_1sw$smolt_k, use = "complete.obs"), 2)
    ),
    x = "Cohorte",
    y = "Coefficient de Fulton (k)"
  ) +
  theme_classic(base_size = 13)

# Attendu : r ≈ -0.58, p < 0.001 -> smolt_k décline significativement
# avec le temps -> risque de covariation temporelle spurieuse avec return_1sw


# =============================================================================
#
#   TEST 2 — CORRÉLATION PARTIELLE (smolt_k ~ return_1sw | cohort)
#
#   On retire l'effet du temps des deux variables et on corrèle les résidus.
#   Si la corrélation résiduelle est non significative, la relation entre
#   smolt_k et return_1sw est spurieuse — elle provient du temps, pas d'un
#   effet causal de la condition des smolts sur la survie en mer.
#
# =============================================================================

df_pcor <- stj_1sw %>%
  dplyr::select(cohort, return_1sw, smolt_k) %>%
  na.omit()

pcor.test(df_pcor$return_1sw, df_pcor$smolt_k, df_pcor$cohort)

# Visualisation côte à côte : corrélation brute vs corrélation partielle
res_smolt_k  <- residuals(lm(smolt_k   ~ cohort, data = df_pcor))
res_return   <- residuals(lm(return_1sw ~ cohort, data = df_pcor))

cohort_mid <- median(stj_1sw$cohort, na.rm = TRUE)

p_brute <- ggplot(df_pcor, aes(x = smolt_k, y = return_1sw)) +
  geom_smooth(method = "lm", se = FALSE, colour = "steelblue4", linewidth = 1.2) +
  geom_point(aes(fill = cohort), shape = 21, size = 3.5, stroke = 0.5) +
  scale_fill_gradient2(low = "#00897B", mid = "#BDBDBD", high = "#8E24AA",
                       midpoint = cohort_mid, name = "Cohort") +
  geom_text_repel(aes(label = cohort), size = 3, max.overlaps = 20) +
  labs(
    title    = "Corrélation brute",
    subtitle = paste0(
      "r = ", round(cor(df_pcor$smolt_k, df_pcor$return_1sw), 2),
      "  p = ", round(cor.test(df_pcor$smolt_k, df_pcor$return_1sw)$p.value, 3)
    ),
    x = "smolt_k",
    y = "return_1sw"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.subtitle = element_text(colour = "steelblue4"))

p_partielle <- ggplot(data.frame(res_smolt_k, res_return, cohort = df_pcor$cohort),
                      aes(x = res_smolt_k, y = res_return)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_smooth(method = "lm", se = FALSE, colour = "grey50",
              linewidth = 1.2, linetype = "dashed") +
  geom_point(aes(fill = cohort), shape = 21, size = 3.5, stroke = 0.5) +
  scale_fill_gradient2(low = "#00897B", mid = "#BDBDBD", high = "#8E24AA",
                       midpoint = cohort_mid, name = "Cohort") +
  geom_text_repel(aes(label = cohort), size = 3, max.overlaps = 20) +
  labs(
    title    = "Corrélation partielle (| cohort)",
    subtitle = paste0(
      "r = ", round(pcor.test(df_pcor$return_1sw, df_pcor$smolt_k, df_pcor$cohort)$estimate, 2),
      "  p = ", round(pcor.test(df_pcor$return_1sw, df_pcor$smolt_k, df_pcor$cohort)$p.value, 3)
    ),
    x = "Résidus smolt_k | cohort",
    y = "Résidus return_1sw | cohort"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.subtitle = element_text(colour = "grey40"))

p_brute + p_partielle + plot_layout(guides = "collect")

# Attendu : r partiel ≈ -0.17, p ≈ 0.38 -> corrélation non significative
# une fois la tendance temporelle contrôlée -> relation brute spurieuse


# =============================================================================
#
#   TEST 3A — COMPARAISON DES MODÈLES : AVEC VS SANS smolt_k
#
#   Si smolt_k masque l'effet d'autres variables, leurs coefficients
#   devraient augmenter substantiellement à son retrait. Si au contraire
#   smolt_k amplifie l'effet apparent d'autres variables, leurs coefficients
#   devraient diminuer (et potentiellement perdre leur significativité).
#
# =============================================================================

# Modèle sans smolt_k — même structure que M1-VD pour comparaison équitable
m_sans_k_vd <- betareg(
  formula = return_1sw ~ BelleIsle_TH + cohort + plankton_index | cohort,
  data    = scaled_stj_1sw,
  link    = "logit"
)

vars_interet <- c("BelleIsle_TH", "cohort", "plankton_index")

coef_avec_k <- coef(stj_1sw_m1_vd, model = "mean")
coef_sans_k <- coef(m_sans_k_vd,   model = "mean")

cat("=== Comparaison des coefficients : avec vs sans smolt_k ===\n")
print(
  data.frame(
    variable    = vars_interet,
    coef_avec_k = round(coef_avec_k[vars_interet], 4),
    coef_sans_k = round(coef_sans_k[vars_interet], 4),
    diff_pct    = round(
      abs(coef_avec_k[vars_interet] - coef_sans_k[vars_interet]) /
        abs(coef_avec_k[vars_interet]) * 100, 1
    )
  )
)

# Significativité dans les deux modèles
summary(stj_1sw_m1_vd)   # avec smolt_k
summary(m_sans_k_vd)     # sans smolt_k


# =============================================================================
#
#   TEST 3B — SÉLECTION DE VARIABLES SANS smolt_k
#
#   Si la structure de sélection devient chaotique sans smolt_k
#   (beaucoup de modèles compétiteurs avec des poids faibles),
#   c'est un signe que smolt_k contribuait à structurer l'information
#   dans les données malgré sa relation potentiellement spurieuse.
#
# =============================================================================

f_nok <- as.formula(
  "return_1sw ~ cohort + delta_T_sst_minus_river + plankton_index +
   StJeanCorridor_TH + BelleIsle_TH + LabradorEntrance_TH +
   LabradorSea_TH + er1NAC_lag1 + pfa_1sw"
)

stj_1sw_nok_glm <- betareg(
  formula   = f_nok,
  data      = scaled_stj_1sw,
  link      = "logit",
  na.action = na.fail
)

stj_1sw_nok_all_models <- dredge(stj_1sw_nok_glm, rank = "AICc")
stj_1sw_nok_top_models <- get.models(stj_1sw_nok_all_models, subset = delta < 2)
stj_1sw_nok_avg_model  <- model.avg(stj_1sw_nok_top_models, full = TRUE)

cat("=== Nombre de modèles retenus sans smolt_k :", length(stj_1sw_nok_top_models), "===\n")
cat("=== Nombre de modèles retenus avec smolt_k :", length(stj_1sw_top_models),     "===\n\n")

cat("=== Importances relatives sans smolt_k ===\n")
print(sw(stj_1sw_nok_avg_model))

cat("\n=== Importances relatives avec smolt_k ===\n")
print(sw(stj_1sw_avg_model))

# Attendu sans smolt_k : ~15 modèles retenus, poids max ≈ 0.11 -> sélection
# chaotique, pas de signal clair au-delà de cohort


# =============================================================================
#
#   VISUALISATION — EFFETS MARGINAUX AVEC VS SANS smolt_k
#
#   Comparaison visuelle des effets marginaux des variables communes
#   (BelleIsle_TH, cohort, plankton_index) dans les deux modèles.
#   Si les courbes sont similaires, smolt_k ne masque pas les autres effets.
#
# =============================================================================

set.seed(42)

vars_communes <- c("BelleIsle_TH", "cohort", "plankton_index")

for (v in vars_communes) {
  
  # Model averaging avec smolt_k
  pd_avg_avec <- build_marginal_plot_data_avg(
    stj_1sw_top_models, v,
    scaled_stj_1sw, stj_1sw,
    stj_1sw_center, stj_1sw_scale,
    n_points = 10, n_boot = 20
  )
  
  # Model averaging sans smolt_k
  pd_avg_sans <- build_marginal_plot_data_avg(
    stj_1sw_nok_top_models, v,
    scaled_stj_1sw, stj_1sw,
    stj_1sw_center, stj_1sw_scale,
    n_points = 10, n_boot = 20
  )
  
  cohort_mid <- median(stj_1sw$cohort, na.rm = TRUE)
  
  p_avec <- build_marginal_ggplot_avg(pd_avg_avec) +
    labs(subtitle = "Avec smolt_k") +
    theme(plot.subtitle = element_text(hjust = 0.5, colour = "steelblue4",
                                       face = "bold"))
  
  # Reconstruire manuellement en orange pour le modèle sans smolt_k
  p_sans <- ggplot(pd_avg_sans$pred_avg, aes(x, fit)) +
    {if (length(pd_avg_sans$pred_models) > 0)
      geom_line(data = dplyr::bind_rows(pd_avg_sans$pred_models),
                aes(x, fit, group = model, linetype = model),
                colour = "grey50", linewidth = 0.8)} +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.20, fill = "darkorange") +
    geom_line(linewidth = 1.4, colour = "darkorange3") +
    geom_point(data = pd_avg_sans$obs_df, aes(x, y, fill = cohort),
               inherit.aes = FALSE,
               size = 3.5, colour = "black", shape = 21, stroke = 0.5) +
    scale_fill_gradient2(low = "#00897B", mid = "#BDBDBD", high = "#8E24AA",
                         midpoint = cohort_mid, name = "Cohort") +
    geom_text_repel(data = pd_avg_sans$obs_df, aes(x, y, label = cohort),
                    inherit.aes = FALSE, size = 3, max.overlaps = 20) +
    labs(x = v, y = "1SW return rate",
         subtitle = "Sans smolt_k") +
    theme_classic(base_size = 13) +
    theme(plot.subtitle = element_text(hjust = 0.5, colour = "darkorange3",
                                       face = "bold"))
  
  print(
    (p_avec + p_sans + plot_layout(guides = "collect")) +
      plot_annotation(title = paste("Variable :", v),
                      subtitle = "Bleu = model averaging avec smolt_k  |  Orange = model averaging sans smolt_k")
  )
}


# =============================================================================
#
#   RÉSUMÉ DE L'ANALYSE
#
#   Test 1 — Tendance temporelle :
#     smolt_k est fortement corrélé avec cohort (r ≈ -0.58, p < 0.001).
#     Les deux variables déclinent ensemble dans le temps.
#
#   Test 2 — Corrélation partielle :
#     Une fois cohort contrôlé, la corrélation smolt_k ~ return_1sw
#     disparaît (r partiel ≈ -0.17, p ≈ 0.38). La relation brute est
#     donc spurieuse — elle provient de la covariation temporelle commune.
#
#   Test 3A — Comparaison des coefficients :
#     BelleIsle_TH et plankton_index changent modérément sans smolt_k
#     (delta ~23-45%) ; cohort absorbe davantage de variance temporelle.
#     Les directions sont conservées — pas d'inversion de signe.
#     BelleIsle_TH perd sa significativité sans smolt_k (p = 0.126 vs
#     p = 0.015), ce qui indique une amplification partielle de son effet
#     apparent par smolt_k.
#
#   Test 3B — Sélection de variables :
#     Sans smolt_k, ~15 modèles compétiteurs avec des poids faibles
#     (max ≈ 0.11) — la sélection devient chaotique. Avec smolt_k,
#     3 modèles retenus avec une structure claire (M1 poids = 0.34).
#
#   QUESTION :
#     Faut-il conserver smolt_k malgré sa relation spurieuse, en
#     l'interprétant avec prudence dans la discussion, ou le retirer
#     en assumant l'instabilité de la sélection qui en découle ?
#
# =============================================================================


