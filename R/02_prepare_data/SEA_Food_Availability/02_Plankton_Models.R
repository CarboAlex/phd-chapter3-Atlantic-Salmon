


# ══════════════════════════════════════════════════════════════════════════════
#
# Fits a GAM on log-transformed monthly CPR zooplankton biomass to model
# temporal, seasonal and spatial variation, then generates an annual median
# biomass index for the Labrador Sea from 100 fixed prediction points, with
# LOO-CV validation. Exports annual plankton indices for STJ and TRI, for both
# 1SW and 2SW return types.
#
# INPUTS:
#   - R/01_derived_data/Food_Availability/cpr_data_raw.fst
#       CPR zooplankton biomass data (Sample_Id, Latitude, Longitude,
#       Midpoint_Date_Local, biomass_g)
#   - R/00_raw_data/ReferenceMigrationLayers/LabradorSea.shp
#       Polygon defining the Labrador Sea spatial domain
#
# OUTPUTS:
#   - R/03_clean_data/Sea/Plankton_LabradorSea_stj_1SW.fst
#   - R/03_clean_data/Sea/Plankton_LabradorSea_tri_1SW.fst
#   - R/03_clean_data/Sea/Plankton_LabradorSea_stj_2SW.fst
#   - R/03_clean_data/Sea/Plankton_LabradorSea_tri_2SW.fst
#       Annual plankton index by river and return type
# ══════════════════════════════════════════════════════════════════════════════


library(fst)
library(dplyr)
library(tidyr)
library(ggplot2)
library(mgcv)
library(DHARMa)
library(sf)
library(rnaturalearth)
library(lubridate)


# ── Working directory ──────────────────────────────────────────────────────
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

# SOURCE: CPR (Continuous Plankton Recorder) zooplankton biomass data
#         provided by the Atlantic Zone Monitoring Program (AZMP / DFO).
cpr <- read.fst(
  file.path(base_path, "R/01_derived_data/Food_Availability/cpr_data_raw.fst")
)

# Aggregate biomass by sample (sum across taxa per sample/date/location)
cpr <- cpr %>%
  dplyr::mutate(Midpoint_Date_Local = as.Date(Midpoint_Date_Local)) %>%
  dplyr::group_by(Sample_Id, Latitude, Longitude, Midpoint_Date_Local) %>%
  dplyr::summarise(sum_biomass_g = sum(biomass_g, na.rm = TRUE),
                   .groups = "drop_last") %>%
  dplyr::mutate(year = year(Midpoint_Date_Local)) %>%
  dplyr::select(-year)


# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

n_points   <- 100        # Number of fixed prediction points in Labrador Sea
output_dir <- file.path(base_path, "R/03_clean_data/Sea")


# ══════════════════════════════════════════════════════════════════════════════
# RESTRICT TO LABRADOR SEA POLYGON
# ══════════════════════════════════════════════════════════════════════════════

# Keep only CPR samples within the Labrador Sea polygon. 

lab_sea <- sf::st_transform(
  sf::st_read(file.path(base_path,
                        "R/00_raw_data/ReferenceMigrationLayers/LabradorSea.shp"),
              quiet = TRUE),
  crs = 4326) |> sf::st_make_valid()

cpr_sf <- sf::st_as_sf(cpr, coords = c("Longitude", "Latitude"),
                       crs = 4326, remove = FALSE) |> sf::st_make_valid()
inside <- sf::st_intersects(cpr_sf, lab_sea, sparse = FALSE)
cpr    <- cpr[which(inside), ]


# ══════════════════════════════════════════════════════════════════════════════
# MONTHLY AGGREGATION
# ══════════════════════════════════════════════════════════════════════════════

# Mean monthly biomass across ALL samples in the region (one value per
# year-month), on log(mean+1) scale. 

cpr_monthly <- cpr %>%
  dplyr::mutate(
    YEAR  = lubridate::year(Midpoint_Date_Local),
    Month = lubridate::month(Midpoint_Date_Local)
  ) %>%
  dplyr::group_by(YEAR, Month) %>%
  dplyr::summarise(
    biomass  = log(mean(sum_biomass_g, na.rm = TRUE) + 1),   # log(mean+1)
    mean_raw = mean(sum_biomass_g, na.rm = TRUE),
    mLat     = mean(Latitude,  na.rm = TRUE),
    mLon     = mean(Longitude, na.rm = TRUE),
    n_samples = dplyr::n(),
    .groups  = "drop"
  )


# ══════════════════════════════════════════════════════════════════════════════
# RESPONSE VARIABLE DISTRIBUTION
# ══════════════════════════════════════════════════════════════════════════════

par(mfrow = c(1, 2))
hist(cpr_monthly$mean_raw, breaks = 40,
     main = "Monthly-mean biomass", xlab = "Biomass (g/sample)")
hist(cpr_monthly$biomass, breaks = 40,
     main = "log(mean biomass + 1)", xlab = "log(mean + 1)")
par(mfrow = c(1, 1))

# Two families are compared below (Tweedie vs log-Gaussian). The response here
# is already log(mean+1); for the Tweedie comparison we refit on the raw
# monthly mean (mean_raw) to keep the family comparison meaningful.


# ══════════════════════════════════════════════════════════════════════════════
# FAMILY SELECTION — CANDIDATE MODELS
# ══════════════════════════════════════════════════════════════════════════════

# [A] Tweedie on raw monthly-mean biomass
m_tw <- gam(mean_raw ~
              s(YEAR,  bs = "tp", k = 10) +
              s(Month, bs = "cc", k = 12) +
              te(mLat, mLon),
            data    = cpr_monthly,
            family  = tw(),
            control = gam.control(nthreads = 4))
summary(m_tw)

# [B] Gaussian on log(mean+1) (final candidate)
m_gauss <- gam(biomass ~
                 s(YEAR,  bs = "tp", k = 10) +
                 s(Month, bs = "cc", k = 12) +
                 te(mLat, mLon),
               data    = cpr_monthly,
               family  = gaussian(),
               control = gam.control(nthreads = 4))
summary(m_gauss)


# ══════════════════════════════════════════════════════════════════════════════
# DHARMA DIAGNOSTICS — FAMILY SELECTION
# ══════════════════════════════════════════════════════════════════════════════

set.seed(123)
sim_tw    <- simulateResiduals(m_tw,    n = 500)
sim_gauss <- simulateResiduals(m_gauss, n = 500)

plot(sim_tw,    title = "Tweedie")
plot(sim_gauss, title = "Gaussian (log)")

resultats_dharma <- data.frame(
  Modele = c("Tweedie", "Gaussian_log"),
  Dispersion_ratio = c(
    testDispersion(sim_tw,    plot = FALSE)$statistic,
    testDispersion(sim_gauss, plot = FALSE)$statistic),
  Dispersion_pvalue = c(
    testDispersion(sim_tw,    plot = FALSE)$p.value,
    testDispersion(sim_gauss, plot = FALSE)$p.value),
  ZeroInflation_ratio = c(
    testZeroInflation(sim_tw,    plot = FALSE)$statistic,
    testZeroInflation(sim_gauss, plot = FALSE)$statistic),
  ZeroInflation_pvalue = c(
    testZeroInflation(sim_tw,    plot = FALSE)$p.value,
    testZeroInflation(sim_gauss, plot = FALSE)$p.value)
)
resultats_dharma[, -1] <- round(resultats_dharma[, -1], 4)
resultats_dharma

# The log-Gaussian model is retained (Imlay/Tara approach); the ZI=Inf for the
# Gaussian is a transformation artefact, not a deficiency.


# ══════════════════════════════════════════════════════════════════════════════
# TEMPORAL RESIDUAL CHECK — FINAL MODEL
# ══════════════════════════════════════════════════════════════════════════════

# A valid annual index requires mean residual bias to be constant over time.

cpr_monthly$resid_tw    <- residuals(m_tw,    type = "deviance")
cpr_monthly$resid_gauss <- residuals(m_gauss, type = "deviance")

resid_temp <- cpr_monthly %>%
  dplyr::group_by(YEAR) %>%
  dplyr::summarise(
    resid_tw_moy    = mean(resid_tw,    na.rm = TRUE),
    resid_gauss_moy = mean(resid_gauss, na.rm = TRUE)
  )

resid_temp %>%
  tidyr::pivot_longer(cols      = c(resid_tw_moy, resid_gauss_moy),
                      names_to  = "modele",
                      values_to = "resid_moy") %>%
  dplyr::mutate(modele = dplyr::recode(modele,
                                       resid_tw_moy    = "Tweedie",
                                       resid_gauss_moy = "Gaussian (log)")) %>%
  ggplot(aes(x = YEAR, y = resid_moy)) +
  geom_point(size = 1.5) + geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  theme_minimal(base_size = 14) +
  facet_wrap(~modele, ncol = 1, scales = "free_y") +
  labs(title    = "Mean annual residuals — residual temporal pattern",
       subtitle = "Loess near y=0 without trend = valid annual index",
       y        = "Mean deviance residual", x = "Year")

lm(resid_tw_moy    ~ YEAR, data = resid_temp) |> summary()
lm(resid_gauss_moy ~ YEAR, data = resid_temp) |> summary()
t.test(resid_temp$resid_tw_moy,    mu = 0)
t.test(resid_temp$resid_gauss_moy, mu = 0)


# ══════════════════════════════════════════════════════════════════════════════
# GAM DIAGNOSTICS — FINAL MODEL
# ══════════════════════════════════════════════════════════════════════════════

par(mfrow = c(2, 2))
gam.check(m_gauss)
par(mfrow = c(1, 1))


# ══════════════════════════════════════════════════════════════════════════════
# FINAL MODEL SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

summary(m_gauss)

par(mfrow = c(1, 3))
plot(m_gauss, select = 1, shade = TRUE,
     main = "Temporal trend — s(YEAR)",
     ylab = "Partial effect on log(mean+1)")
plot(m_gauss, select = 2, shade = TRUE,
     main = "Seasonality — s(Month)", ylab = "Partial effect")
plot(m_gauss, select = 3, scheme = 2,
     main = "Spatial surface — te(mLat, mLon)")
par(mfrow = c(1, 1))


# ══════════════════════════════════════════════════════════════════════════════
# PREDICTION POINT SELECTION
# ══════════════════════════════════════════════════════════════════════════════

# 100 fixed prediction points sampled from CPR locations within the Labrador Sea
# polygon (refinement over Imlay's single most-sampled station).

world <- ne_countries(scale = "medium", returnclass = "sf")

cpr_sf_in <- sf::st_as_sf(cpr, coords = c("Longitude", "Latitude"),
                          crs = 4326, remove = FALSE) |> sf::st_make_valid()

set.seed(10)
ss    <- sample(seq_len(nrow(cpr_sf_in)), n_points)
ss_df <- as.data.frame(cpr_sf_in[ss, ])[, c("Longitude", "Latitude")]

# ── Map of prediction points ──────────────────────────────────────────────────
lab_sea_water <- sf::st_difference(lab_sea, sf::st_union(world) |>
                                     sf::st_make_valid()) |> sf::st_make_valid()

ggplot() +
  geom_sf(data = world, fill = "gray89", color = "#AAAAAA", linewidth = 0.25) +
  geom_sf(data = lab_sea_water, fill = "#4B9FD4", alpha = 0.20,
          color = "#2C6E9E", linewidth = 0.3, linetype = "dashed") +
  geom_point(data = cpr, aes(Longitude, Latitude),
             color = "red", size = 0.8, alpha = 0.4) +
  geom_point(data = ss_df, aes(Longitude, Latitude),
             shape = 8, color = "black", size = 2.5, stroke = 1) +
  coord_sf(xlim = c(-71, -38), ylim = c(44, 62), expand = FALSE) +
  theme_bw(base_size = 13) +
  labs(x = "Longitude", y = "Latitude",
       title = "Labrador Sea — CPR samples (red) and 100 prediction points")


# ══════════════════════════════════════════════════════════════════════════════
# ANNUAL BIOMASS PREDICTION — 100 FIXED POINTS
# ══════════════════════════════════════════════════════════════════════════════

# For each fixed point, predict monthly biomass for every year x month, sum the
# 12 monthly values within each year (Imlay: "summed the monthly values across
# each year"), then take the median across the 100 points.
#
# NOTE: the GAM was fit on monthly-mean positions, so te(mLat, mLon) is only
# well-constrained near the mean sampling locations. Predicting at 100 dispersed
# points partly extrapolates the spatial term; expect the 100 predictions to be
# similar to one another (the index is driven by s(YEAR)/s(Month), not space).
# The spread across points is therefore best read as spatial uncertainty.
#
# Model predicts log(mean+1); back-transform via exp(pred) - 1.

# Continuous year grid (min..max) so that unsampled years (e.g. 1987-1990)
# are predicted by the GAM rather than skipped. The s(YEAR) smooth interpolates
# these gaps, as in Imlay et al. (2024) who predict over a fixed 1970-2021 grid.
years_all <- seq(min(cpr_monthly$YEAR), max(cpr_monthly$YEAR), by = 1L)


newdata <- expand.grid(
  YEAR  = years_all,
  Month = 1:12,
  it    = seq_len(nrow(ss_df))
) %>%
  dplyr::mutate(
    mLat = ss_df$Latitude[it],
    mLon = ss_df$Longitude[it]
  )

pred_all <- predict(m_gauss, newdata = newdata, type = "link", se.fit = TRUE)

plankton <- newdata %>%
  dplyr::mutate(
    fit    = exp(pred_all$fit) - 1,
    fit_lo = exp(pred_all$fit - 1.96 * pred_all$se.fit) - 1,
    fit_hi = exp(pred_all$fit + 1.96 * pred_all$se.fit) - 1
  ) %>%
  dplyr::group_by(it, YEAR) %>%
  dplyr::summarise(
    mLat        = first(mLat),
    mLon        = first(mLon),
    Plankton    = sum(fit,    na.rm = TRUE),   # sum of 12 monthly values
    Plankton_lo = sum(fit_lo, na.rm = TRUE),
    Plankton_hi = sum(fit_hi, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  dplyr::rename(year = YEAR)

summ_plankton <- plankton %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(
    zoo_median = median(Plankton),
    zoo_mean   = mean(Plankton),
    zoo_sd     = sd(Plankton),
    zoo_q25    = quantile(Plankton, 0.25),
    zoo_q75    = quantile(Plankton, 0.75),
    .groups    = "drop"
  )


# ── Visualisation ─────────────────────────────────────────────────────────────
gap_years <- setdiff(min(plankton$year):max(plankton$year),
                     unique(lubridate::year(cpr$Midpoint_Date_Local)))
if (length(gap_years) > 0) {
  gap_xmin <- min(gap_years) - 0.5
  gap_xmax <- max(gap_years) + 0.5
}

ggplot() +
  geom_line(data = plankton,
            aes(x = year, y = Plankton, group = it),
            linewidth = 0.5, alpha = 0.3) +
  geom_ribbon(data = summ_plankton,
              aes(x = year, ymin = zoo_q25, ymax = zoo_q75),
              fill = "red", alpha = 0.2) +
  geom_line(data = summ_plankton,
            aes(x = year, y = zoo_median),
            linewidth = 1.2, color = "red") +
  geom_point(data = summ_plankton,
             aes(x = year, y = zoo_median),
             size = 2, color = "red") +
  {if (length(gap_years) > 0)
    annotate("rect", xmin = gap_xmin, xmax = gap_xmax,
             ymin = -Inf, ymax = Inf, fill = "grey40", alpha = 0.3)} +
  theme_bw(base_size = 14) +
  labs(x = "Year", y = "Plankton biomass (summed monthly means)",
       title = "Annual biomass index — Labrador Sea (Imlay/Tara-style)",
       subtitle = "Red: median of 100 points | Ribbon: IQR")


# ══════════════════════════════════════════════════════════════════════════════
# TEMPORAL LOO CROSS-VALIDATION — ANNUAL MEDIAN INDEX
# ══════════════════════════════════════════════════════════════════════════════

ref_index <- summ_plankton %>%
  dplyr::select(year, zoo_median_ref = zoo_median)

set.seed(42)
loo_index <- lapply(years_all, function(yr) {
  cat("  Annual index for", yr, "...\n")
  
  train <- cpr_monthly %>% dplyr::filter(YEAR != yr)
  
  m_cv <- tryCatch(
    gam(biomass ~
          s(YEAR,  bs = "tp", k = 10) +
          s(Month, bs = "cc", k = 12) +
          te(mLat, mLon),
        data    = train,
        family  = gaussian(),
        control = gam.control(nthreads = 4)),
    error = function(e) NULL
  )
  if (is.null(m_cv)) return(NULL)
  
  nd <- expand.grid(Month = 1:12, it = seq_len(nrow(ss_df))) %>%
    dplyr::mutate(
      YEAR = yr,
      mLat = ss_df$Latitude[it],
      mLon = ss_df$Longitude[it]
    )
  
  pred <- tryCatch(
    predict(m_cv, newdata = nd, type = "link", se.fit = FALSE),
    error = function(e) NULL
  )
  if (is.null(pred)) return(NULL)
  
  plankton_cv <- nd %>%
    dplyr::mutate(fit = pmax(exp(pred) - 1, 0)) %>%
    dplyr::group_by(it) %>%
    dplyr::summarise(Plankton = sum(fit, na.rm = TRUE), .groups = "drop")
  
  data.frame(
    year_removed = yr,
    index_cv     = median(plankton_cv$Plankton),
    index_ref    = ref_index$zoo_median_ref[ref_index$year == yr]
  )
})

loo_index <- dplyr::bind_rows(loo_index) %>%
  dplyr::mutate(diff_pct = (index_cv - index_ref) / index_ref * 100)

metrics_index <- loo_index %>%
  dplyr::summarise(
    R2   = cor(index_cv, index_ref)^2,
    RMSE = sqrt(mean((index_cv - index_ref)^2)),
    Bias = mean(index_cv - index_ref),
    MAPE = mean(abs(diff_pct))
  )

cat("\n=== TEMPORAL LOO — ANNUAL MEDIAN INDEX ===\n")
cat(sprintf("  R2   : %.3f\n", metrics_index$R2))
cat(sprintf("  RMSE : %.1f\n", metrics_index$RMSE))
cat(sprintf("  Bias : %.1f\n", metrics_index$Bias))
cat(sprintf("  MAPE : %.1f%%\n", metrics_index$MAPE))

ggplot(loo_index, aes(x = index_ref, y = index_cv)) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_point(size = 3, color = "steelblue") +
  geom_text(aes(label = year_removed),
            hjust = -0.2, vjust = 0.5, size = 3.5, color = "grey40") +
  theme_bw(base_size = 13) +
  labs(title    = "Temporal LOO — Annual median index",
       subtitle = "Predicted index (year removed) vs full-model reference index",
       x        = "Reference index", y = "LOO predicted index")



# ══════════════════════════════════════════════════════════════════════════════
# DECISION TEST: PCA multi-year index  vs  year-t only
#
# If consecutive sea-years (t, t+1, t+2) are too strongly correlated,
# a PCA on them is meaningless (PC1 would just replicate a single variable).
# ══════════════════════════════════════════════════════════════════════════════

# ── Build the annual series with lagged years ─────────────────────────────────
lag_tbl <- summ_plankton %>%
  dplyr::arrange(year) %>%
  dplyr::transmute(
    year,
    t   = zoo_median,
    tp1 = dplyr::lag(zoo_median, 1),   # year t+1 relative to entry year
    tp2 = dplyr::lag(zoo_median, 2)    # year t+2
  )

# ── 1SW relevant pair: t, t+1 ─────────────────────────────────────────────────
pair_1sw <- lag_tbl %>% dplyr::select(t, tp1) %>% na.omit()
r_1sw    <- cor(pair_1sw$t, pair_1sw$tp1)

# 1SW: consecutive-year correlation (t, t+1)
r_1sw


# ── 2SW relevant triplet: t, t+1, t+2 ─────────────────────────────────────────
trip_2sw <- lag_tbl %>% dplyr::select(t, tp1, tp2) %>% na.omit()
cor_2sw  <- cor(trip_2sw)

# 2SW: consecutive-year correlations (t, t+1, t+2)
round(cor_2sw, 3)


# Only year t will be used, since years t, t+1 and t+2 are too strongly correlated. 


# ══════════════════════════════════════════════════════════════════════════════
# EXPORT
# ══════════════════════════════════════════════════════════════════════════════

plankton_base <- summ_plankton %>%
  dplyr::select(year, zoo_median)

# 1SW: plankton in year t
plankton_1sw <- plankton_base %>%
  dplyr::rename(Plankton_t = zoo_median) %>%
  dplyr::mutate(cohort = year) %>%
  dplyr::select(cohort, Plankton_t)

for (rv in c("stj", "tri")) {
  out <- plankton_1sw %>%
    dplyr::mutate(river = rv, sw_type = "1SW") %>%
    dplyr::select(river, sw_type, cohort, Plankton_t)
  fname <- paste0("Plankton_LabradorSea_", rv, "_1SW.fst")
  write_fst(out, path = file.path(output_dir, fname))
  message("Saved: ", fname)
}

# 2SW: plankton in year t and year t+1
plankton_t <- plankton_base %>%
  dplyr::rename(Plankton_t = zoo_median) %>%
  dplyr::mutate(cohort = year) %>%
  dplyr::select(cohort, Plankton_t)

plankton_tp1 <- plankton_base %>%
  dplyr::rename(Plankton_tp1 = zoo_median) %>%
  dplyr::mutate(cohort = year - 1L) %>%
  dplyr::select(cohort, Plankton_tp1)

plankton_2sw <- plankton_t %>%
  dplyr::left_join(plankton_tp1, by = "cohort")

for (rv in c("stj", "tri")) {
  out <- plankton_2sw %>%
    dplyr::mutate(river = rv, sw_type = "2SW") %>%
    dplyr::select(river, sw_type, cohort, Plankton_t, Plankton_tp1)
  fname <- paste0("Plankton_LabradorSea_", rv, "_2SW.fst")
  write_fst(out, path = file.path(output_dir, fname))
  message("Saved: ", fname)
}




