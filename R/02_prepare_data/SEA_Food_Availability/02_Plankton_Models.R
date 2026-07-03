


# ══════════════════════════════════════════════════════════════════════════════
# SCRIPT: Zooplankton Biomass Index — Labrador Sea
#
# Fits a GAM on log-transformed CPR zooplankton biomass to model temporal,
# seasonal, spatial, and SST effects. Generates an annual median biomass
# index for the Labrador Sea from 100 fixed prediction points, with LOO-CV
# validation. Exports annual plankton indices for STJ and TRI, for both
# 1SW and 2SW return types.
#
# INPUTS:
#   - R/01_derived_data/Food_Availability/cpr_data_raw.fst
#       CPR zooplankton biomass data (Sample_Id, Latitude, Longitude,
#       Midpoint_Date_Local, biomass_g)
#   - R/00_raw_data/SST_ORAS5_Labrador/ESACCI_SST_L4_YYYY_MM.nc
#       Daily SST rasters (ESA CCI L4, zipped .nc files)
#   - R/00_raw_data/ReferenceMigrationLayers/LabradorSea.shp
#       Polygon defining the Labrador Sea spatial domain
#
# OUTPUTS:
#   - R/01_derived_data/Food_Availability/CPR_SST_extracted.fst
#       CPR dataset with SST_14d, SST_30d, SST_60d columns appended
#   - R/01_derived_data/Food_Availability/SST_prediction_100pts.rds
#       Daily SST_30d matrix for 100 fixed prediction points
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
library(terra)
library(sf)
library(rnaturalearth)
library(pbapply)
library(lubridate)
library(parallel)


# ── Working directory ──────────────────────────────────────────────────────
# Automatically detects the repository root based on the script location.
# No hardcoded path — works on any machine.
current_file <- rstudioapi::getActiveDocumentContext()$path
path_parts   <- unlist(strsplit(normalizePath(current_file, winslash = "/"), "/"))
target_index <- which(path_parts == "phd-chapter3-Atlantic-Salmon")
base_path    <- paste(path_parts[1:target_index], collapse = "/")
setwd(base_path)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: DATA IMPORT
# ══════════════════════════════════════════════════════════════════════════════

# SOURCE: CPR (Continuous Plankton Recorder) zooplankton biomass data
#         provided by the Atlantic Zone Monitoring Program (AZMP / DFO).
# STATUS: Confidential — not included in this repository.
#         An example file with identical structure and simulated data is
#         available at the same path with the suffix _EXAMPLE.fst.
#         To request access, contact the corresponding author once the
#         associated manuscript is published.
cpr <- read.fst(
  file.path(base_path, "R/01_derived_data/Food_Availability/cpr_data_raw.fst")
)

# Aggregate biomass by sample (sum across taxa per sample/date/location)
# and restrict to years >= 1982 to match the SST data availability.
cpr <- cpr %>%
  dplyr::mutate(Midpoint_Date_Local = as.Date(Midpoint_Date_Local)) %>%
  dplyr::group_by(Sample_Id, Latitude, Longitude, Midpoint_Date_Local) %>%
  dplyr::summarise(sum_biomass_g = sum(biomass_g, na.rm = TRUE),
                   .groups = "drop_last") %>%
  dplyr::mutate(year = year(Midpoint_Date_Local)) %>%
  dplyr::filter(year >= 1982) %>%
  dplyr::select(-year)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

sst_dir        <- file.path(base_path, "R/00_raw_data/SST_ORAS5_Labrador")
radius_m       <- 30000     # Spatial buffer radius for SST extraction (30 km)
windows        <- c(SST_14d = 14L, SST_30d = 30L, SST_60d = 60L)
n_cores        <- 4L        # Parallel workers for SST extraction
window_days    <- 30L       # Rolling window for prediction SST (days)
n_points       <- 100       # Number of fixed prediction points in Labrador Sea

sst_output_path  <- file.path(base_path,
                              "R/01_derived_data/Food_Availability/CPR_SST_extracted.fst")
checkpoint_dir   <- file.path(base_path,
                              "R/01_derived_data/Food_Availability/SST_prediction_checkpoints")
final_sst_path   <- file.path(base_path,
                              "R/01_derived_data/Food_Availability/SST_prediction_100pts.rds")
output_dir       <- file.path(base_path, "R/01_derived_data/Food_Availability")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3: RESPONSE VARIABLE DISTRIBUTION
# ══════════════════════════════════════════════════════════════════════════════

# ── Zero frequency assessment ────────────────────────────────────────────────
# Evaluate whether zeros are frequent enough and structurally distinct from
# positive values to justify a zero-inflated model (e.g., ZI-Tweedie, hurdle).

prop_zeros <- mean(cpr$sum_biomass_g == 0, na.rm = TRUE)
round(prop_zeros * 100, 2)

# ── Distribution plots ───────────────────────────────────────────────────────
par(mfrow = c(1, 3))

hist(cpr$sum_biomass_g, breaks = 100,
     main = "Raw distribution",
     xlab = "Biomass (g/sample)")

hist(log(cpr$sum_biomass_g + 1), breaks = 100,
     main = "Distribution log(biomass + 1)",
     xlab = "log(Biomass + 1)")

hist(log(cpr$sum_biomass_g[cpr$sum_biomass_g > 0]), breaks = 100,
     main = "Distribution log(biomass) | biomass > 0",
     xlab = "log(Biomass)")

par(mfrow = c(1, 1))

# Two families are retained for comparison:
#   [A] Tweedie           — natively handles zeros + continuous positive values
#   [B] Gaussian log(y+1) — justified if log(y|y>0) is approximately normal
#                           and zeros likely represent very low biomass rather
#                           than true absences (CPR tows span 18 km; complete
#                           plankton absence over 18 km is biologically
#                           implausible). The objective is a relative annual
#                           index at Labrador Sea scale, so the distinction
#                           between near-zero and zero biomass is negligible.


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4: SST EXTRACTION FOR CPR POINTS
# ══════════════════════════════════════════════════════════════════════════════

# ── Overview ─────────────────────────────────────────────────────────────────
# For each CPR sample, extract mean SST within a 30 km buffer over three
# rolling windows ending on the sampling date (14, 30, and 60 days prior).
#
# Justification for a 30 km radius: pixel resolution is ~43 x 46 km at ~50 N.
# A radius of 18 km (CPR tow length) risks capturing only a partial pixel.
# A 30 km radius guarantees at least one complete pixel and typically 2-3,
# consistent with the data resolution.

if (file.exists(sst_output_path)) {
  
  message("SST file already available — loading without recalculation.")
  cpr <- read_fst(sst_output_path)
  message(sprintf("  %d rows loaded.", nrow(cpr)))
  
} else {
  
  # ── Build daily NC file index ──────────────────────────────────────────────
  # Each .nc is a ZIP containing one NC per day of the month.
  # Build a data.frame (date | zip_path | nc_name) covering all dates
  # needed for the three windows across all CPR points.
  
  cpr_dates     <- unique(cpr$Midpoint_Date_Local)
  earliest_need <- min(cpr_dates) - max(windows)
  latest_need   <- max(cpr_dates)
  
  zip_files <- list.files(sst_dir,
                          pattern    = "ESACCI_SST_L4_\\d{4}_\\d{2}\\.nc$",
                          full.names = TRUE)
  
  file_index <- lapply(zip_files, function(f) {
    bname <- basename(f)
    yr_k  <- as.integer(sub("ESACCI_SST_L4_(\\d{4})_(\\d{2})\\.nc", "\\1", bname))
    mo_k  <- as.integer(sub("ESACCI_SST_L4_(\\d{4})_(\\d{2})\\.nc", "\\2", bname))
    if (is.na(yr_k) || is.na(mo_k)) return(NULL)
    
    month_start <- as.Date(sprintf("%04d-%02d-01", yr_k, mo_k))
    month_end   <- as.Date(sprintf("%04d-%02d-%02d", yr_k, mo_k,
                                   lubridate::days_in_month(month_start)))
    if (month_end < earliest_need || month_start > latest_need) return(NULL)
    
    contents <- tryCatch(unzip(f, list = TRUE), error = function(e) NULL)
    if (is.null(contents) || nrow(contents) == 0L) return(NULL)
    
    nc_names  <- contents$Name
    date_strs <- substr(basename(nc_names), 1, 8)
    dates_i   <- as.Date(date_strs, format = "%Y%m%d")
    valid     <- !is.na(dates_i) & dates_i >= earliest_need & dates_i <= latest_need
    
    if (!any(valid)) return(NULL)
    
    data.frame(date     = dates_i[valid],
               zip_path = f,
               nc_name  = nc_names[valid],
               stringsAsFactors = FALSE)
  })
  
  file_index <- do.call(rbind, file_index[!sapply(file_index, is.null)])
  file_index <- file_index[order(file_index$date), ]
  rownames(file_index) <- NULL
  
  # ── Prepare CPR points and spatial buffers ─────────────────────────────────
  # Project to metric CRS, buffer at 30 km, reproject to WGS84.
  # Buffers are computed once to reduce computation time.
  
  cpr_vect <- terra::vect(
    as.data.frame(cpr),
    geom = c("Longitude", "Latitude"),
    crs  = "EPSG:4326"
  )
  
  cpr_m    <- terra::project(cpr_vect, "EPSG:32621")
  bufs_m   <- terra::buffer(cpr_m, width = radius_m)
  bufs_wgs <- terra::project(bufs_m, "EPSG:4326")
  
  ext_global <- terra::ext(bufs_wgs)
  
  # ── Accumulation matrices ──────────────────────────────────────────────────
  # Accumulate SST sum and valid-day count per point per window.
  # Final mean = sum / n_valid_days.
  
  n_pts       <- nrow(cpr)
  sst_sum     <- matrix(0,  nrow = n_pts, ncol = length(windows))
  sst_ndays   <- matrix(0L, nrow = n_pts, ncol = length(windows))
  colnames(sst_sum)   <- names(windows)
  colnames(sst_ndays) <- names(windows)
  
  window_starts <- lapply(windows, function(w) cpr$Midpoint_Date_Local - (w - 1L))
  window_ends   <- cpr$Midpoint_Date_Local
  
  # ── Main extraction loop: one pass per daily NC file ──────────────────────
  # For each day, identify which points require that day in at least one
  # window, load the SST raster once, and extract the spatial mean per buffer.
  
  n_files <- nrow(file_index)
  t_start <- proc.time()
  
  for (fi in seq_len(n_files)) {
    
    day_date <- file_index$date[fi]
    
    # A point needs day d if: window_start[w][i] <= d <= window_end[i]
    # Use the widest window (60d) as the outer mask.
    in_max_window <- (day_date >= window_starts[[length(windows)]]) &
      (day_date <= window_ends)
    pts_needed <- which(in_max_window)
    
    if (length(pts_needed) == 0L) next
    
    tmp_dir <- tempfile()
    dir.create(tmp_dir)
    
    r_sst <- tryCatch({
      extracted <- unzip(file_index$zip_path[fi],
                         files = file_index$nc_name[fi],
                         exdir = tmp_dir)
      if (length(extracted) == 0L) stop("unzip empty")
      
      r_all   <- terra::rast(extracted[1])
      sst_idx <- grep("analysed_sst$", names(r_all), value = FALSE)
      if (length(sst_idx) == 0L) sst_idx <- 1L
      r_sst_k <- r_all[[sst_idx[1]]] - 273.15   # Kelvin to Celsius
      
      terra::crs(r_sst_k) <- "EPSG:4326"
      terra::crop(r_sst_k, ext_global)
      
    }, error = function(e) {
      message(sprintf("  [SKIP] %s : %s",
                      basename(file_index$zip_path[fi]),
                      conditionMessage(e)))
      NULL
    }, finally = {
      unlink(tmp_dir, recursive = TRUE)
    })
    
    if (is.null(r_sst)) next
    
    extracted_vals <- tryCatch(
      terra::extract(r_sst, bufs_wgs[pts_needed, ],
                     fun = mean, na.rm = TRUE)[, 2],
      error = function(e) rep(NA_real_, length(pts_needed))
    )
    
    for (w in seq_along(windows)) {
      in_w   <- day_date >= window_starts[[w]][pts_needed]
      pts_w  <- pts_needed[in_w]
      vals_w <- extracted_vals[in_w]
      
      valid_v <- !is.na(vals_w)
      if (any(valid_v)) {
        sst_sum  [pts_w[valid_v], w] <- sst_sum  [pts_w[valid_v], w] + vals_w[valid_v]
        sst_ndays[pts_w[valid_v], w] <- sst_ndays[pts_w[valid_v], w] + 1L
      }
    }
    
    if (fi %% 100 == 0 || fi == n_files) {
      elapsed <- (proc.time() - t_start)["elapsed"]
      message(sprintf("  [%d/%d] %s — %.0fs elapsed", fi, n_files, day_date, elapsed))
    }
  }
  
  # ── Compute means and append to CPR dataset ────────────────────────────────
  
  sst_means <- sst_sum / ifelse(sst_ndays > 0L, sst_ndays, NA_real_)
  
  cpr$SST_14d <- sst_means[, "SST_14d"]
  cpr$SST_30d <- sst_means[, "SST_30d"]
  cpr$SST_60d <- sst_means[, "SST_60d"]
  
  message("\n=== MISSING VALUE SUMMARY ===")
  for (col in c("SST_14d", "SST_30d", "SST_60d")) {
    n_na  <- sum(is.na(cpr[[col]]))
    n_tot <- nrow(cpr)
    message(sprintf("  %s : %d NA / %d (%.1f%%)", col, n_na, n_tot,
                    n_na / n_tot * 100))
  }
  
  message("\n=== STATISTICAL SUMMARY ===")
  print(summary(cpr[, c("SST_14d", "SST_30d", "SST_60d")]))
  
  message("\n=== CROSS-WINDOW CORRELATIONS ===")
  print(round(cor(cpr[, c("SST_14d", "SST_30d", "SST_60d")],
                  use = "complete.obs"), 4))
  
  write_fst(cpr, sst_output_path)
  message(sprintf("File exported: %s", sst_output_path))
}

# ── SST distribution by window ────────────────────────────────────────────────
cpr %>%
  tidyr::pivot_longer(cols = c(SST_14d, SST_30d, SST_60d),
                      names_to  = "fenetre",
                      values_to = "SST") %>%
  ggplot(aes(x = SST, fill = fenetre)) +
  geom_histogram(bins = 80, alpha = 0.6, position = "identity") +
  facet_wrap(~fenetre, ncol = 1) +
  theme_bw(base_size = 13) +
  labs(title = "SST distribution by temporal window",
       x = "SST (degrees C)", y = "Frequency")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5: FAMILY SELECTION AND SST COVARIATE TEST
# ══════════════════════════════════════════════════════════════════════════════

# ── Prepare modelling dataset ─────────────────────────────────────────────────
cpr_sst <- cpr %>%
  dplyr::filter(!is.na(SST_30d)) %>%
  dplyr::mutate(doy             = yday(Midpoint_Date_Local),
                year            = year(Midpoint_Date_Local),
                sum_biomass_log = log(sum_biomass_g + 1))

# ── Fit candidate models ──────────────────────────────────────────────────────

# [A] Tweedie on raw biomass
m_tw <- gam(sum_biomass_g ~
              s(year, bs = "tp") +
              s(doy,  bs = "cc") +
              te(Longitude, Latitude) +
              s(SST_30d, bs = "tp", k = 5),
            data    = cpr_sst,
            family  = tw(),
            control = gam.control(nthreads = 4))

# [B] Gaussian on log(biomass+1), without SST (reference model)
m_ref <- gam(sum_biomass_log ~
               s(year, bs = "tp") +
               s(doy,  bs = "cc") +
               te(Longitude, Latitude),
             data    = cpr_sst,
             family  = gaussian(),
             control = gam.control(nthreads = 4))

# [C] Gaussian on log(biomass+1), SST as a smooth
m_sst_smooth <- gam(sum_biomass_log ~
                      s(year, bs = "tp") +
                      s(doy,  bs = "cc") +
                      te(Longitude, Latitude) +
                      s(SST_30d, bs = "tp", k = 5),
                    data    = cpr_sst,
                    family  = gaussian(),
                    control = gam.control(nthreads = 4))

# [D] Gaussian on log(biomass+1), SST as a linear term (more parsimonious)
m_sst_linear <- gam(sum_biomass_log ~
                      s(year, bs = "tp") +
                      s(doy,  bs = "cc") +
                      te(Longitude, Latitude) +
                      SST_30d,
                    data    = cpr_sst,
                    family  = gaussian(),
                    control = gam.control(nthreads = 4))


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6: DHARMA DIAGNOSTICS — FAMILY SELECTION
# ══════════════════════════════════════════════════════════════════════════════

# DHARMa generates simulated quantile residuals that should be uniformly
# distributed on [0,1] if the model is correctly specified.
#
# Tests interpreted:
#   testUniformity  : KS test — are residuals U[0,1]?
#   testDispersion  : ratio ~1 = dispersion well captured
#                     < 1 = underdispersion; > 1 = overdispersion
#   testZeroInflation: ratio ~1 = zeros well modelled

set.seed(123)
sim_tw    <- simulateResiduals(m_tw,         n = 500)
sim_gauss <- simulateResiduals(m_sst_smooth, n = 500)

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

# Tweedie: dispersion ~1.93 (overdispersed); ZI ratio ~0.93 (zeros OK).
# Gaussian: dispersion ~1.00 (near perfect); ZI = Inf (artefact).
#
# The Inf ZI for the Gaussian is an artefact of the log(y+1) transformation:
# log(0+1) = 0 is the minimum of the continuous transformed response, but a
# Gaussian model predicts values near zero, never exactly zero. DHARMa counts
# exact zeros in simulations, which essentially never occur for a continuous
# distribution, yielding Inf. This reflects a property of the transformation,
# not a model deficiency.
#
# The Tweedie overdispersion (ratio ~1.93) indicates unmodelled variability,
# likely driven by the combined contribution of organism count and individual
# body mass varying with species composition, season, and region. This
# overdispersion is difficult to address without additional covariates.
# The log-Gaussian model is retained.


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7: SST COVARIATE EVALUATION
# ══════════════════════════════════════════════════════════════════════════════

# ── AIC comparison ────────────────────────────────────────────────────────────
# DeltaAIC > 2 = notable improvement; DeltaAIC > 10 = strong improvement.
aic_table <- data.frame(
  Modele = c("Reference (no SST)", "SST smooth", "SST linear"),
  AIC    = c(AIC(m_ref), AIC(m_sst_smooth), AIC(m_sst_linear)),
  DevianceExpliquee = c(
    summary(m_ref)$dev.expl,
    summary(m_sst_smooth)$dev.expl,
    summary(m_sst_linear)$dev.expl) * 100
)
aic_table$DeltaAIC <- aic_table$AIC - min(aic_table$AIC)
aic_table[, -1]    <- round(aic_table[, -1], 2)
aic_table[order(aic_table$AIC), ]

# ── SST ~ Year correlation ─────────────────────────────────────────────────────
# If |r| > 0.3 and p < 0.05, confounding risk — exclude SST.
cor.test(cpr_sst$SST_30d, cpr_sst$year)
# r = 0.188, p < 2.2e-16: weak correlation, no meaningful confounding.

# ── Stability of s(Year) with and without SST ─────────────────────────────────
# Strongly divergent curves would indicate confounding.
par(mfrow = c(1, 2))
plot(m_ref,        select = 1, shade = TRUE,
     main = "s(Year) — without SST", ylim = c(-1, 1))
plot(m_sst_smooth, select = 1, shade = TRUE,
     main = "s(Year) — with SST",    ylim = c(-1, 1))
par(mfrow = c(1, 1))

# ── Temporal residuals: with and without SST ──────────────────────────────────
cpr_sst$resid_ref <- residuals(m_ref,        type = "deviance")
cpr_sst$resid_sst <- residuals(m_sst_smooth, type = "deviance")

resid_comp <- cpr_sst %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(
    resid_ref = mean(resid_ref, na.rm = TRUE),
    resid_sst = mean(resid_sst, na.rm = TRUE)
  )

resid_comp %>%
  tidyr::pivot_longer(cols      = c(resid_ref, resid_sst),
                      names_to  = "modele",
                      values_to = "resid_moy") %>%
  dplyr::mutate(modele = dplyr::recode(modele,
                                       resid_ref = "Without SST",
                                       resid_sst = "With SST (smooth)")) %>%
  ggplot(aes(x = year, y = resid_moy)) +
  geom_point() + geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  facet_wrap(~modele, ncol = 1, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(title = "Temporal residuals — with and without SST",
       y     = "Mean deviance residual",
       x     = "Year")

# No systematic temporal pattern in residuals with or without SST: both curves
# oscillate randomly around zero. The three checks converge: SST improves fit
# (AIC), its correlation with year is too weak to constitute major confounding
# (r = 0.188), and it alters neither the shape of s(Year) nor the annual
# residual pattern. SST is included in the final model.


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 8: TEMPORAL RESIDUAL CHECK — FINAL MODEL
# ══════════════════════════════════════════════════════════════════════════════

cpr_sst$resid_tw    <- residuals(m_tw,         type = "deviance")
cpr_sst$resid_gauss <- residuals(m_sst_smooth, type = "deviance")

resid_temp <- cpr_sst %>%
  dplyr::group_by(year) %>%
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
  ggplot(aes(x = year, y = resid_moy)) +
  geom_point(size = 1.5) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  theme_minimal(base_size = 14) +
  facet_wrap(~modele, ncol = 1, scales = "free_y") +
  labs(title    = "Mean annual residuals — residual temporal pattern",
       subtitle = "Loess near y=0 without trend = valid annual index",
       y        = "Mean deviance residual",
       x        = "Year")

# Near-zero slope and non-significant temporal trend = mean residual bias is
# constant over time (no systematic temporal trend).
lm(resid_tw_moy    ~ year, data = resid_temp) |> summary()
lm(resid_gauss_moy ~ year, data = resid_temp) |> summary()

# The Tweedie systematically overpredicts biomass (t-test significantly != 0).
# The Gaussian residuals are centred on zero (t-test non-significant).
t.test(resid_temp$resid_tw_moy,    mu = 0)
t.test(resid_temp$resid_gauss_moy, mu = 0)

# The Gaussian ZI=Inf is an artefact (see Section 6). What matters is that the
# model predicts low biomass values where observations are zero — biologically
# equivalent at Labrador Sea scale. The log-Gaussian model is retained.


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 9: GAM DIAGNOSTICS — FINAL MODEL
# ══════════════════════════════════════════════════════════════════════════════

par(mfrow = c(2, 2))
gam.check(m_sst_smooth)
par(mfrow = c(1, 1))

# Heteroscedasticity is present: residual variance is highest for low predicted
# values (near-zero biomass observations) and decreases for higher predictions.
# For an annual index, this means low-biomass years are estimated with greater
# uncertainty than high-biomass years — an acceptable limitation at the spatial
# scale of this analysis. To be noted as a limitation in the manuscript.


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 10: FINAL MODEL SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

summary(m_sst_smooth)

par(mfrow = c(2, 2))
plot(m_sst_smooth, select = 1, shade = TRUE,
     main = "Temporal trend — s(Year)",
     ylab = "Partial effect on log(biomass+1)")
plot(m_sst_smooth, select = 2, shade = TRUE,
     main = "Seasonality — s(DOY)",
     ylab = "Partial effect")
plot(m_sst_smooth, select = 3, scheme = 2,
     main = "Spatial surface — te(Lon, Lat)")
plot(m_sst_smooth, select = 4, shade = TRUE,
     main = "SST effect — s(SST)",
     ylab = "Partial effect on log(biomass+1)")
par(mfrow = c(1, 1))


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 11: PREDICTION POINT SELECTION
# ══════════════════════════════════════════════════════════════════════════════

# SOURCE: Labrador Sea polygon used as spatial domain for prediction.
# STATUS: Included in this repository.
lab_sea <- sf::st_transform(
  sf::st_read(file.path(base_path,
                        "R/00_raw_data/ReferenceMigrationLayers/LabradorSea.shp"),
              quiet = TRUE),
  crs = 4326) |> sf::st_make_valid()

world     <- ne_countries(scale = "medium", returnclass = "sf")
provinces <- ne_states(country = "Canada", returnclass = "sf")
states    <- sf::st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))

# Convert CPR data to sf and retain only points within the Labrador Sea polygon
cpr_ordered <- cpr[order(cpr$Midpoint_Date_Local), ]
cpr_sf <- sf::st_as_sf(
  cpr_ordered,
  coords = c("Longitude", "Latitude"),
  crs    = 4326,
  remove = FALSE
) |> sf::st_make_valid()

intersections <- sf::st_intersects(cpr_sf, lab_sea, sparse = FALSE)
cpr_sf        <- cpr_sf[which(intersections), ]

# Sample 100 random prediction points from within CPR coverage in Labrador Sea
set.seed(10)
ss    <- sample(1:nrow(cpr_sf), n_points)
ss_sf <- cpr_sf[ss, ]

land_mask     <- sf::st_union(world) |> sf::st_make_valid()
lab_sea_water <- sf::st_difference(lab_sea, land_mask) |> sf::st_make_valid()

# ── Map of prediction points ───────────────────────────────────────────────────
plot_pred_points <- ggplot() +
  geom_sf(data = sf::st_as_sfc(sf::st_bbox(
    c(xmin = -80, xmax = -30, ymin = 40, ymax = 70), crs = 4326)),
    fill = "#D6E8F5", color = NA) +
  geom_sf(data = world,     fill = "gray89", color = "#AAAAAA", linewidth = 0.25) +
  geom_sf(data = provinces, fill = "gray89", color = "#AAAAAA", linewidth = 0.15) +
  geom_sf(data = states,    fill = NA,       color = "#CCCCCC", linewidth = 0.15) +
  geom_sf(data = lab_sea_water,
          fill      = "#4B9FD4", alpha    = 0.20,
          color     = "#2C6E9E", linewidth = 0.3, linetype = "dashed") +
  geom_sf(data = cpr_sf, color = "red",   size  = 1,   alpha  = 0.5) +
  geom_sf(data = ss_sf,  shape = 8,       color = "black", size = 3, stroke = 1) +
  annotate("text", x = -51.5, y = 59.0,
           label    = "Labrador Sea", fontface = "bold",
           size     = 3.8, color = "#1a4e7a") +
  ggspatial::annotation_north_arrow(
    location    = "tl", which_north = "true",
    pad_x       = unit(0.4, "cm"), pad_y = unit(0.4, "cm"),
    height      = unit(1.8, "cm"), width = unit(1.8, "cm"),
    style       = ggspatial::north_arrow_nautical(
      fill      = c("black", "white"),
      line_col  = "grey20", text_size = 9)) +
  ggspatial::annotation_scale(
    location    = "bl", width_hint = 0.18,
    text_cex    = 0.70,
    line_col    = "#444444", text_col   = "#444444",
    bar_cols    = c("black", "white"),
    pad_x       = unit(0.4, "cm"), pad_y = unit(0.4, "cm")) +
  coord_sf(xlim = c(-71, -38), ylim = c(44, 62), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme(
    panel.background  = element_rect(fill = "#D6E8F5", colour = "white"),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major  = element_line(color = "#DDDDDD", linewidth = 0.2),
    axis.text         = element_text(size = 10, color = "black"),
    axis.title        = element_text(size = 12, color = "black"),
    axis.ticks        = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    plot.margin       = margin(2, 6, 2, 2)
  )
plot_pred_points

ss_df <- as.data.frame(ss_sf)
ss_df <- ss_df[, -which(colnames(ss_df) == "geometry")]


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 12: SST EXTRACTION FOR 100 PREDICTION POINTS
# ══════════════════════════════════════════════════════════════════════════════

# For each combination of point (100) x day (~15,500 days over ~42 years),
# compute the 30-day rolling mean SST ending on that day.
# Results are checkpointed monthly to allow resumption after interruption.

dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

if (file.exists(final_sst_path)) {
  
  message("Prediction SST already available, loading without recalculation.")
  sst_prediction <- readRDS(final_sst_path)
  message(sprintf("  Matrix loaded: %d points x %d days.",
                  nrow(sst_prediction$sst_matrix),
                  ncol(sst_prediction$sst_matrix)))
  
} else {
  
  # ── Build daily NC file index ──────────────────────────────────────────────
  
  series_start  <- as.Date(sprintf("%d-01-01",
                                   min(lubridate::year(cpr$Midpoint_Date_Local))))
  series_end    <- as.Date(sprintf("%d-12-31",
                                   max(lubridate::year(cpr$Midpoint_Date_Local))))
  earliest_need <- series_start - window_days
  
  message(sprintf("Target series  : %s to %s", series_start, series_end))
  message(sprintf("SST required   : %s to %s (30d buffer included)",
                  earliest_need, series_end))
  
  zip_files <- list.files(sst_dir,
                          pattern    = "ESACCI_SST_L4_\\d{4}_\\d{2}\\.nc$",
                          full.names = TRUE)
  message(sprintf("%d SST ZIP files found.", length(zip_files)))
  
  file_index <- lapply(zip_files, function(f) {
    bname <- basename(f)
    yr_k  <- as.integer(sub("ESACCI_SST_L4_(\\d{4})_(\\d{2})\\.nc", "\\1", bname))
    mo_k  <- as.integer(sub("ESACCI_SST_L4_(\\d{4})_(\\d{2})\\.nc", "\\2", bname))
    if (is.na(yr_k) || is.na(mo_k)) return(NULL)
    
    month_start <- as.Date(sprintf("%04d-%02d-01", yr_k, mo_k))
    month_end   <- as.Date(sprintf("%04d-%02d-%02d", yr_k, mo_k,
                                   lubridate::days_in_month(month_start)))
    if (month_end < earliest_need || month_start > series_end) return(NULL)
    
    contents <- tryCatch(unzip(f, list = TRUE), error = function(e) NULL)
    if (is.null(contents) || nrow(contents) == 0L) return(NULL)
    
    nc_names  <- contents$Name
    date_strs <- substr(basename(nc_names), 1, 8)
    dates_i   <- as.Date(date_strs, format = "%Y%m%d")
    valid     <- !is.na(dates_i) & dates_i >= earliest_need & dates_i <= series_end
    
    if (!any(valid)) return(NULL)
    
    data.frame(
      date      = dates_i[valid],
      zip_path  = f,
      nc_name   = nc_names[valid],
      yr        = yr_k,
      mo        = mo_k,
      stringsAsFactors = FALSE
    )
  })
  
  file_index <- do.call(rbind, file_index[!sapply(file_index, is.null)])
  file_index <- file_index[order(file_index$date), ]
  rownames(file_index) <- NULL
  
  months_index <- unique(file_index[, c("yr", "mo")])
  months_index <- months_index[order(months_index$yr, months_index$mo), ]
  
  message(sprintf("SST index: %d days across %d months.",
                  nrow(file_index), nrow(months_index)))
  
  # ── Spatial buffers for 100 prediction points ──────────────────────────────
  pts_vect <- terra::vect(
    data.frame(lon = ss_df$Longitude, lat = ss_df$Latitude),
    geom = c("lon", "lat"),
    crs  = "EPSG:4326"
  )
  pts_m    <- terra::project(pts_vect, "EPSG:32621")
  bufs_m   <- terra::buffer(pts_m, width = radius_m)
  bufs_wgs <- terra::project(bufs_m, "EPSG:4326")
  ext_pts  <- terra::ext(bufs_wgs)
  
  bufs_wgs_file <- tempfile(fileext = ".rds")
  saveRDS(bufs_wgs, bufs_wgs_file)
  
  message(sprintf("Buffers created for %d points (radius = %d km).",
                  nrow(ss_df), radius_m / 1000))
  
  # ── Accumulation matrix ────────────────────────────────────────────────────
  # sst_daily [n_pts x n_days]: daily SST extracted at each point,
  # including the 30d pre-series buffer needed for rolling window computation.
  
  n_pts      <- nrow(ss_df)
  all_dates  <- seq.Date(earliest_need, series_end, by = "day")
  n_days_tot <- length(all_dates)
  
  sst_daily <- matrix(NA_real_, nrow = n_pts, ncol = n_days_tot)
  colnames(sst_daily) <- as.character(all_dates)
  
  message(sprintf("SST matrix: %d points x %d days (%.1f MB estimated).",
                  n_pts, n_days_tot,
                  n_pts * n_days_tot * 8 / 1e6))
  
  # ── Resume from last checkpoint if available ───────────────────────────────
  last_checkpoint <- file.path(checkpoint_dir, "last_completed_month.rds")
  
  if (file.exists(last_checkpoint)) {
    ckpt      <- readRDS(last_checkpoint)
    sst_daily <- ckpt$sst_daily
    last_yr   <- ckpt$last_yr
    last_mo   <- ckpt$last_mo
    
    months_todo <- months_index[
      months_index$yr > last_yr |
        (months_index$yr == last_yr & months_index$mo > last_mo), ]
    
    message(sprintf("Resuming from checkpoint: last completed month = %04d-%02d.",
                    last_yr, last_mo))
    message(sprintf("Remaining months: %d / %d.",
                    nrow(months_todo), nrow(months_index)))
  } else {
    months_todo <- months_index
    message(sprintf("Starting from beginning: %d months to process.",
                    nrow(months_todo)))
  }
  
  # ── Main extraction loop: by month ────────────────────────────────────────
  # For each month, process daily NC files. Extraction across the 100 points
  # is parallelised within each day. Checkpoint saved after each month.
  
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterEvalQ(cl, library(terra))
  parallel::clusterExport(cl, c("bufs_wgs_file", "ext_pts"), envir = environment())
  parallel::clusterEvalQ(cl, {
    bufs_wgs_worker <- readRDS(bufs_wgs_file)
  })
  
  t_start <- proc.time()
  
  for (mi in seq_len(nrow(months_todo))) {
    
    yr_m <- months_todo$yr[mi]
    mo_m <- months_todo$mo[mi]
    
    fi_month <- file_index[file_index$yr == yr_m & file_index$mo == mo_m, ]
    if (nrow(fi_month) == 0L) next
    
    for (di in seq_len(nrow(fi_month))) {
      
      day_date  <- fi_month$date[di]
      date_idx  <- which(all_dates == day_date)
      if (length(date_idx) == 0L) next
      
      zip_path_d <- fi_month$zip_path[di]
      nc_name_d  <- fi_month$nc_name[di]
      
      tmp_dir <- tempfile()
      dir.create(tmp_dir)
      
      r_sst <- tryCatch({
        extracted <- unzip(zip_path_d, files = nc_name_d, exdir = tmp_dir)
        if (length(extracted) == 0L) stop("unzip empty")
        
        r_all   <- terra::rast(extracted[1])
        sst_idx <- grep("analysed_sst$", names(r_all), value = FALSE)
        if (length(sst_idx) == 0L) sst_idx <- 1L
        
        r_k <- r_all[[sst_idx[1]]] - 273.15
        terra::crs(r_k) <- "EPSG:4326"
        terra::crop(r_k, ext_pts)
        
      }, error = function(e) {
        message(sprintf("  [SKIP] %s/%s : %s",
                        basename(zip_path_d), basename(nc_name_d),
                        conditionMessage(e)))
        NULL
      }, finally = {
        unlink(tmp_dir, recursive = TRUE)
      })
      
      if (is.null(r_sst)) next
      
      r_tmp_file <- tempfile(fileext = ".tif")
      terra::writeRaster(r_sst, r_tmp_file, overwrite = TRUE)
      parallel::clusterExport(cl, "r_tmp_file", envir = environment())
      
      chunks <- parallel::splitIndices(n_pts, n_cores)
      
      extracted_vals <- unlist(
        parallel::parLapply(cl, chunks, function(idx) {
          r_worker <- terra::rast(r_tmp_file)
          vals <- tryCatch(
            terra::extract(r_worker, bufs_wgs_worker[idx, ],
                           fun = mean, na.rm = TRUE)[, 2],
            error = function(e) rep(NA_real_, length(idx))
          )
          vals
        })
      )
      
      unlink(r_tmp_file)
      
      sst_daily[, date_idx] <- extracted_vals
    }
    
    saveRDS(
      list(sst_daily = sst_daily,
           last_yr   = yr_m,
           last_mo   = mo_m),
      last_checkpoint
    )
    
    elapsed <- (proc.time() - t_start)["elapsed"]
    message(sprintf("  [%d/%d] %04d-%02d completed — %.0fs elapsed (%.1f min).",
                    mi, nrow(months_todo), yr_m, mo_m,
                    elapsed, elapsed / 60))
  }
  
  parallel::stopCluster(cl)
  message("Daily extraction complete.")
  
  # ── Compute 30-day rolling mean SST ───────────────────────────────────────
  # For each point and each target day (series_start to series_end),
  # SST_30d = mean of the 30 days ending on that day (inclusive).
  # Buffer days (before series_start) are used in calculations but excluded
  # from the output matrix.
  
  message("Computing 30-day rolling mean SST...")
  
  target_idx   <- which(all_dates >= series_start)
  target_dates <- all_dates[target_idx]
  n_target     <- length(target_dates)
  
  sst_30d <- matrix(NA_real_, nrow = n_pts, ncol = n_target)
  colnames(sst_30d) <- as.character(target_dates)
  
  for (di in seq_len(n_target)) {
    end_idx     <- target_idx[di]
    start_idx   <- max(1L, end_idx - window_days + 1L)
    window_vals <- sst_daily[, start_idx:end_idx, drop = FALSE]
    sst_30d[, di] <- rowMeans(window_vals, na.rm = TRUE)
  }
  
  sst_30d[is.nan(sst_30d)] <- NA_real_
  
  message(sprintf("Final SST_30d matrix: %d points x %d days.", n_pts, n_target))
  message(sprintf("Residual NAs: %.1f%%", mean(is.na(sst_30d)) * 100))
  
  # ── Export ─────────────────────────────────────────────────────────────────
  sst_prediction <- list(
    sst_matrix   = sst_30d,
    dates        = target_dates,
    pts          = ss_df[, c("Longitude", "Latitude")],
    series_start = series_start,
    series_end   = series_end,
    window_days  = window_days,
    radius_m     = radius_m
  )
  
  saveRDS(sst_prediction, final_sst_path)
  message(sprintf("Exported: %s", final_sst_path))
  message("Intermediate checkpoints retained in: ", checkpoint_dir)
}

# ── SST distribution check ────────────────────────────────────────────────────
sst_medians <- apply(sst_prediction$sst_matrix, 1, median, na.rm = TRUE)
hist(sst_medians, breaks = 30,
     main = "Median SST distribution across prediction points",
     xlab = "Median SST_30d (degrees C)")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 13: ANNUAL BIOMASS PREDICTION — 100 FIXED POINTS
# ══════════════════════════════════════════════════════════════════════════════

# For each fixed point, predict daily biomass over the full time series,
# aggregate to monthly means, then to annual sums, and compute the median
# across the 100 points as the annual index.
#
# The model predicts log(biomass+1); back-transform via exp(pred) - 1.

all_dates_pred <- sst_prediction$dates

newdata_list <- lapply(seq_len(nrow(ss_df)), function(i) {
  data.frame(
    Midpoint_Date_Local = all_dates_pred,
    year                = lubridate::year(all_dates_pred),
    doy                 = lubridate::yday(all_dates_pred),
    Longitude           = ss_df$Longitude[i],
    Latitude            = ss_df$Latitude[i],
    SST_30d             = sst_prediction$sst_matrix[i, ]
  ) %>%
    dplyr::filter(!is.na(SST_30d))
})

all_newdata <- dplyr::bind_rows(
  lapply(seq_along(newdata_list), function(i) {
    newdata_list[[i]] %>% dplyr::mutate(it = i)
  })
)

pred_all <- predict(m_sst_smooth,
                    newdata = all_newdata,
                    type    = "link",
                    se.fit  = TRUE)

plankton <- all_newdata %>%
  dplyr::mutate(
    fit    = exp(pred_all$fit) - 1,
    fit_lo = exp(pred_all$fit - 1.96 * pred_all$se.fit) - 1,
    fit_hi = exp(pred_all$fit + 1.96 * pred_all$se.fit) - 1,
    month  = lubridate::month(Midpoint_Date_Local)
  ) %>%
  dplyr::group_by(it, year, month) %>%
  dplyr::summarise(
    Latitude       = first(Latitude),
    Longitude      = first(Longitude),
    fit_monthly    = mean(fit,    na.rm = TRUE),
    fit_lo_monthly = mean(fit_lo, na.rm = TRUE),
    fit_hi_monthly = mean(fit_hi, na.rm = TRUE),
    .groups        = "drop"
  ) %>%
  dplyr::group_by(it, year) %>%
  dplyr::summarise(
    Latitude       = first(Latitude),
    Longitude      = first(Longitude),
    Plankton       = sum(fit_monthly,    na.rm = TRUE),
    Plankton_lo    = sum(fit_lo_monthly, na.rm = TRUE),
    Plankton_hi    = sum(fit_hi_monthly, na.rm = TRUE),
    n_months_valid = sum(!is.na(fit_monthly)),
    .groups        = "drop"
  )

plankton <- dplyr::bind_rows(plankton) %>%
  dplyr::filter(n_months_valid == 12)


# ── Annual median index across 100 points ─────────────────────────────────────
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
gap_xmin  <- min(gap_years) - 0.5
gap_xmax  <- max(gap_years) + 0.5

ggplot(data = plankton, aes(x = year, y = Plankton, group = it)) +
  geom_line(linewidth = 0.5, alpha = 0.4) +
  annotate("rect", xmin = gap_xmin, xmax = gap_xmax,
           ymin = -Inf, ymax = Inf, fill = "grey40", alpha = 0.3) +
  theme_bw(base_size = 14) +
  labs(x = "Year", y = "Plankton biomass (g/sample)")

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
  annotate("rect", xmin = gap_xmin, xmax = gap_xmax,
           ymin = -Inf, ymax = Inf, fill = "grey40", alpha = 0.3) +
  theme_bw(base_size = 14) +
  labs(x        = "Year",
       y        = "Plankton biomass (g/sample)",
       title    = "Annual biomass index — Labrador Sea",
       subtitle = "Red: median of 100 points | Ribbon: IQR")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 14: TEMPORAL LOO CROSS-VALIDATION — ANNUAL MEDIAN INDEX
# ══════════════════════════════════════════════════════════════════════════════

# For each year removed, refit the GAM and compare the predicted annual median
# index (100 points) against the full-model reference index.

pseudo_gap_years <- sort(unique(cpr_sst$year))

ref_index <- summ_plankton %>%
  dplyr::select(year, zoo_median_ref = zoo_median)

set.seed(42)
loo_index <- lapply(pseudo_gap_years, function(yr) {
  cat("  Annual index for", yr, "...\n")
  
  train <- cpr_sst %>% dplyr::filter(year != yr)
  
  m_cv <- tryCatch(
    gam(sum_biomass_log ~
          s(year,    bs = "tp", k = 5) +
          s(doy,     bs = "cc") +
          te(Longitude, Latitude, k = c(15, 15)) +
          s(SST_30d, bs = "tp", k = 5),
        data    = train,
        family  = gaussian(),
        control = gam.control(nthreads = 4)),
    error = function(e) NULL
  )
  if (is.null(m_cv)) return(NULL)
  
  plankton_cv <- lapply(seq_len(nrow(ss_df)), function(i) {
    nd <- data.frame(
      Midpoint_Date_Local = sst_prediction$dates,
      year                = lubridate::year(sst_prediction$dates),
      doy                 = lubridate::yday(sst_prediction$dates),
      Longitude           = ss_df$Longitude[i],
      Latitude            = ss_df$Latitude[i],
      SST_30d             = sst_prediction$sst_matrix[i, ]
    ) %>%
      dplyr::filter(year == yr, !is.na(SST_30d))
    
    if (nrow(nd) == 0L) return(NULL)
    
    pred <- tryCatch(
      predict(m_cv, newdata = nd, type = "link", se.fit = FALSE),
      error = function(e) NULL
    )
    if (is.null(pred)) return(NULL)
    
    nd %>%
      dplyr::mutate(
        fit   = pmax(exp(pred) - 1, 0),
        month = lubridate::month(Midpoint_Date_Local)
      ) %>%
      dplyr::group_by(month) %>%
      dplyr::summarise(fit_monthly = mean(fit, na.rm = TRUE), .groups = "drop") %>%
      dplyr::summarise(Plankton = sum(fit_monthly, na.rm = TRUE)) %>%
      dplyr::mutate(it = i, year = yr)
  })
  
  plankton_cv <- dplyr::bind_rows(plankton_cv)
  if (nrow(plankton_cv) == 0L) return(NULL)
  
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
  annotate("text",
           x     = min(loo_index$index_ref) +
             0.01 * diff(range(loo_index$index_ref)),
           y     = max(loo_index$index_cv)  -
             0.01 * diff(range(loo_index$index_cv)),
           label = paste0("R2 = ",   round(metrics_index$R2,   3),
                          "\nMAPE = ", round(metrics_index$MAPE, 1), "%"),
           hjust = 0, vjust = 1, size = 4, family = "mono") +
  theme_bw(base_size = 13) +
  labs(title    = "Temporal LOO — Annual median index",
       subtitle = "Predicted index (year removed) vs full-model reference index",
       x        = "Reference index (g/sample per day)",
       y        = "LOO predicted index (g/sample per day)")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 15: EXPORT
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





