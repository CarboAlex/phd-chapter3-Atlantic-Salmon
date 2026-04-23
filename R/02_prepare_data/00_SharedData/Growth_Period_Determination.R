

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: HEADER
# ══════════════════════════════════════════════════════════════════════════════
#
# Identifies the annual growth period for Atlantic salmon in the St. Jean
# and Trinite rivers by fitting a cyclic GAM to daily water temperatures and
# detecting when the smoothed curve crosses the 6 degree C threshold.
#
# The 6 degree C growth threshold is based on:
#   Elliott, J.M. & Hurley, M.A. (1997). A functional model for maximum
#     growth of Atlantic salmon parr, Salmo salar, from two populations in
#     northwest England. Functional Ecology, 11: 592-603.
#   Elliott, J.M. & Elliott, J.A. (2010). Temperature requirements of
#     Atlantic salmon Salmo salar, brown trout Salmo trutta and Arctic charr
#     Salvelinus alpinus: predicting the effects of climate change. Journal of
#     Fish Biology, 77: 1793-1817.
#   Ouellet-Proulx, S. et al. (2023). A potential growth thermal index for
#     estimating juvenile Atlantic salmon (Salmo salar) size-at-age across
#     geographical scales. Journal of Fish Biology,
#     doi: 10.1111/jfb.15535.
#   Hendry, K. & Cragg-Hine, D. (2003). Ecology of the Atlantic salmon.
#     Conserving Natura 2000 Rivers Ecology Series No. 7. English Nature,
#     Peterborough.
#   Symons, P.E.K. (1979). Estimated escapement of Atlantic salmon
#     (Salmo salar) for maximum smolt production in rivers of different
#     productivity. Journal of the Fisheries Research Board of Canada,
#     36: 132-140.
#
# Inputs:
#   - R/00_raw_data/Temperature_Discharge_River/
#       CEQEAU_Temperature_Discharge_1979_2024_SJ.csv
#       CEQEAU_Temperature_Discharge_1979_2024_TR.csv
#       Simulated daily water temperature for Saint-Jean and Trinite (CSV)
#
# Outputs:
#   - R/01_derived_data/GrowthPeriods/Growth_Periods_STJ.fst
#       Annual growth period (start/end DOY) for Saint-Jean (FST)
#   - R/01_derived_data/GrowthPeriods/Growth_Periods_TRI.fst
#       Annual growth period (start/end DOY) for Trinite (FST)
#
# ══════════════════════════════════════════════════════════════════════════════


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: LIBRARIES
# ══════════════════════════════════════════════════════════════════════════════

library(fst)
library(dplyr)
library(lubridate)
library(ggplot2)
library(mgcv)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3: WORKING DIRECTORY
# ══════════════════════════════════════════════════════════════════════════════

# ── Working directory ──────────────────────────────────────────────────────
# Automatically detects the repository root based on the script location.
# No hardcoded path — works on any machine.
current_file <- rstudioapi::getActiveDocumentContext()$path
path_parts   <- unlist(strsplit(normalizePath(current_file, winslash = "/"), "/"))
target_index <- which(path_parts == "phd-chapter3-Atlantic-Salmon")
base_path    <- paste(path_parts[1:target_index], collapse = "/")
setwd(base_path)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4: DATA IMPORT
# ══════════════════════════════════════════════════════════════════════════════

# ── Daily water temperature — Saint-Jean ──────────────────────────────────────
# SOURCE: CEQEAU hydrological model outputs
# STATUS: Confidential — not included in this repository.
#         An example file with identical structure and simulated data is
#         available at the same path with the suffix _EXAMPLE.csv.
#         To request access, contact the data provider once the thesis
#         is published.
temp_stj <- read.csv(
  file.path(base_path,
            "R/00_raw_data/Temperature_Discharge_River",
            "CEQEAU_Temperature_Discharge_1979_2024_SJ.csv"))


# ── Daily water temperature — Trinite ────────────────────────────────────────
# SOURCE: CEQEAU hydrological model outputs
# STATUS: Confidential — not included in this repository.
#         An example file with identical structure and simulated data is
#         available at the same path with the suffix _EXAMPLE.csv.
#         To request access, contact the data provider once the thesis
#         is published.
temp_tri <- read.csv(
  file.path(base_path,
            "R/00_raw_data/Temperature_Discharge_River",
            "CEQEAU_Temperature_Discharge_1979_2024_TR.csv"))


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5: CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

temp_threshold <- 6  # Growth temperature threshold (degrees C); see references
# in the header for justification.


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6: FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

# ── plot_gam_year ─────────────────────────────────────────────────────────────
# Purpose: Fits a cyclic GAM for a single year at a given k and plots the
#   smoothed curve over raw data, with the 6 °C threshold and detected
#   growth period boundaries. Used to visually select an appropriate k.
# Parameters:
#   data            data frame with columns year, doy, temperature_c
#   year_plot       integer; year to plot
#   temp_threshold  numeric; temperature threshold in degrees C (default: 6)
#   k               integer; number of basis functions for the cyclic spline
#   river_name      character; river name used in the plot title
# Returns: a ggplot object.

plot_gam_year <- function(data, year_plot, temp_threshold = 6, k = 20,
                          river_name = "River") {
  
  data_i <- data %>% filter(year == year_plot)
  
  gam_i <- gam(
    temperature_c ~ s(doy, bs = "cc", k = k),
    data   = data_i,
    method = "REML",
    knots  = list(doy = c(0.5, 366.5))
  )
  
  newdata_i      <- data.frame(doy = 1:366)
  pred           <- predict(gam_i, newdata_i, se.fit = TRUE)
  newdata_i$pred <- pred$fit
  newdata_i$lwr  <- pred$fit - 1.96 * pred$se.fit
  newdata_i$upr  <- pred$fit + 1.96 * pred$se.fit
  newdata_i$above <- newdata_i$pred >= temp_threshold
  
  above      <- newdata_i$above
  cross_up   <- above & !c(FALSE, head(above, -1))
  cross_down <- !above & c(FALSE, head(above, -1))
  start_doy  <- newdata_i$doy[min(which(cross_up))]
  end_doy    <- newdata_i$doy[max(which(cross_down))]
  
  dev_expl <- round(summary(gam_i)$dev.expl * 100, 1)
  
  ribbon_data <- newdata_i %>% filter(doy >= start_doy & doy <= end_doy)
  
  ggplot() +
    geom_ribbon(data = ribbon_data,
                aes(x = doy, ymin = -Inf, ymax = Inf),
                fill = "#4CAF50", alpha = 0.08) +
    geom_ribbon(data = newdata_i,
                aes(x = doy, ymin = lwr, ymax = upr),
                fill = "#1565C0", alpha = 0.15) +
    geom_point(data = data_i,
               aes(x = doy, y = temperature_c),
               size = 0.8, alpha = 0.4, color = "grey40") +
    geom_line(data = newdata_i,
              aes(x = doy, y = pred),
              color = "#1565C0", linewidth = 1) +
    geom_hline(yintercept = temp_threshold,
               linetype = "dashed", color = "#E53935", linewidth = 0.7) +
    geom_vline(xintercept = start_doy,
               linetype = "dotted", color = "#4CAF50", linewidth = 0.8) +
    geom_vline(xintercept = end_doy,
               linetype = "dotted", color = "#4CAF50", linewidth = 0.8) +
    annotate("text", x = start_doy + 3,
             y = max(data_i$temperature_c, na.rm = TRUE),
             label = paste0("Start: DOY ", start_doy),
             hjust = 0, size = 3.2, color = "#4CAF50") +
    annotate("text", x = end_doy - 3,
             y = max(data_i$temperature_c, na.rm = TRUE),
             label = paste0("End: DOY ", end_doy),
             hjust = 1, size = 3.2, color = "#4CAF50") +
    annotate("text", x = 5,
             y = max(data_i$temperature_c, na.rm = TRUE) - 1,
             label = paste0("Deviance explained: ", dev_expl, "%"),
             hjust = 0, size = 3, color = "grey30", fontface = "italic") +
    labs(
      title    = paste0(river_name, " — ", year_plot, " (k = ", k, ")"),
      subtitle = "Cyclic GAM (bs = 'cc') — 6 °C threshold in red",
      x        = "Day of year (DOY)",
      y        = "Water temperature (°C)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey40", size = 10)
    )
}

# ── yearly_growth_period ──────────────────────────────────────────────────────
# Purpose: Fits a cyclic GAM to daily water temperature for each year and
#   returns the first and last DOY at which the smoothed curve crosses
#   temp_threshold, defining the annual growth period.
# Parameters:
#   data            data frame with columns year, doy, temperature_c
#   temp_threshold  numeric; temperature threshold in degrees C (default: 6)
# Returns: data frame with columns year, start_doy, end_doy, start_date,
#   end_date; one row per year with a valid crossing.

yearly_growth_period <- function(data, temp_threshold = 6, k_value = 30) {
  
  years_all    <- unique(data$year)
  results_list <- list()
  
  for (i in seq_along(years_all)) {
    
    data_i <- data %>%
      dplyr::filter(year == years_all[i])
    
    # Cyclic cubic spline ensures smooth continuity across the year boundary
    gam_i <- gam(
      temperature_c ~ s(doy, bs = "cc", k = k_value),
      data   = data_i,
      method = "REML",
      knots  = list(doy = c(0.5, 366.5))
    )
    
    # Retrieve the goodness-of-fit metrics from the GAM
    gam_summary_i <- summary(gam_i)
    dev_expl_i    <- gam_summary_i$dev.expl # An Explanation of Deviance (0–1)
    reml_score_i  <- gam_i$gcv.ubre         # REML Score
    
    newdata_i          <- data.frame(year = years_all[i], doy = 1:366)
    newdata_i$pred_temp <- predict(gam_i, newdata_i)
    
    # Identify upward and downward crossings of the temperature threshold
    newdata_i <- newdata_i %>%
      mutate(
        above      = pred_temp >= temp_threshold,
        cross_up   = above & !lag(above, default = FALSE),
        cross_down = !above & lag(above, default = FALSE)
      )
    
    first_day <- newdata_i$doy[min(which(newdata_i$cross_up))]
    last_day  <- newdata_i$doy[max(which(newdata_i$cross_down))]
    
    if (is.infinite(first_day) | is.infinite(last_day)) {
      periods_above_threshold <- data.frame()
    } else {
      periods_above_threshold <- data.frame(
        year      = years_all[i],
        start_doy = first_day,
        end_doy   = last_day,
        dev_expl  = round(dev_expl_i, 3),
        reml      = round(reml_score_i, 3)
      ) %>%
        mutate(
          start_date = as.Date(start_doy - 1, origin = paste0(year, "-01-01")),
          end_date   = as.Date(end_doy   - 1, origin = paste0(year, "-01-01"))
        )
    }
    
    results_list[[i]] <- periods_above_threshold
  }
  
  bind_rows(results_list)
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7: MAIN SCRIPT
# ══════════════════════════════════════════════════════════════════════════════

# ── Reformat temperature data ─────────────────────────────────────────────────

# Standardize column names across both rivers before further processing
colnames(temp_stj)[
  which(colnames(temp_stj) == "Observed.Discharge..m3.s.")] <-
  "Observed_Discharge_m3s"
colnames(temp_stj)[
  which(colnames(temp_stj) == "Simulated.Discharge..m3.s.")] <-
  "Simulated_Discharge_m3s"
colnames(temp_stj)[
  which(colnames(temp_stj) == "Tw..degree.C.")] <-
  "Temperature_C"

colnames(temp_tri)[
  which(colnames(temp_tri) == "Observed.Discharge..RatingCurve...m3.s.")] <-
  "Observed_Discharge_m3s"
colnames(temp_tri)[
  which(colnames(temp_tri) == "Simulated.Discharge..m3.s.")] <-
  "Simulated_Discharge_m3s"
colnames(temp_tri)[
  which(colnames(temp_tri) == "Simulated.Tw..degree.C.")] <-
  "Temperature_C"

colnames(temp_stj) <- tolower(colnames(temp_stj))
colnames(temp_tri) <- tolower(colnames(temp_tri))

# Parse date strings — temporarily set locale to C to ensure English month
# abbreviations are recognised regardless of the user's system locale
old_locale <- Sys.getlocale("LC_TIME")
Sys.setlocale("LC_TIME", "C")
temp_stj$date <- as.Date(temp_stj$date, format = "%d-%b-%Y")
temp_tri$date <- as.Date(temp_tri$date, format = "%d-%b-%Y")
Sys.setlocale("LC_TIME", old_locale)

# Keep only the temperature column needed for GAM fitting
temp_stj <- subset(temp_stj, select = c("date", "temperature_c"))
temp_tri <- subset(temp_tri, select = c("date", "temperature_c"))

# Derive temporal columns used for year-by-year filtering inside the function
temp_stj$year  <- year(temp_stj$date)
temp_stj$month <- month(temp_stj$date)
temp_stj$day   <- day(temp_stj$date)
temp_stj$doy   <- yday(temp_stj$date)

temp_tri$year  <- year(temp_tri$date)
temp_tri$month <- month(temp_tri$date)
temp_tri$day   <- day(temp_tri$date)
temp_tri$doy   <- yday(temp_tri$date)


# ── Visual selection of k ─────────────────────────────────────────────────────
# Compare GAM smooths for several k values on representative years before
# committing to a final k in yearly_growth_period().
# Adjust year_plot and k values as needed, then comment out this block.

set.seed(125)
years_to_check <- sample(unique(temp_stj$year), 5)
k_to_check     <- c(10, 15, 20, 30, 40, 50)

# St. jean
plots_stj <- lapply(years_to_check, function(yr) {
  lapply(k_to_check, function(k_val) {
    plot_gam_year(temp_stj, year_plot = yr, k = k_val, river_name = "St. Jean")
  })
})

# Trinite
plots_tri <- lapply(years_to_check, function(yr) {
  lapply(k_to_check, function(k_val) {
    plot_gam_year(temp_tri, year_plot = yr, k = k_val, river_name = "Trinite")
  })
})

# Display one year at a time across k values (St. Jean)
patchwork::wrap_plots(plots_stj[[1]], ncol = 2)
patchwork::wrap_plots(plots_stj[[2]], ncol = 2)
patchwork::wrap_plots(plots_stj[[3]], ncol = 2)
patchwork::wrap_plots(plots_stj[[4]], ncol = 2)
patchwork::wrap_plots(plots_stj[[5]], ncol = 2)

# Display one year at a time across k values (Trinite)
patchwork::wrap_plots(plots_tri[[1]], ncol = 2)
patchwork::wrap_plots(plots_tri[[2]], ncol = 2)
patchwork::wrap_plots(plots_tri[[3]], ncol = 2)
patchwork::wrap_plots(plots_tri[[4]], ncol = 2)
patchwork::wrap_plots(plots_tri[[5]], ncol = 2)


# ── Detect annual growth periods ──────────────────────────────────────────────

range_stj <- yearly_growth_period(temp_stj, temp_threshold = temp_threshold, k_value = 30)
range_tri <- yearly_growth_period(temp_tri, temp_threshold = temp_threshold, k_value = 30)


# ── Summary of growth periods and quality of GAMs ────────────────────

summary_stj <- range_stj %>%
  summarise(
    mean_start_doy  = round(mean(start_doy), 1),
    mean_end_doy    = round(mean(end_doy), 1),
    mean_start_date = format(as.Date(mean(start_doy) - 1, origin = "2001-01-01"), "%d %b"),
    mean_end_date   = format(as.Date(mean(end_doy)   - 1, origin = "2001-01-01"), "%d %b"),
    mean_dev_expl   = round(mean(dev_expl, na.rm = TRUE), 3),
    min_dev_expl    = round(min(dev_expl,  na.rm = TRUE), 3)
  )

summary_tri <- range_tri %>%
  summarise(
    mean_start_doy  = round(mean(start_doy), 1),
    mean_end_doy    = round(mean(end_doy), 1),
    mean_start_date = format(as.Date(mean(start_doy) - 1, origin = "2001-01-01"), "%d %b"),
    mean_end_date   = format(as.Date(mean(end_doy)   - 1, origin = "2001-01-01"), "%d %b"),
    mean_dev_expl   = round(mean(dev_expl, na.rm = TRUE), 3),
    min_dev_expl    = round(min(dev_expl,  na.rm = TRUE), 3)
  )

summary_stj
summary_tri


# ── Export ────────────────────────────────────────────────────────────────────

write_fst(x    = range_stj,
          path = file.path(base_path,
                           "R/01_derived_data",
                           "Growth_Periods_StJean.fst"))

write_fst(x    = range_tri,
          path = file.path(base_path,
                           "R/01_derived_data",
                           "Growth_Periods_Trinite.fst"))






