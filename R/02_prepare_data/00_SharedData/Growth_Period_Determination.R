

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: HEADER
# ══════════════════════════════════════════════════════════════════════════════
#
# Identifies the annual growth period for Atlantic salmon in the Saint-Jean
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

# ── yearly_growth_period ──────────────────────────────────────────────────────
# Purpose: Fits a cyclic GAM to daily water temperature for each year and
#   returns the first and last DOY at which the smoothed curve crosses
#   temp_threshold, defining the annual growth period.
# Parameters:
#   data            data frame with columns year, doy, temperature_c
#   temp_threshold  numeric; temperature threshold in degrees C (default: 6)
# Returns: data frame with columns year, start_doy, end_doy, start_date,
#   end_date; one row per year with a valid crossing.

yearly_growth_period <- function(data, temp_threshold = 6) {
  
  years_all    <- unique(data$year)
  results_list <- list()
  
  for (i in seq_along(years_all)) {
    
    data_i <- data %>%
      dplyr::filter(year == years_all[i])
    
    # Cyclic cubic spline ensures smooth continuity across the year boundary
    gam_i <- gam(
      temperature_c ~ s(doy, bs = "cc", k = 30),
      data   = data_i,
      method = "REML",
      knots  = list(doy = c(0.5, 366.5))
    )
    
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
      periods_above6 <- data.frame()
    } else {
      periods_above6 <- data.frame(
        year      = years_all[i],
        start_doy = first_day,
        end_doy   = last_day
      ) %>%
        mutate(
          start_date = as.Date(start_doy - 1,
                               origin = paste0(year, "-01-01")),
          end_date   = as.Date(end_doy   - 1,
                               origin = paste0(year, "-01-01"))
        )
    }
    
    results_list[[i]] <- periods_above6
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


# ── Detect annual growth periods ──────────────────────────────────────────────

range_stj <- yearly_growth_period(temp_stj, temp_threshold = temp_threshold)
range_tri <- yearly_growth_period(temp_tri, temp_threshold = temp_threshold)


# ── Export ────────────────────────────────────────────────────────────────────

write_fst(x    = range_stj,
          path = file.path(base_path,
                           "R/01_derived_data",
                           "Growth_Periods_StJean.fst"))

write_fst(x    = range_tri,
          path = file.path(base_path,
                           "R/01_derived_data",
                           "Growth_Periods_Trinite.fst"))






