


# ══════════════════════════════════════════════════════════════════════════════
# Computes the thermal gradient (delta T) between sea surface temperature (SST)
# during the smolt emigration window and mean river water temperature for the
# Saint-Jean and Trinite rivers, then saves one file per river x sea-age class.
#
# Inputs:
#   - R/00_raw_data/SST_ORAS5_Labrador/         : monthly NetCDF (.nc) files
#                                                  of ESA CCI L4 SST
#   - R/00_raw_data/Temperature_Discharge_River/ : CEQEAU daily temperature
#                                                  and discharge CSVs for
#                                                  Saint-Jean and Trinite
#   - R/00_raw_data/ReferenceMigrationLayers/   : migration corridor shapefiles
#                                                  (Parts) for both rivers
#   - R/01_derived_data/Year_Median_Date_StJean.fst   : smolt phenology (fst)
#   - R/01_derived_data/Year_Median_Date_Trinite.fst  : smolt phenology (fst)
#
# Outputs:
#   - R/03_clean_data/Sea/DeltaT_SST_River_stj_1SW.fst
#   - R/03_clean_data/Sea/DeltaT_SST_River_stj_2SW.fst
#   - R/03_clean_data/Sea/DeltaT_SST_River_tri_1SW.fst
#   - R/03_clean_data/Sea/DeltaT_SST_River_tri_2SW.fst
# ══════════════════════════════════════════════════════════════════════════════


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: LIBRARIES
# ══════════════════════════════════════════════════════════════════════════════

library(fst)
library(data.table)
library(dplyr)
library(lubridate)
library(slider)
library(tidyr)
library(sf)
library(terra)
library(ncdf4)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)




# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: WORKING DIRECTORY
# ══════════════════════════════════════════════════════════════════════════════

# ── Working directory ─────────────────────────────────────────────────────────
# Automatically detects the repository root based on the script location.
# No hardcoded path — works on any machine.
current_file <- rstudioapi::getActiveDocumentContext()$path
path_parts   <- unlist(strsplit(normalizePath(current_file, winslash = "/"), "/"))
target_index <- which(path_parts == "phd-chapter3-Atlantic-Salmon")
base_path    <- paste(path_parts[1:target_index], collapse = "/")
setwd(base_path)




# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3: DATA IMPORT
# ══════════════════════════════════════════════════════════════════════════════

# ── River temperature and discharge (Saint-Jean) ──────────────────────────────
# SOURCE: CEQEAU hydrological model output provided by MELCCFP.
# STATUS: Confidential — not included in this repository.
#         An example file with identical structure and simulated data is
#         available at the same path with the suffix _EXAMPLE.csv.
#         To request access, contact MELCCFP once the thesis is published.
temp_sj_raw <- read.csv(
  file.path(base_path,
            "R/00_raw_data/Temperature_Discharge_River",
            "CEQEAU_Temperature_Discharge_1979_2024_SJ.csv")
)

# ── River temperature and discharge (Trinite) ─────────────────────────────────
# SOURCE: CEQEAU hydrological model output provided by MELCCFP.
# STATUS: Confidential — not included in this repository.
#         An example file with identical structure and simulated data is
#         available at the same path with the suffix _EXAMPLE.csv.
#         To request access, contact MELCCFP once the thesis is published.
temp_tr_raw <- read.csv(
  file.path(base_path,
            "R/00_raw_data/Temperature_Discharge_River",
            "CEQEAU_Temperature_Discharge_1979_2024_TR.csv")
)

# ── Smolt emigration phenology (Saint-Jean) ───────────────────────────────────
# SOURCE: PIT-tag capture data processed in smolt phenology scripts.
# STATUS: Confidential — not included in this repository.
#         An example file with identical structure and simulated data is
#         available at the same path with the suffix _EXAMPLE.fst.
#         To request access, contact MELCCFP once the thesis is published.
median_smolt_date_stj <- read_fst(
  file.path(base_path,
            "R/01_derived_data",
            "Year_Median_Date_StJean.fst")
)

# ── Smolt emigration phenology (Trinite) ──────────────────────────────────────
# SOURCE: PIT-tag capture data processed in smolt phenology scripts.
# STATUS: Confidential — not included in this repository.
#         An example file with identical structure and simulated data is
#         available at the same path with the suffix _EXAMPLE.fst.
#         To request access, contact MELCCFP once the thesis is published.
median_smolt_date_tri <- read_fst(
  file.path(base_path,
            "R/01_derived_data",
            "Year_Median_Date_Trinite.fst")
)




# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4: CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

# Directory containing monthly ESA CCI SST NetCDF files
sst_dir <- file.path(base_path, "R/00_raw_data/SST_ORAS5_Labrador")

# Output directory for delta T results
output_dir <- file.path(base_path, "R/03_clean_data/Sea")

# SST directories per river (currently identical; separated for flexibility)
sst_dir_STJ <- sst_dir
sst_dir_TRI <- sst_dir

# Minimum number of valid pixels required to compute a spatial SST mean
MIN_VALID_PIXELS <- 3L

# FID value identifying the corridor polygon subset to use for each river
FID_STJ <- 1L
FID_TRI <- 1L




# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5: FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

# ── prepare_smolt_dates ───────────────────────────────────────────────────────
# Purpose: Selects relevant columns from the smolt phenology table and imputes
#          missing median/quantile dates using a centred 11-year rolling mean
#          (5 years before and 5 years after the missing value).
# Parameters:
#   smolt_data : data.frame — smolt phenology table with columns river, cohort,
#                median_date_doy, median_date, quantile_5_date, quantile_95_date
# Returns: data.frame with imputed median_date_doy and quantile date columns.

prepare_smolt_dates <- function(smolt_data) {
  smolt_data %>%
    dplyr::select(river, cohort, median_date_doy, median_date,
                  quantile_5_date, quantile_95_date) %>%
    arrange(river, cohort) %>%
    group_by(river) %>%
    mutate(
      median_date_doy_imputed = if_else(
        is.na(median_date_doy),
        round(slider::slide_dbl(
          median_date_doy,
          ~ mean(.x, na.rm = TRUE),
          .before = 5, .after = 5, .complete = TRUE
        )),
        median_date_doy
      ),
      median_date = if_else(
        is.na(median_date_doy),
        as.Date(median_date_doy_imputed - 1,
                origin = paste0(cohort, "-01-01")),
        median_date
      ),
      quantile_5_doy  = yday(quantile_5_date),
      quantile_95_doy = yday(quantile_95_date),
      quantile_5_doy_imputed = if_else(
        is.na(median_date_doy),
        round(slider::slide_dbl(
          quantile_5_doy,
          ~ mean(.x, na.rm = TRUE),
          .before = 5, .after = 5, .complete = TRUE
        )),
        quantile_5_doy
      ),
      quantile_95_doy_imputed = if_else(
        is.na(median_date_doy),
        round(slider::slide_dbl(
          quantile_95_doy,
          ~ mean(.x, na.rm = TRUE),
          .before = 5, .after = 5, .complete = TRUE
        )),
        quantile_95_doy
      ),
      quantile_5_date = if_else(
        is.na(median_date_doy),
        as.Date(quantile_5_doy_imputed - 1,
                origin = paste0(cohort, "-01-01")),
        quantile_5_date
      ),
      quantile_95_date = if_else(
        is.na(median_date_doy),
        as.Date(quantile_95_doy_imputed - 1,
                origin = paste0(cohort, "-01-01")),
        quantile_95_date
      ),
      median_date_doy = median_date_doy_imputed
    ) %>%
    dplyr::select(-median_date_doy_imputed, -quantile_5_doy,
                  -quantile_95_doy, -quantile_5_doy_imputed,
                  -quantile_95_doy_imputed) %>%
    ungroup()
}


# ── prepare_river_temp ────────────────────────────────────────────────────────
# Purpose: Standardises column names in a raw CEQEAU output table and parses
#          the date column, returning a tidy daily temperature series.
# Parameters:
#   raw     : data.frame — raw CEQEAU CSV loaded with read.csv()
#   obs_col : character  — column name for observed discharge
#   sim_col : character  — column name for simulated discharge
#   tw_col  : character  — column name for water temperature
# Returns: data.frame with columns date, temperature_c, year, month, day, doy.

prepare_river_temp <- function(raw, obs_col, sim_col, tw_col) {
  
  df <- raw
  colnames(df)[colnames(df) == obs_col] <- "observed_discharge_m3s"
  colnames(df)[colnames(df) == sim_col] <- "simulated_discharge_m3s"
  colnames(df)[colnames(df) == tw_col]  <- "temperature_c"
  colnames(df) <- tolower(colnames(df))
  
  # Force English locale to parse abbreviated month names (e.g., "Jan", "Feb")
  old_locale <- Sys.getlocale("LC_TIME")
  Sys.setlocale("LC_TIME", "C")
  df$date <- as.Date(df$date, format = "%d-%b-%Y")
  Sys.setlocale("LC_TIME", old_locale)
  
  df %>%
    dplyr::select(date, temperature_c) %>%
    mutate(
      year  = year(date),
      month = month(date),
      day   = day(date),
      doy   = yday(date)
    )
}


# ── mean_river_temp_smolting ──────────────────────────────────────────────────
# Purpose: Computes the yearly mean river water temperature across the smolt
#          emigration window (Q5 to Q95 of captures, inclusive).
# Parameters:
#   temperature_data : data.frame — daily river temperature with columns
#                      year, doy, temperature_c
#   smolt_dates      : data.frame — smolt phenology with columns cohort,
#                      quantile_5_date, quantile_95_date
# Returns: data.frame with columns year and mean_temperature_during_smolting.

mean_river_temp_smolting <- function(temperature_data, smolt_dates) {
  
  years <- unique(temperature_data$year)
  res   <- data.frame(year = years,
                      mean_temperature_during_smolting = NA_real_)
  
  for (i in seq_along(years)) {
    
    yr    <- years[i]
    win_i <- smolt_dates %>% filter(cohort == yr)
    
    if (nrow(win_i) == 0L) next
    
    temp_i <- temperature_data %>%
      filter(
        year == yr,
        doy  >= yday(win_i$quantile_5_date),
        doy  <= yday(win_i$quantile_95_date)
      )
    
    res$mean_temperature_during_smolting[res$year == yr] <-
      mean(temp_i$temperature_c, na.rm = TRUE)
  }
  
  res
}


# ── build_file_index_sst ──────────────────────────────────────────────────────
# Purpose: Scans a directory of monthly ESA CCI SST NetCDF archives and builds
#          a date-indexed table of individual daily files needed to cover the
#          smolt emigration windows.
# Parameters:
#   sst_dir     : character  — path to directory containing .nc archive files
#   smolt_dates : data.frame — smolt phenology with quantile_5_date and
#                              quantile_95_date columns
# Returns: data.frame with columns date, zip_path, nc_name, sorted by date.

build_file_index_sst <- function(sst_dir, smolt_dates) {
  
  # Determine the unique year-month combinations that need to be covered
  all_dates_needed <- smolt_dates %>%
    rowwise() %>%
    mutate(d = list(seq.Date(quantile_5_date, quantile_95_date,
                             by = "day"))) %>%
    tidyr::unnest(d) %>%
    mutate(yr = lubridate::year(d), mo = lubridate::month(d)) %>%
    dplyr::distinct(yr, mo) %>%
    dplyr::arrange(yr, mo)
  
  zip_files <- list.files(sst_dir, pattern = "\\.nc$", full.names = TRUE)
  
  if (length(zip_files) == 0L) {
    warning("No ZIP (.nc) files found in: ", sst_dir)
    return(data.frame(date     = as.Date(character()),
                      zip_path = character(),
                      nc_name  = character()))
  }
  
  rows <- vector("list", length(zip_files))
  
  for (k in seq_along(zip_files)) {
    
    f     <- zip_files[k]
    bname <- basename(f)
    
    yr_k <- as.integer(
      sub("ESACCI_SST_L4_(\\d{4})_(\\d{2})\\.nc", "\\1", bname)
    )
    mo_k <- as.integer(
      sub("ESACCI_SST_L4_(\\d{4})_(\\d{2})\\.nc", "\\2", bname)
    )
    
    if (is.na(yr_k) || is.na(mo_k)) next
    if (!any(all_dates_needed$yr == yr_k &
             all_dates_needed$mo == mo_k)) next
    
    contents <- tryCatch(unzip(f, list = TRUE), error = function(e) NULL)
    if (is.null(contents) || nrow(contents) == 0L) next
    
    nc_names  <- contents$Name
    date_strs <- substr(basename(nc_names), 1, 8)
    dates_i   <- as.Date(date_strs, format = "%Y%m%d")
    
    valid <- !is.na(dates_i)
    if (!any(valid)) next
    
    rows[[k]] <- data.frame(
      date     = dates_i[valid],
      zip_path = f,
      nc_name  = nc_names[valid],
      stringsAsFactors = FALSE
    )
  }
  
  idx <- do.call(rbind, rows[!sapply(rows, is.null)])
  
  if (is.null(idx) || nrow(idx) == 0L) {
    warning("File index is empty — no matching ZIP files processed.")
    return(data.frame(date     = as.Date(character()),
                      zip_path = character(),
                      nc_name  = character()))
  }
  
  idx[order(idx$date), ]
}


# ── load_daily_sst ────────────────────────────────────────────────────────────
# Purpose: Extracts one daily ESA CCI SST raster from a NetCDF archive,
#          corrects the spatial extent from the lon/lat metadata, and converts
#          values from Kelvin to Celsius.
# Parameters:
#   zip_path : character — path to the monthly .nc archive file
#   nc_name  : character — internal name of the daily NetCDF file to extract
# Returns: SpatRaster in EPSG:4326 with SST values in degrees Celsius,
#          or NULL on failure.

load_daily_sst <- function(zip_path, nc_name) {
  
  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  
  extracted <- tryCatch(
    unzip(zip_path, files = nc_name, exdir = tmp_dir),
    error = function(e) {
      warning("unzip failed: ", basename(zip_path), " / ", nc_name,
              " — ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(extracted) || length(extracted) == 0L) return(NULL)
  
  tmp_nc <- extracted[1]
  if (!file.exists(tmp_nc)) return(NULL)
  
  r_all <- tryCatch(terra::rast(tmp_nc), error = function(e) NULL)
  if (is.null(r_all)) return(NULL)
  
  # Select the analysed_sst layer; fall back to layer 1 if not found
  sst_idx <- grep("analysed_sst$", names(r_all), value = FALSE)
  if (length(sst_idx) == 0L) sst_idx <- 1L
  r_sst <- r_all[[sst_idx[1]]]
  
  # Recompute extent from NetCDF coordinate variables to avoid half-pixel shift
  nc_tmp <- tryCatch(ncdf4::nc_open(tmp_nc), error = function(e) NULL)
  if (!is.null(nc_tmp)) {
    lons_nc <- tryCatch(ncdf4::ncvar_get(nc_tmp, "lon"), error = function(e) NULL)
    lats_nc <- tryCatch(ncdf4::ncvar_get(nc_tmp, "lat"), error = function(e) NULL)
    ncdf4::nc_close(nc_tmp)
    
    if (!is.null(lons_nc) && !is.null(lats_nc)) {
      res_x <- mean(diff(sort(unique(lons_nc))), na.rm = TRUE) / 2
      res_y <- mean(diff(sort(unique(lats_nc))), na.rm = TRUE) / 2
      terra::ext(r_sst) <- terra::ext(
        min(lons_nc) - res_x, max(lons_nc) + res_x,
        min(lats_nc) - res_y, max(lats_nc) + res_y
      )
    }
  }
  
  terra::crs(r_sst) <- "EPSG:4326"
  
  # Convert from Kelvin to Celsius
  r_sst - 273.15
}


# ── spatial_mean_sst_one_day ──────────────────────────────────────────────────
# Purpose: Computes the spatial mean SST (Celsius) within a polygon for a
#          single daily raster extracted from an archive.
# Parameters:
#   zip_path    : character   — path to the monthly .nc archive
#   nc_name     : character   — internal name of the daily NetCDF to extract
#   polygon_vect: SpatVector  — study area polygon in EPSG:4326
#   min_pixels  : integer     — minimum number of valid pixels required
# Returns: numeric scalar (mean SST in Celsius), or NA_real_ on failure.

spatial_mean_sst_one_day <- function(zip_path, nc_name, polygon_vect,
                                     min_pixels = MIN_VALID_PIXELS) {
  
  r_sst_c <- load_daily_sst(zip_path, nc_name)
  if (is.null(r_sst_c)) return(NA_real_)
  
  r_masked <- tryCatch(
    terra::mask(terra::crop(r_sst_c, polygon_vect), polygon_vect),
    error = function(e) NULL
  )
  if (is.null(r_masked)) return(NA_real_)
  
  vals <- terra::values(r_masked)[, 1]
  
  mean(vals, na.rm = TRUE)
}


# ── mean_sst_window ───────────────────────────────────────────────────────────
# Purpose: Computes the mean daily spatial SST within a polygon across a date
#          range by iterating over the pre-built file index.
# Parameters:
#   file_index   : data.frame  — output of build_file_index_sst()
#   start_date   : Date        — first day of the window (inclusive)
#   end_date     : Date        — last day of the window (inclusive)
#   polygon_vect : SpatVector  — study area polygon in EPSG:4326
#   min_pixels   : integer     — minimum valid pixels passed downstream
# Returns: list with mean_sst (numeric) and n_days_valid (integer).

mean_sst_window <- function(file_index, start_date, end_date,
                            polygon_vect,
                            min_pixels = MIN_VALID_PIXELS) {
  
  idx <- file_index[
    !is.na(file_index$date) &
      file_index$date >= start_date &
      file_index$date <= end_date, ]
  
  if (nrow(idx) == 0L) return(list(mean_sst = NA_real_, n_days_valid = 0L))
  
  daily_means <- vapply(seq_len(nrow(idx)), function(i) {
    spatial_mean_sst_one_day(
      zip_path     = idx$zip_path[i],
      nc_name      = idx$nc_name[i],
      polygon_vect = polygon_vect,
      min_pixels   = min_pixels
    )
  }, FUN.VALUE = numeric(1L))
  
  list(
    mean_sst     = mean(daily_means),
    n_days_valid = length(daily_means)
  )
}


# ── compute_mean_sst_smolting ─────────────────────────────────────────────────
# Purpose: Iterates over cohort years and computes the mean SST within the
#          smolt emigration window (Q5-Q95) for a given river polygon.
# Parameters:
#   smolt_dates   : data.frame  — smolt phenology with cohort, quantile_5_date,
#                                 quantile_95_date
#   file_index    : data.frame  — output of build_file_index_sst()
#   polygon_wgs84 : sf object   — migration corridor polygon in EPSG:4326
#   min_pixels    : integer     — minimum valid pixels passed downstream
# Returns: data.frame with year, n_days_window, n_days_valid,
#          mean_sst_during_smolting.

compute_mean_sst_smolting <- function(smolt_dates, file_index,
                                      polygon_wgs84,
                                      min_pixels = MIN_VALID_PIXELS) {
  
  poly_v <- terra::vect(polygon_wgs84)
  
  res <- smolt_dates %>%
    dplyr::select(cohort, quantile_5_date, quantile_95_date) %>%
    dplyr::rename(year = cohort) %>%
    dplyr::mutate(
      n_days_window            = as.integer(
        difftime(quantile_95_date, quantile_5_date, units = "days")
      ) + 1L,
      n_days_valid             = NA_integer_,
      mean_sst_during_smolting = NA_real_
    )
  
  for (i in seq_len(nrow(res))) {
    
    out <- mean_sst_window(
      file_index   = file_index,
      start_date   = res$quantile_5_date[i],
      end_date     = res$quantile_95_date[i],
      polygon_vect = poly_v,
      min_pixels   = min_pixels
    )
    
    res$mean_sst_during_smolting[i] <- out$mean_sst
    res$n_days_valid[i]             <- out$n_days_valid
  }
  
  res %>%
    dplyr::select(year, n_days_window, n_days_valid, mean_sst_during_smolting)
}


# ── build_delta_T ─────────────────────────────────────────────────────────────
# Purpose: Joins the SST and river temperature tables and computes the thermal
#          gradient delta_T = mean SST (sea) - mean river temperature.
#          Positive values indicate warmer sea than river at smolting.
# Parameters:
#   sst_data       : data.frame — output of compute_mean_sst_smolting()
#   river_temp_data: data.frame — output of mean_river_temp_smolting()
#   river_label    : character  — river identifier ("stj" or "tri")
# Returns: data.frame with river, year, n_days_window, n_days_valid,
#          mean_sst_during_smolting, mean_river_temp, delta_T_sst_minus_river.

build_delta_T <- function(sst_data, river_temp_data, river_label) {
  
  sst_data %>%
    dplyr::left_join(
      river_temp_data %>%
        dplyr::rename(mean_river_temp = mean_temperature_during_smolting),
      by = "year"
    ) %>%
    dplyr::mutate(
      river                   = river_label,
      delta_T_sst_minus_river = mean_sst_during_smolting - mean_river_temp
    ) %>%
    dplyr::select(
      river, year,
      n_days_window, n_days_valid,
      mean_sst_during_smolting,
      mean_river_temp,
      delta_T_sst_minus_river
    )
}




# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6: MAIN SCRIPT
# ══════════════════════════════════════════════════════════════════════════════

# ── Prepare smolt emigration date windows ─────────────────────────────────────
# Impute missing phenology values with a centred 11-year rolling mean.

median_smolt_date_stj <- prepare_smolt_dates(median_smolt_date_stj)
median_smolt_date_tri <- prepare_smolt_dates(median_smolt_date_tri)


# ── Prepare daily river temperature series ────────────────────────────────────
# Standardise column names and parse dates for each river.

temp_sj <- prepare_river_temp(
  raw     = temp_sj_raw,
  obs_col = "Observed.Discharge..m3.s.",
  sim_col = "Simulated.Discharge..m3.s.",
  tw_col  = "Tw..degree.C."
)

temp_tr <- prepare_river_temp(
  raw     = temp_tr_raw,
  obs_col = "Observed.Discharge..RatingCurve...m3.s.",
  sim_col = "Simulated.Discharge..m3.s.",
  tw_col  = "Simulated.Tw..degree.C."
)


# ── Compute mean river temperature during the smolting window ─────────────────
# Window defined by Q5 to Q95 of smolt captures for each cohort year.

temp_during_smolting_stj <- mean_river_temp_smolting(temp_sj,
                                                     median_smolt_date_stj)
temp_during_smolting_tri <- mean_river_temp_smolting(temp_tr,
                                                     median_smolt_date_tri)


# ── Load migration corridor polygons ─────────────────────────────────────────
# Only the corridor segment identified by FID == 1 is used for each river.

poly_tri_wgs84 <- read_sf(
  file.path(base_path,
            "R/00_raw_data/ReferenceMigrationLayers",
            "TriniteMigrationCorridor_Parts.shp")
) %>%
  filter(FID == FID_TRI) %>%
  st_transform(crs = 4326)

poly_stj_wgs84 <- read_sf(
  file.path(base_path,
            "R/00_raw_data/ReferenceMigrationLayers",
            "StJeanMigrationCorridor_Parts.shp")
) %>%
  filter(FID == FID_STJ) %>%
  st_transform(crs = 4326)


# ── Build SST file indexes ────────────────────────────────────────────────────
# Index only the archive files that overlap with the smolting windows.

file_index_stj <- build_file_index_sst(sst_dir_STJ, median_smolt_date_stj)
file_index_tri <- build_file_index_sst(sst_dir_TRI, median_smolt_date_tri)


# ── Compute mean SST during smolting windows ──────────────────────────────────
# Spatially averaged over the corridor polygon for each cohort year.

sst_during_smolting_stj <- compute_mean_sst_smolting(
  smolt_dates   = median_smolt_date_stj,
  file_index    = file_index_stj,
  polygon_wgs84 = poly_stj_wgs84
)

sst_during_smolting_tri <- compute_mean_sst_smolting(
  smolt_dates   = median_smolt_date_tri,
  file_index    = file_index_tri,
  polygon_wgs84 = poly_tri_wgs84
)


# ── Verification plot — polygon subsets on SST map ───────────────────────────
# Visual check that both corridor polygons fall within the expected SST extent.

r_verify <- load_daily_sst(
  zip_path = file_index_stj$zip_path[1],
  nc_name  = file_index_stj$nc_name[1]
)

bbox_stj  <- sf::st_bbox(poly_stj_wgs84)
bbox_tri  <- sf::st_bbox(poly_tri_wgs84)
xlim_plot <- c(min(bbox_stj["xmin"], bbox_tri["xmin"]) - 2,
               max(bbox_stj["xmax"], bbox_tri["xmax"]) + 2)
ylim_plot <- c(min(bbox_stj["ymin"], bbox_tri["ymin"]) - 2,
               max(bbox_stj["ymax"], bbox_tri["ymax"]) + 2)

r_crop <- terra::crop(
  r_verify,
  terra::ext(xlim_plot[1], xlim_plot[2], ylim_plot[1], ylim_plot[2])
)

r_df <- as.data.frame(r_crop, xy = TRUE)
colnames(r_df)[3] <- "sst"

world     <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
provinces <- rnaturalearth::ne_states(country = "Canada", returnclass = "sf")

ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = sst), alpha = 0.85) +
  scale_fill_distiller(palette = "RdYlBu", name = "SST (degC)",
                       na.value = "transparent") +
  geom_sf(data = world,     fill = "grey85", colour = "grey60",
          linewidth = 0.3) +
  geom_sf(data = provinces, fill = NA,       colour = "grey50",
          linewidth = 0.2) +
  geom_sf(data = poly_stj_wgs84, fill = "#E69F00", colour = "#E69F00",
          alpha = 0.55, linewidth = 0.9) +
  geom_sf(data = poly_tri_wgs84, fill = "#56B4E9", colour = "#56B4E9",
          alpha = 0.55, linewidth = 0.9) +
  coord_sf(xlim = xlim_plot, ylim = ylim_plot, expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")


# ── Compute delta T (sea SST minus river temperature) ────────────────────────
# delta_T > 0 : sea is warmer than river at smolting (thermal shock to salmon)
# delta_T < 0 : sea is cooler than river at smolting

delta_T_stj <- build_delta_T(sst_during_smolting_stj,
                             temp_during_smolting_stj, "stj") %>%
  as.data.frame()

delta_T_tri <- build_delta_T(sst_during_smolting_tri,
                             temp_during_smolting_tri, "tri") %>%
  as.data.frame()

delta_T_all <- dplyr::bind_rows(delta_T_stj, delta_T_tri)


# ── Split by sea-age class and save ──────────────────────────────────────────
# delta_T is measured at the smolting year and is identical for 1SW and 2SW
# fish from the same cohort; rows are duplicated to match the survival table
# structure used downstream.

delta_T_stj_1SW <- delta_T_stj %>% dplyr::mutate(sw_type = "1SW")
delta_T_stj_2SW <- delta_T_stj %>% dplyr::mutate(sw_type = "2SW")
delta_T_tri_1SW <- delta_T_tri %>% dplyr::mutate(sw_type = "1SW")
delta_T_tri_2SW <- delta_T_tri %>% dplyr::mutate(sw_type = "2SW")

write_fst(delta_T_stj_1SW,
          file.path(output_dir, "DeltaT_SST_River_stj_1SW.fst"),
          compress = 100)
write_fst(delta_T_stj_2SW,
          file.path(output_dir, "DeltaT_SST_River_stj_2SW.fst"),
          compress = 100)
write_fst(delta_T_tri_1SW,
          file.path(output_dir, "DeltaT_SST_River_tri_1SW.fst"),
          compress = 100)
write_fst(delta_T_tri_2SW,
          file.path(output_dir, "DeltaT_SST_River_tri_2SW.fst"),
          compress = 100)



