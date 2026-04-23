

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: HEADER
# ══════════════════════════════════════════════════════════════════════════════
#
# Imports, merges, and cleans daily air temperature records from multiple
# climate stations near the Saint-Jean and Trinite rivers, tests inter-station
# correlations, and exports daily averages across stations for each river.
#
# Note: mean temperature (degrees C) is defined as the average of the daily
# maximum and minimum temperatures at a given location.
#
# Data source: ClimateData.ca (https://climatedata.ca/)
#
# Inputs:
#   - R/00_raw_data/Air_Temperature_Stations/*Trinite*.csv
#       Daily air temperature records for stations near the Trinite River (CSV)
#   - R/00_raw_data/Air_Temperature_Stations/*StJean*.csv
#       Daily air temperature records for stations near the Saint-Jean River
#       (CSV)
#
# Outputs:
#   - R/01_derived_data/Winter_Temperature/tri_air_temp_daily_mean_stations.fst
#       Daily mean of mean/min/max air temperature across Trinite stations (FST)
#   - R/01_derived_data/Winter_Temperature/stj_air_temp_daily_mean_stations.fst
#       Daily mean of mean/min/max air temperature across Saint-Jean stations
#       (FST)
#
# ══════════════════════════════════════════════════════════════════════════════


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: LIBRARIES
# ══════════════════════════════════════════════════════════════════════════════

library(data.table)
library(ggplot2)
library(ggpubr)
library(lubridate)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(fst)


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

# ── Air temperature station files — Trinite ───────────────────────────────────
# SOURCE: ClimateData.ca — Canadian climate station records
# STATUS: Included in this repository

tri_files_path <- list.files(
  file.path(base_path, "R/00_raw_data/Air_Temperature_Stations"),
  pattern    = "Trinite",
  full.names = TRUE)

tri_temp_list <- lapply(tri_files_path, function(x) { fread(x) })
tri_temp      <- rbindlist(tri_temp_list)


# ── Air temperature station files — Saint-Jean ────────────────────────────────
# SOURCE: ClimateData.ca — Canadian climate station records
# STATUS: Included in this repository

stj_files_path <- list.files(
  file.path(base_path, "R/00_raw_data/Air_Temperature_Stations"),
  pattern    = "StJean",
  full.names = TRUE)

stj_temp_list <- lapply(stj_files_path, function(x) { fread(x) })
stj_temp      <- rbindlist(stj_temp_list)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5: CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

# Columns retained from the raw station files
cols_keep <- c("x", "y", "STATION_NAME", "LOCAL_DATE",
               "LOCAL_YEAR", "LOCAL_MONTH", "LOCAL_DAY",
               "MEAN_TEMPERATURE", "MIN_TEMPERATURE", "MAX_TEMPERATURE")

# Output directory for derived temperature files
out_path <- file.path(base_path, "R/01_derived_data/Winter_Temperature")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6: FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

# ── cor_test_plot ─────────────────────────────────────────────────────────────
# Purpose: Runs a Pearson correlation test between two numeric vectors and
#   returns a scatter plot with a regression line and r, p, and n annotations.
# Parameters:
#   var_1   numeric vector (x-axis variable)
#   var_2   numeric vector (y-axis variable)
# Returns: a ggplot object.

cor_test_plot <- function(var_1, var_2) {
  
  cor_res   <- cor.test(var_1, var_2)
  r_val     <- round(cor_res$estimate, 3)
  p_val     <- signif(cor_res$p.value, 3)
  n_val     <- cor_res$parameter + 2  # degrees of freedom + 2 = n
  label_txt <- paste0(" r = ", r_val, ",\n p = ", p_val, ",\n n = ", n_val)
  
  ggplot(mapping = aes(x = var_1, y = var_2)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "red") +
    annotate("text",
             x     = min(var_1, na.rm = TRUE),
             y     = max(var_2, na.rm = TRUE),
             label = label_txt,
             hjust = 0, vjust = 1, size = 5) +
    theme_bw()
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7: MAIN SCRIPT
# ══════════════════════════════════════════════════════════════════════════════

# ── Clean and reformat Trinite station data ───────────────────────────────────

tri_temp <- tri_temp[, ..cols_keep]

tri_temp[, `:=`(
  x                = as.numeric(x),
  y                = as.numeric(y),
  STATION_NAME     = as.factor(STATION_NAME),
  LOCAL_DATE       = as.IDate(LOCAL_DATE, format = "%Y-%m-%d"),
  MEAN_TEMPERATURE = as.numeric(MEAN_TEMPERATURE),
  MIN_TEMPERATURE  = as.numeric(MIN_TEMPERATURE),
  MAX_TEMPERATURE  = as.numeric(MAX_TEMPERATURE)
)]

# Standardize station names for safe use as column names
tri_temp$STATION_NAME <- gsub(" ", "_", tri_temp$STATION_NAME)
tri_temp$STATION_NAME <- gsub("-", "_", tri_temp$STATION_NAME)


# ── Clean and reformat Saint-Jean station data ────────────────────────────────

stj_temp <- stj_temp[, ..cols_keep]

stj_temp[, `:=`(
  x                = as.numeric(x),
  y                = as.numeric(y),
  STATION_NAME     = as.factor(STATION_NAME),
  LOCAL_DATE       = as.IDate(LOCAL_DATE, format = "%Y-%m-%d"),
  MEAN_TEMPERATURE = as.numeric(MEAN_TEMPERATURE),
  MIN_TEMPERATURE  = as.numeric(MIN_TEMPERATURE),
  MAX_TEMPERATURE  = as.numeric(MAX_TEMPERATURE)
)]

stj_temp$STATION_NAME <- gsub(" ", "_", stj_temp$STATION_NAME)


# ── Station location maps ─────────────────────────────────────────────────────

canada <- ne_countries(country = "canada", returnclass = "sf", scale = 10)

location_station_tri <- ggplot() +
  geom_sf(data = canada, fill = "grey95", color = "grey50") +
  geom_sf(data = st_as_sf(unique(tri_temp[, .(x, y, STATION_NAME)]),
                          coords = c("x", "y"), crs = 4326),
          size = 3) +
  coord_sf(xlim = c(-68.5, -66.5), ylim = c(48.5, 50.0)) +
  theme_minimal()

location_station_stj <- ggplot() +
  geom_sf(data = canada, fill = "grey95", color = "grey50") +
  geom_sf(data = st_as_sf(unique(stj_temp[, .(x, y, STATION_NAME)]),
                          coords = c("x", "y"), crs = 4326),
          size = 3) +
  coord_sf(xlim = c(-65.5, -64), ylim = c(48.25, 49.21)) +
  theme_minimal()

ggarrange(location_station_tri, location_station_stj, ncol = 2)


# ── Data availability by station and year ────────────────────────────────────

avail_stj <- stj_temp[, .(AVAILABLE = 1), by = .(STATION_NAME, LOCAL_DATE)]
avail_tri <- tri_temp[, .(AVAILABLE = 1), by = .(STATION_NAME, LOCAL_DATE)]

breaks_2yr <- seq(
  from = as.Date("1975-01-01"),
  to   = as.Date("2025-01-01"),
  by   = "2 years"
)

ggarrange(
  ggplot(avail_stj, aes(x = LOCAL_DATE, y = STATION_NAME, fill = AVAILABLE)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "black", guide = "none") +
    scale_x_date(breaks = breaks_2yr, labels = format(breaks_2yr, "%Y")) +
    labs(x = "Year", y = "Station", title = "St-Jean") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0),
          axis.text.y = element_text(size = 8)),
  
  ggplot(avail_tri, aes(x = LOCAL_DATE, y = STATION_NAME, fill = AVAILABLE)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "black", guide = "none") +
    scale_x_date(breaks = breaks_2yr, labels = format(breaks_2yr, "%Y")) +
    labs(x = "Year", y = "Station", title = "Trinite") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0),
          axis.text.y = element_text(size = 8)),
  
  ncol = 1
)



# ── Exploring missing data ────────────────────────────────────

summarize_coverage_v2 <- function(avail_dt, site_name) {
  date_range <- seq(min(avail_dt$LOCAL_DATE), max(avail_dt$LOCAL_DATE), by = "day")
  total_days <- length(date_range)
  
  # Couverture par station individuelle
  per_station <- avail_dt[, .(
    covered_days = sum(AVAILABLE == 1),
    first_year = year(min(LOCAL_DATE)),
    last_year  = year(max(LOCAL_DATE))
  ), by = STATION_NAME]
  
  per_station[, coverage_pct := round(covered_days / total_days * 100, 1)]
  per_station[, n_years      := last_year - first_year + 1]
  
  data.table(
    Site                  = site_name,
    Periode               = paste(min(per_station$first_year), max(per_station$last_year), sep = "–"),
    N_stations            = nrow(per_station),
    Couverture_moy_pct    = round(mean(per_station$coverage_pct), 1),
    Couverture_sd_pct     = round(sd(per_station$coverage_pct), 1),
    Duree_moy_ans         = round(mean(per_station$n_years), 1),
    Duree_sd_ans          = round(sd(per_station$n_years), 1)
  )
}

rbind(
  summarize_coverage_v2(avail_stj %>% 
                          dplyr::filter(LOCAL_DATE >= "1980-01-01"), 
                        "St-Jean"),
  summarize_coverage_v2(avail_tri %>% 
                          dplyr::filter(LOCAL_DATE >= "1980-01-01"),
                        "Trinite")
)


check_combined_coverage <- function(avail_dt, site_name, start_year = 1980) {
  
  avail_dt <- avail_dt[year(LOCAL_DATE) >= start_year]
  date_range <- seq(as.IDate(paste0(start_year, "-01-01")), 
                    max(avail_dt$LOCAL_DATE), by = "day")
  total_days <- length(date_range)
  
  covered <- avail_dt[AVAILABLE == 1, .(n_stations = .N), by = LOCAL_DATE]
  
  covered_days    <- nrow(covered)
  missing_days    <- total_days - covered_days
  single_days     <- covered[n_stations == 1, .N]
  
  cat(site_name, ":\n")
  cat("  Couverture totale       :", round(covered_days / total_days * 100, 1), "%\n")
  cat("  Jours non disponibles   :", missing_days, "\n")
  cat("  Jours 1 seule station   :", round(single_days / covered_days * 100, 1), "% des jours couverts\n")
}

check_combined_coverage(avail_stj, "St-Jean")
check_combined_coverage(avail_tri, "Trinite")


# ── Inter-station temperature correlations ────────────────────────────────────

# Reshape to wide format for pairwise correlation plotting
tri_mean_wide <- dcast(tri_temp, LOCAL_DATE ~ STATION_NAME,
                       value.var = "MEAN_TEMPERATURE", fun.aggregate = mean)
tri_min_wide  <- dcast(tri_temp, LOCAL_DATE ~ STATION_NAME,
                       value.var = "MIN_TEMPERATURE",  fun.aggregate = mean)
tri_max_wide  <- dcast(tri_temp, LOCAL_DATE ~ STATION_NAME,
                       value.var = "MAX_TEMPERATURE",  fun.aggregate = mean)

stj_mean_wide <- dcast(stj_temp, LOCAL_DATE ~ STATION_NAME,
                       value.var = "MEAN_TEMPERATURE", fun.aggregate = mean)
stj_min_wide  <- dcast(stj_temp, LOCAL_DATE ~ STATION_NAME,
                       value.var = "MIN_TEMPERATURE",  fun.aggregate = mean)
stj_max_wide  <- dcast(stj_temp, LOCAL_DATE ~ STATION_NAME,
                       value.var = "MAX_TEMPERATURE",  fun.aggregate = mean)

# Mean temperature correlations — Trinite
ggarrange(
  cor_test_plot(tri_mean_wide$GODBOUT,      tri_mean_wide$PENTECOTE),
  cor_test_plot(tri_mean_wide$GODBOUT,      tri_mean_wide$POINTE_DES_MONTS),
  cor_test_plot(tri_mean_wide$PENTECOTE,    tri_mean_wide$POINTE_DES_MONTS),
  ncol = 2, nrow = 2)

# Minimum temperature correlations — Trinite
ggarrange(
  cor_test_plot(tri_min_wide$GODBOUT,       tri_min_wide$PENTECOTE),
  cor_test_plot(tri_min_wide$GODBOUT,       tri_min_wide$POINTE_DES_MONTS),
  cor_test_plot(tri_min_wide$PENTECOTE,     tri_min_wide$POINTE_DES_MONTS),
  ncol = 2, nrow = 2)

# Maximum temperature correlations — Trinite
ggarrange(
  cor_test_plot(tri_max_wide$GODBOUT,       tri_max_wide$PENTECOTE),
  cor_test_plot(tri_max_wide$GODBOUT,       tri_max_wide$POINTE_DES_MONTS),
  cor_test_plot(tri_max_wide$PENTECOTE,     tri_max_wide$POINTE_DES_MONTS),
  ncol = 2, nrow = 2)

# Mean temperature correlations — Saint-Jean
ggarrange(
  cor_test_plot(stj_mean_wide$FAREWELL_COVE, stj_mean_wide$FONTENELLE),
  cor_test_plot(stj_mean_wide$FAREWELL_COVE, stj_mean_wide$GASPE),
  cor_test_plot(stj_mean_wide$FAREWELL_COVE, stj_mean_wide$GASPE_A),
  cor_test_plot(stj_mean_wide$FONTENELLE,    stj_mean_wide$GASPE),
  cor_test_plot(stj_mean_wide$FONTENELLE,    stj_mean_wide$GASPE_A),
  cor_test_plot(stj_mean_wide$GASPE,         stj_mean_wide$GASPE_A),
  cor_test_plot(stj_mean_wide$GASPE_A,       stj_mean_wide$GASPE_AIRPORT),
  ncol = 3, nrow = 3)

# Minimum temperature correlations — Saint-Jean
ggarrange(
  cor_test_plot(stj_min_wide$FAREWELL_COVE, stj_min_wide$FONTENELLE),
  cor_test_plot(stj_min_wide$FAREWELL_COVE, stj_min_wide$GASPE),
  cor_test_plot(stj_min_wide$FAREWELL_COVE, stj_min_wide$GASPE_A),
  cor_test_plot(stj_min_wide$FONTENELLE,    stj_min_wide$GASPE),
  cor_test_plot(stj_min_wide$FONTENELLE,    stj_min_wide$GASPE_A),
  cor_test_plot(stj_min_wide$GASPE,         stj_min_wide$GASPE_A),
  cor_test_plot(stj_min_wide$GASPE_A,       stj_min_wide$GASPE_AIRPORT),
  ncol = 3, nrow = 3)

# Maximum temperature correlations — Saint-Jean
ggarrange(
  cor_test_plot(stj_max_wide$FAREWELL_COVE, stj_max_wide$FONTENELLE),
  cor_test_plot(stj_max_wide$FAREWELL_COVE, stj_max_wide$GASPE),
  cor_test_plot(stj_max_wide$FAREWELL_COVE, stj_max_wide$GASPE_A),
  cor_test_plot(stj_max_wide$FONTENELLE,    stj_max_wide$GASPE),
  cor_test_plot(stj_max_wide$FONTENELLE,    stj_max_wide$GASPE_A),
  cor_test_plot(stj_max_wide$GASPE,         stj_max_wide$GASPE_A),
  cor_test_plot(stj_max_wide$GASPE_A,       stj_max_wide$GASPE_AIRPORT),
  ncol = 3, nrow = 3)


# ── Average daily temperature across stations ─────────────────────────────────

# For each date, average all available stations to obtain a single daily value
tri_temp_summary <- tri_temp[, .(
  Station_Mean_T = mean(MEAN_TEMPERATURE, na.rm = TRUE),
  Station_Min_T  = mean(MIN_TEMPERATURE,  na.rm = TRUE),
  Station_Max_T  = mean(MAX_TEMPERATURE,  na.rm = TRUE)
), by = LOCAL_DATE]
tri_temp_summary <- tri_temp_summary[order(tri_temp_summary$LOCAL_DATE), ]

stj_temp_summary <- stj_temp[, .(
  Station_Mean_T = mean(MEAN_TEMPERATURE, na.rm = TRUE),
  Station_Min_T  = mean(MIN_TEMPERATURE,  na.rm = TRUE),
  Station_Max_T  = mean(MAX_TEMPERATURE,  na.rm = TRUE)
), by = LOCAL_DATE]
stj_temp_summary <- stj_temp_summary[order(stj_temp_summary$LOCAL_DATE), ]


# ── Export ────────────────────────────────────────────────────────────────────

write_fst(tri_temp_summary,
          file.path(out_path, "tr_air_temp_daily_mean_stations.fst"))

write_fst(stj_temp_summary,
          file.path(out_path, "sj_air_temp_daily_mean_stations.fst"))

