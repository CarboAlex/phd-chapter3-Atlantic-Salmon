

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: HEADER
# ══════════════════════════════════════════════════════════════════════════════
#
# Imports ERA5 air temperature raster data (GRIB format) from the Copernicus
# Climate Data Store, converts hourly values from Kelvin to Celsius, and
# exports daily mean temperatures for the Saint-Jean and Trinite rivers.
#
# Data source: Copernicus Climate Data Store (CDS)
#   https://cds.climate.copernicus.eu/requests?tab=all
#
# Daily mean temperature is defined as the average of the daily minimum and
# maximum temperatures, consistent with the convention used by Canadian
# climate stations.
#
# Inputs:
#   - R/00_raw_data/Air_Temperature_Raster/tri_1975_2008.zip
#   - R/00_raw_data/Air_Temperature_Raster/tri_2009_2025.zip
#       ERA5 hourly air temperature rasters for the Trinite area (GRIB in ZIP)
#   - R/00_raw_data/Air_Temperature_Raster/stj_1978_2024.zip
#       ERA5 hourly air temperature rasters for the Saint-Jean area (GRIB in
#       ZIP)
#
# Outputs:
#   - R/01_derived_data/Winter_Temperature/tri_air_temp_daily_mean_rasters.fst
#       Daily mean air temperature for the Trinite River (FST)
#   - R/01_derived_data/Winter_Temperature/stj_air_temp_daily_mean_rasters.fst
#       Daily mean air temperature for the Saint-Jean River (FST)
#
# ══════════════════════════════════════════════════════════════════════════════


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: LIBRARIES
# ══════════════════════════════════════════════════════════════════════════════

library(data.table)
library(terra)
library(ggplot2)
library(ggpubr)
library(rnaturalearth)
library(rnaturalearthdata)
library(fst)
library(sf)


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

# ── ERA5 raster files — Trinite ───────────────────────────────────────────────
# SOURCE: Copernicus Climate Data Store (ERA5 reanalysis)
# STATUS: Too large for GitHub — not included in this repository.
#         Download from: https://cds.climate.copernicus.eu/requests?tab=all

# GRIB files are extracted to a temporary directory and deleted after import
tmpdir <- tempdir()

unzip(file.path(base_path,
                "R/00_raw_data/Air_Temperature_Raster/tr_1975_2008.zip"),
      exdir = tmpdir)
unzip(file.path(base_path,
                "R/00_raw_data/Air_Temperature_Raster/tr_2009_2025.zip"),
      exdir = tmpdir)

gribfile_tri_1975_2008 <- list.files(tmpdir, pattern = "tr_1975_2008",
                                     full.names = TRUE)
gribfile_tri_2009_2025 <- list.files(tmpdir, pattern = "tr_2009_2025",
                                     full.names = TRUE)

tri_1975_2008 <- rast(gribfile_tri_1975_2008)
tri_2009_2025 <- rast(gribfile_tri_2009_2025)


# ── ERA5 raster files — Saint-Jean ────────────────────────────────────────────
# SOURCE: Copernicus Climate Data Store (ERA5 reanalysis)
# STATUS: Too large for GitHub — not included in this repository.
#         Download from: https://cds.climate.copernicus.eu/requests?tab=all

unzip(file.path(base_path,
                "R/00_raw_data/Air_Temperature_Raster/sj_1978_2024.zip"),
      exdir = tmpdir)

gribfile_stj_1978_2024 <- list.files(tmpdir, pattern = "sj_1978_2024",
                                     full.names = TRUE)

stj_1978_2024 <- rast(gribfile_stj_1978_2024)

# Schedule GRIB deletion after the session ends to free disk space
on.exit({
  if (file.exists(gribfile)) {
    file.remove(gribfile)
  }
}, add = TRUE)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5: CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

# Output directory for derived temperature files
out_path <- file.path(base_path, "R/01_derived_data/Winter_Temperature")

# River mouth coordinates used to verify raster spatial coverage
tri_mouth_lon <- -67.305774
tri_mouth_lat  <-  49.418871
stj_mouth_lon <- -64.382411
stj_mouth_lat  <-  48.780331


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6: FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

# ── extract_air_temp_raster ───────────────────────────────────────────────────
# Purpose: Converts a multi-layer terra SpatRaster of hourly air temperatures
#   into a long-format data.table with date, hour, coordinates, spatial
#   resolution, and temperature in Celsius.
# Parameters:
#   data   SpatRaster with layers named by datetime (time() must be set)
# Returns: data.table with columns lon, lat, date, hour, resolution_x,
#   resolution_y, temp (degrees C).

extract_air_temp_raster <- function(data) {
  
  data_dt <- as.data.table(as.data.frame(data, xy = TRUE, na.rm = FALSE))
  
  # Rename columns: first two are coordinates, remainder are timestamped layers
  setnames(data_dt,
           c("lon", "lat",
             paste("temp", format(time(data), "%Y%m%d%H%M"), sep = "_")))
  
  # Pivot from wide to long so each row is one cell-timestep observation
  data_dt_long <- melt(data_dt, id.vars = c("lon", "lat"),
                       variable.name = "layer", value.name = "temp")
  data_dt_long[, layer := as.character(layer)]
  
  # Attach spatial resolution from the original raster
  data_dt_long$resolution_x <- res(data)[1]
  data_dt_long$resolution_y <- res(data)[2]
  
  # Build a lookup table mapping each layer name to its parsed datetime
  layers     <- unique(data_dt_long$layer)
  time_strs  <- sub(".*_", "", layers)
  datetimes  <- as.POSIXct(time_strs, format = "%Y%m%d%H%M", tz = "UTC")
  lookup     <- data.table(layer = layers, datetime = datetimes)
  lookup[, date := as.Date(datetime)]
  lookup[, hour := format(datetime, "%H:%M")]
  
  data_dt_long <- merge(data_dt_long, lookup, by = "layer", all.x = TRUE)
  
  data_dt_long <- data_dt_long[,
                               c("lon", "lat", "date", "hour", "resolution_x", "resolution_y", "temp")]
  
  # Convert from Kelvin to Celsius
  data_dt_long[, temp := temp - 273.15]
  
  return(data_dt_long)
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7: MAIN SCRIPT
# ══════════════════════════════════════════════════════════════════════════════

# ── Verify raster spatial coverage ───────────────────────────────────────────

canada <- ne_countries(scale = "large", country = "canada", returnclass = "sf")

# Use only the first layer for plotting to keep memory usage low
tri_rast <- tri_1975_2008[[1]]
stj_rast <- stj_1978_2024[[1]]

tri_df        <- as.data.frame(tri_rast, xy = TRUE, na.rm = FALSE)
names(tri_df)[3] <- "value"
stj_df        <- as.data.frame(stj_rast, xy = TRUE, na.rm = FALSE)
names(stj_df)[3] <- "value"

tri_mouth <- st_as_sf(data.frame(lon = tri_mouth_lon, lat = tri_mouth_lat),
                      coords = c("lon", "lat"), crs = 4326)
stj_mouth <- st_as_sf(data.frame(lon = stj_mouth_lon, lat = stj_mouth_lat),
                      coords = c("lon", "lat"), crs = 4326)

ggarrange(
  ggplot() +
    geom_sf(data = canada, fill = "grey90", color = "black") +
    geom_raster(data = tri_df, aes(x = x, y = y, fill = value),
                alpha = 0.4) +
    geom_sf(data = tri_mouth, color = "red", size = 2) +
    scale_fill_viridis_c(option = "plasma", na.value = NA, name = "Temp") +
    coord_sf(xlim = c(-69, -66.5), ylim = c(48.5, 50.5), expand = FALSE) +
    theme_minimal(),
  
  ggplot() +
    geom_sf(data = canada, fill = "grey90", color = "black") +
    geom_raster(data = stj_df, aes(x = x, y = y, fill = value),
                alpha = 0.4) +
    geom_sf(data = stj_mouth, color = "red", size = 2) +
    scale_fill_viridis_c(option = "plasma", na.value = NA, name = "Temp") +
    coord_sf(xlim = c(-66, -63.7), ylim = c(48, 49.5), expand = FALSE) +
    theme_minimal(),
  
  ncol = 2)


# ── Extract hourly temperature time series from rasters ───────────────────────

tri_1975_2008_dt <- extract_air_temp_raster(data = tri_1975_2008)
tri_2009_2025_dt <- extract_air_temp_raster(data = tri_2009_2025)
tri_dt           <- rbind(tri_1975_2008_dt, tri_2009_2025_dt)

stj_dt <- extract_air_temp_raster(data = stj_1978_2024)


# ── Compute daily mean temperature ────────────────────────────────────────────

# Daily mean is computed as the average of the daily min and max,
# consistent with the convention used by Canadian climate stations
tri_daily_mean <- tri_dt[, .(
  Raster_Mean_T = mean(
    (max(temp, na.rm = TRUE) + min(temp, na.rm = TRUE)) / 2)
), by = date]

stj_daily_mean <- stj_dt[, .(
  Raster_Mean_T = mean(
    (max(temp, na.rm = TRUE) + min(temp, na.rm = TRUE)) / 2)
), by = date]


# ── Export ────────────────────────────────────────────────────────────────────

write_fst(tri_daily_mean,
          file.path(out_path, "tr_air_temp_daily_mean_rasters.fst"))

write_fst(stj_daily_mean,
          file.path(out_path, "sj_air_temp_daily_mean_rasters.fst"))

