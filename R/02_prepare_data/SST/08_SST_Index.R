


# ── Working directory ──────────────────────────────────────────────────────
# Automatically detects the repository root based on the script location.
# No hardcoded path — works on any machine.
current_file <- rstudioapi::getActiveDocumentContext()$path
path_parts   <- unlist(strsplit(normalizePath(current_file, winslash = "/"), "/"))
target_index <- which(path_parts == "phd-chapter3-Atlantic-Salmon")
base_path    <- paste(path_parts[1:target_index], collapse = "/")
setwd(base_path)


# |=======================|
# |    SETTINGS           |
# |=======================|

download_dir <- file.path(base_path, "R/00_raw_data/SST_ORAS5_Labrador")
output_base  <- file.path(base_path, "R/01_derived_data/SST_Index")
shp_path_dir <- file.path(base_path, "R/00_raw_data/ReferenceMigrationLayers")

# Temporary folder for unzipped files
tmp_dir      <- file.path(download_dir, "tmp_processing")

# Parallel workers (used for daily files within each month)
n_workers <- 4

# Minimum valid cache size (bytes)
min_cache_size <- 1000


# |========================|
# |     LOAD LIBRARIES     |
# |========================|

library(data.table)
library(fst)
library(terra)
library(sf)
library(lubridate)
library(future.apply)
library(ncdf4)

if (!dir.exists(output_base)) dir.create(output_base, recursive = TRUE)
if (!dir.exists(tmp_dir))     dir.create(tmp_dir, recursive = TRUE)


# |========================|
# |     LOAD FUNCTIONS     |
# |========================|

thermal_habitat <- function(sst, ice, year_sw, lon, lat, res_x, res_y) {
  
  # Description:
  #   Calculates daily thermal habitat index based on SST, ice cover and sea-winter year.
  #   Also returns the area (km²) of each grid cell for weighting.
  #
  # Arguments:
  #   - sst: numeric vector. Sea surface temperature values.
  #   - ice: numeric vector. Sea ice concentration (0–1).
  #   - year_sw: integer. Sea-winter category (1 or 2).
  #   - lon, lat: numeric vectors. Coordinates of points (° WGS84).
  #   - res_x, res_y: numeric. Resolution of cells (°).
  #
  # Returns:
  #   - A list with:
  #       * growth_weighted: growth index weighted by cell area.
  #       * cell_area_km2: area of each grid cell in km².
  
  # Approximate cell size in meters & compute area
  dx_m <- res_x * cos(lat * pi / 180) * 111320
  dy_m <- res_y * 111320
  cell_area_km2 <- (dx_m * dy_m) / 1e6
  
  # Set thermal parameters based on sea-winter year
  if (year_sw == 1L) {
    t_min <- 1; t_opt1 <- 6; t_opt2 <- 14; t_max <- 18
  } else if (year_sw == 2L) {
    t_min <- 1; t_opt1 <- 6; t_opt2 <- 12; t_max <- 15
  } else stop("year_sw must be 1 or 2")
  
  # Precompute slopes and intercepts
  a_low  <- 1 / (t_opt1 - t_min)
  b_low  <- 1 - a_low * t_opt1
  a_high <- 1 / (t_opt2 - t_max)
  b_high <- 1 - a_high * t_opt2
  
  # Compute growth index per cell
  growth <- numeric(length(sst))
  growth[] <- NA_real_
  
  valid <- ice < 0.5 & !is.na(sst)
  
  # Growth areas
  growth[valid & sst >= t_opt1 & sst <= t_opt2] <- 1
  growth[valid & sst >= t_min & sst < t_opt1]   <- a_low * sst[valid & sst >= t_min & sst < t_opt1] + b_low
  growth[valid & sst > t_opt2 & sst <= t_max]   <- a_high * sst[valid & sst > t_opt2 & sst <= t_max] + b_high
  growth[!valid | sst < 1 | sst > t_max]        <- 0
  
  
  # Weight by cell area
  growth_weighted <- growth * cell_area_km2
  
  return(list(growth_weighted = growth_weighted, cell_area_km2 = cell_area_km2))
}


# |============================================================|
# |     DEFINE POLYGONS                                       |
# |============================================================|

shape_files_standard <- c(
  LabradorSea      = "LabradorSea.shp",
  LabradorEntrance = "LabradorSeaEntrance.shp",
  BelleIsle        = "BelleIsleStrait.shp",
  TriniteCorridor  = "TriniteMigrationCorridor_Parts.shp",
  StJeanCorridor   = "StJeanMigrationCorridor_Parts.shp"
)

load_standard_shp <- function(fname) {
  shp <- st_read(file.path(shp_path_dir, fname), quiet = TRUE) |>
    st_transform(4326) |>
    st_make_valid()
  if ("name" %in% names(shp)) names(shp)[names(shp) == "name"] <- "FID"
  if (!"FID" %in% names(shp)) shp$FID <- 0L
  shp$FID <- as.integer(shp$FID)
  shp
}

all_polygons  <- lapply(shape_files_standard,  load_standard_shp)

# Pre-compute per-polygon, per-FID: output dir + WKT geometry + FID type
poly_meta <- lapply(names(all_polygons), function(pname) {
  shp      <- all_polygons[[pname]]
  out_dir  <- file.path(output_base, paste0("SST_Index_", pname))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  fid_list <- lapply(seq_len(nrow(shp)), function(i) {
    list(
      fid = shp$FID[i],
      wkt = st_as_text(st_geometry(shp[i, ]))
    )
  })
  
  list(
    name       = pname,
    out_dir    = out_dir,
    fid_list   = fid_list,
    fid_is_char = is.character(shp$FID)
  )
})
names(poly_meta) <- names(all_polygons)



# |============================================================|
# |     LIST AVAILABLE ZIP FILES                              |
# |============================================================|

zip_files <- list.files(download_dir,
                        pattern = "ESACCI_SST_L4_\\d{4}_\\d{2}\\.nc$",
                        full.names = TRUE)
if (length(zip_files) == 0) stop("No ESACCI ZIP files found.")

dt_zips <- data.table(path = zip_files, fname = basename(zip_files))
dt_zips[, year  := as.integer(sub("ESACCI_SST_L4_(\\d{4})_(\\d{2})\\.nc",
                                  "\\1", fname))]
dt_zips[, month := as.integer(sub("ESACCI_SST_L4_(\\d{4})_(\\d{2})\\.nc",
                                  "\\2", fname))]
setorder(dt_zips, year, month)

message("ZIP files: ", nrow(dt_zips),
        " (", min(dt_zips$year), "-01 to ",
        max(dt_zips$year), "-12)")


# |============================================================|
# |     HELPER: process one daily .nc for all polygons        |
# |============================================================|
# Returns a named list: poly_name -> data.table (one row per FID)

# f_nc <- nc_files[1]
# poly_meta_local = poly_meta_todo

# future_lapply(nc_files, function(f_nc) {
#   process_day(f_nc, poly_meta_todo)
# }, future.seed = TRUE)


process_day <- function(f_nc, poly_meta_local) {
  
  library(terra)
  library(data.table)
  library(ncdf4)
  library(sf)
  library(lubridate)
  
  # Extract date from filename
  bname    <- basename(f_nc)
  date_str <- substr(bname, 1, 8)
  date_i   <- tryCatch(as.Date(date_str, format = "%Y%m%d"),
                       error = function(e) NA_Date_)
  if (is.na(date_i)) return(NULL)
  
  # Load raster ONCE for this day
  r_all <- tryCatch(rast(f_nc), error = function(e) NULL)
  if (is.null(r_all)) return(NULL)
  
  sst_idx <- grep("analysed_sst$", names(r_all), value = FALSE)
  if (length(sst_idx) == 0) sst_idx <- 1L
  r_sst <- r_all[[sst_idx[1]]]
  
  ice_idx <- grep("sea_ice_fraction", names(r_all), value = FALSE)
  has_ice <- length(ice_idx) > 0
  r_ice   <- if (has_ice) r_all[[ice_idx[1]]] else NULL
  
  # Fix extent using true coordinates
  nc_tmp  <- nc_open(f_nc)
  lons_nc <- ncvar_get(nc_tmp, "lon")
  lats_nc <- ncvar_get(nc_tmp, "lat")
  nc_close(nc_tmp)
  
  true_ext <- ext(min(lons_nc) - 0.025, max(lons_nc) + 0.025,
                  min(lats_nc) - 0.025, max(lats_nc) + 0.025)
  
  ext(r_sst) <- true_ext
  crs(r_sst) <- "EPSG:4326"
  r_sst_c    <- r_sst - 273.15   # Kelvin -> Celsius
  
  if (has_ice) {
    ext(r_ice) <- true_ext
    crs(r_ice) <- "EPSG:4326"
  }
  
  # Process each polygon
  results <- lapply(poly_meta_local, function(pm) {
    
    rows_poly <- lapply(pm$fid_list, function(fid_item) {
      
      shp_fid <- tryCatch(
        vect(st_as_sf(st_as_sfc(fid_item$wkt, crs = 4326))),
        error = function(e) NULL
      )
      if (is.null(shp_fid)) return(NULL)
      
      # Crop + mask SST
      r_sst_m <- tryCatch(
        mask(crop(r_sst_c, shp_fid), shp_fid),
        error = function(e) NULL
      )
      if (is.null(r_sst_m)) return(NULL)
      
      sst_vals <- values(r_sst_m)[, 1]
      
      # Ice
      if (has_ice) {
        r_ice_m  <- tryCatch(
          mask(crop(r_ice, shp_fid), shp_fid),
          error = function(e) NULL
        )
        ice_vals <- if (!is.null(r_ice_m)) {
          v <- values(r_ice_m)[, 1]; v[is.na(v)] <- 0; v
        } else rep(0, length(sst_vals))
      } else {
        ice_vals <- rep(0, length(sst_vals))
      }
      
      # Build pixel table
      xy    <- xyFromCell(r_sst_m, seq_len(ncell(r_sst_m)))
      dt_px <- data.table(x = xy[,1], y = xy[,2],
                          sst = sst_vals, ice_frac = ice_vals)
      dt_px <- dt_px[!is.na(sst)]
      if (nrow(dt_px) == 0) return(NULL)
      
      # TH indices (Earth curvature correction inside thermal_habitat via cos(lat))
      res_1sw <- tryCatch(
        thermal_habitat(sst=dt_px$sst, ice=dt_px$ice_frac, year_sw=1L,
                        lon=dt_px$x, lat=dt_px$y, res_x=0.05, res_y=0.05),
        error = function(e) NULL
      )
      res_2sw <- tryCatch(
        thermal_habitat(sst=dt_px$sst, ice=dt_px$ice_frac, year_sw=2L,
                        lon=dt_px$x, lat=dt_px$y, res_x=0.05, res_y=0.05),
        error = function(e) NULL
      )
      
      data.table(
        date        = date_i,
        year        = year(date_i),
        month       = month(date_i),
        doy         = yday(date_i),
        FID         = as.character(fid_item$fid),
        n_pixels    = nrow(dt_px),
        TH_1SW_sum  = if (!is.null(res_1sw)) sum(res_1sw$growth_weighted,  na.rm=TRUE) else NA_real_,
        TH_1SW_mean = if (!is.null(res_1sw)) mean(res_1sw$growth_weighted, na.rm=TRUE) else NA_real_,
        TH_2SW_sum  = if (!is.null(res_2sw)) sum(res_2sw$growth_weighted,  na.rm=TRUE) else NA_real_,
        TH_2SW_mean = if (!is.null(res_2sw)) mean(res_2sw$growth_weighted, na.rm=TRUE) else NA_real_,
        sst_mean    = mean(dt_px$sst, na.rm=TRUE),
        sst_median  = median(dt_px$sst, na.rm=TRUE),
        sst_sd      = sd(dt_px$sst, na.rm=TRUE),
        sst_min     = min(dt_px$sst, na.rm=TRUE),
        sst_max     = max(dt_px$sst, na.rm=TRUE),
        ice_mean    = mean(ice_vals[!is.na(sst_vals)], na.rm=TRUE)
      )
    })
    
    rbindlist(rows_poly, fill = TRUE)
  })
  
  results   # named list: poly_name -> data.table
}


# |============================================================|
# |     MAIN LOOP: month -> days (parallel) -> polygons       |
# |============================================================|

plan(multisession, workers = n_workers)

# iz = 1
for (iz in seq_len(nrow(dt_zips))) {
  
  yr    <- dt_zips$year[iz]
  mo    <- dt_zips$month[iz]
  zpath <- dt_zips$path[iz]
  
  # Check if ALL polygons already have a cache for this month
  all_cached <- all(sapply(names(poly_meta), function(pname) {
    f <- file.path(poly_meta[[pname]]$out_dir,
                   paste0("TH_", yr, "_", sprintf("%02d", mo), ".fst"))
    file.exists(f) && file.size(f) > min_cache_size
  }))
  
  if (all_cached) {
    message(sprintf("[%d/%d] %d-%02d : all polygons cached, skipping.",
                    iz, nrow(dt_zips), yr, mo))
    next
  }
  
  message(sprintf("[%d/%d] %d-%02d : unzipping...",
                  iz, nrow(dt_zips), yr, mo))
  
  # Unzip monthly file ONCE
  tmp_month <- file.path(tmp_dir, paste0(yr, "_", sprintf("%02d", mo)))
  if (!dir.exists(tmp_month)) dir.create(tmp_month, recursive = TRUE)
  
  nc_names <- tryCatch(
    unzip(zpath, list = TRUE)$Name,
    error = function(e) { message("  ERROR listing ZIP: ", conditionMessage(e)); character(0) }
  )
  
  if (length(nc_names) == 0) {
    message("  Skipping — empty or corrupt ZIP.")
    unlink(tmp_month, recursive = TRUE)
    next
  }
  
  tryCatch(
    unzip(zpath, files = nc_names, exdir = tmp_month, overwrite = TRUE),
    error = function(e) message("  ERROR unzipping: ", conditionMessage(e))
  )
  
  nc_files <- file.path(tmp_month, nc_names)
  nc_files <- nc_files[file.exists(nc_files)]
  
  if (length(nc_files) == 0) {
    message("  No .nc files found after unzip.")
    unlink(tmp_month, recursive = TRUE)
    next
  }
  
  message(sprintf("  %d daily files | processing all polygons in parallel...",
                  length(nc_files)))
  
  # Identify which polygons still need processing this month
  poly_todo <- names(poly_meta)[sapply(names(poly_meta), function(pname) {
    f <- file.path(poly_meta[[pname]]$out_dir,
                   paste0("TH_", yr, "_", sprintf("%02d", mo), ".fst"))
    !file.exists(f) || file.size(f) <= min_cache_size
  })]
  
  # Subset poly_meta to only those needed
  poly_meta_todo <- poly_meta[poly_todo]
  
  # Process all daily files in parallel
  # Each worker returns a named list (poly -> data.table) for one day
  daily_list <- future_lapply(nc_files, function(f_nc) {
    process_day(f_nc, poly_meta_todo)
  }, future.seed = TRUE)
  
  # Assemble and save per polygon
  for (pname in poly_todo) {
    
    pm <- poly_meta[[pname]]
    f_cache <- file.path(pm$out_dir,
                         paste0("TH_", yr, "_", sprintf("%02d", mo), ".fst"))
    
    # Collect all days for this polygon
    poly_days <- rbindlist(
      lapply(daily_list, function(day_res) {
        if (is.null(day_res) || is.null(day_res[[pname]])) return(NULL)
        day_res[[pname]]
      }),
      fill = TRUE
    )
    
    if (nrow(poly_days) == 0) {
      message(sprintf("  [%s] No valid data for %d-%02d.", pname, yr, mo))
      next
    }
    
    # Restore original FID type
    if (!pm$fid_is_char) {
      poly_days[, FID := suppressWarnings(as.integer(FID))]
    }
    
    setorder(poly_days, date, FID)
    write_fst(poly_days, f_cache, compress = 50)
    message(sprintf("  [%s] %d rows saved (%d days x %d FIDs).",
                    pname,
                    nrow(poly_days),
                    uniqueN(poly_days$date),
                    uniqueN(poly_days$FID)))
  }
  
  # Delete unzipped files — ALWAYS keep original ZIPs
  unlink(tmp_month, recursive = TRUE)
  gc()
  
  message(sprintf("[%d/%d] %d-%02d done.", iz, nrow(dt_zips), yr, mo))
}

plan(sequential)





