


# ==============================================================================
# Download daily ESA CCI Level-4 SST data for the Labrador Sea region from
# the Copernicus Climate Data Store (CDS) API. One NetCDF file is saved per
# month, covering the period 1982-2024.
#
# Inputs:
#   - CDS API token stored in .Renviron at the repository root (CDS_TOKEN=...)
#
# Outputs:
#   - R/00_raw_data/SST_ORAS5_Labrador/ESACCI_SST_L4_<YYYY>_<MM>.nc
#     One NetCDF file per month (~516 files for 1982-2024).
# ==============================================================================


# ── Working directory ──────────────────────────────────────────────────────
# Automatically detects the repository root based on the script location.
# No hardcoded path — works on any machine.
current_file <- rstudioapi::getActiveDocumentContext()$path
path_parts   <- unlist(strsplit(normalizePath(current_file, winslash = "/"), "/"))
target_index <- which(path_parts == "phd-chapter3-Atlantic-Salmon")
base_path    <- paste(path_parts[1:target_index], collapse = "/")
setwd(base_path)


library(ecmwfr)


# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

# ── Authentication ─────────────────────────────────────────────────────────────
# This script requires a Copernicus Climate Data Store (CDS) API token.
# To set it up:
#   1. Register at https://cds.climate.copernicus.eu and retrieve your token
#      from your profile page.
#   2. Create a file named .Renviron at the root of this repository
#      (phd-chapter3-Atlantic-Salmon/.Renviron) with the following line:
#         CDS_TOKEN=your_token_here
# The .Renviron file is listed in .gitignore and will not be committed.

readRenviron(".Renviron")

download_dir <- file.path(base_path, "R/00_raw_data/SST_ORAS5_Labrador")  # Output directory for downloaded files

years_download  <- 1982:2023  # Years to download
months_download <- 1:12       # Months to download (all months)

# Geographic bounding box — Labrador Sea with margin
# CDS format: c(North, West, South, East)
area_labrador <- c(65, -68, 43, -38)

# Pause between requests (seconds) — avoids overloading the server
sleep_between_requests <- 30


# ══════════════════════════════════════════════════════════════════════════════
# MAIN SCRIPT
# ══════════════════════════════════════════════════════════════════════════════

# ── Initialize output directory and CDS credentials ──────────────────────────

if (!dir.exists(download_dir)) dir.create(download_dir, recursive = TRUE)

wf_set_key(key = Sys.getenv("CDS_TOKEN"))


# ── Build list of months to download ─────────────────────────────────────────

# One file per month (~516 files for 1982-2024); each file contains all days
# of that month. Files already on disk (> 1 MB) are skipped automatically.

months_todo <- expand.grid(year = years_download, month = months_download)
months_todo <- months_todo[order(months_todo$year, months_todo$month), ]

n_total   <- nrow(months_todo)
n_done    <- 0
n_skipped <- 0
n_error   <- 0


# ── Download loop — one request per month ─────────────────────────────────────

for (i in seq_len(n_total)) {
  
  yr <- months_todo$year[i]
  mo <- months_todo$month[i]
  
  fname <- paste0("ESACCI_SST_L4_", yr, "_", sprintf("%02d", mo), ".nc")
  f_out <- file.path(download_dir, fname)
  
  # Skip if file already exists and is not empty (> 1 MB)
  if (file.exists(f_out) && file.size(f_out) > 1e6) {
    n_skipped <- n_skipped + 1
    if (n_skipped %% 12 == 0)
      message("Skipped ", n_skipped, " already-downloaded months...")
    next
  }
  
  message(sprintf("[%d/%d] Downloading %d-%02d...", i, n_total, yr, mo))
  
  # Request all possible days; the server ignores days that do not exist
  all_days <- sprintf("%02d", 1:31)
  
  request <- list(
    dataset_short_name  = "satellite-sea-surface-temperature",
    variable            = "all",
    processinglevel     = "level_4",
    sensor_on_satellite = "combined_product",
    version             = "3_0",
    temporal_resolution = "daily",
    year                = as.character(yr),
    month               = ifelse(nchar(mo) == 2, as.character(mo), as.character(paste0(0, mo))),
    day                 = sprintf("%02d", 1:31),
    data_format         = "netcdf",
    download_format     = "unarchived",
    area                = c(65, -68, 43, -38)
  )
  
  tryCatch({
    
    # Retry up to 3 times if a rate-limit error is encountered
    success  <- FALSE
    attempts <- 0
    
    while (!success && attempts < 3) {
      attempts <- attempts + 1
      tryCatch({
        wf_request(
          request  = request,
          transfer = TRUE,
          path     = download_dir,
          verbose  = FALSE
        )
        
        # Find the most recently created file in the directory and rename it
        new_files <- list.files(download_dir, full.names = TRUE)
        new_files <- new_files[!file.exists(paste0(new_files, ".lock"))]
        latest    <- new_files[which.max(file.mtime(new_files))]
        
        if (file.exists(latest) && latest != f_out) {
          file.rename(latest, f_out)
        }
        
        success <- TRUE
      }, error = function(e) {
        msg <- conditionMessage(e)
        if (grepl("rate limit|429|Rate limit", msg, ignore.case = TRUE)) {
          wait_time <- 60 * attempts  # 60s, 120s, 180s for successive retries
          message(sprintf("  -> Rate limit hit. Waiting %ds before retry %d/3...",
                          wait_time, attempts))
          Sys.sleep(wait_time)
        } else {
          stop(e)  # Propagate non-rate-limit errors immediately
        }
      })
    }
    
    if (!success) stop("Max retries reached")
    
    n_done <- n_done + 1
    message(sprintf("  -> Done. (%d downloaded, %d skipped, %d errors)",
                    n_done, n_skipped, n_error))
    Sys.sleep(60)  # 60-second pause between requests
    
  }, error = function(e) {
    n_error <<- n_error + 1
    message(sprintf("  -> ERROR : %s", conditionMessage(e)))
    if (file.exists(f_out) && file.size(f_out) < 1e6) file.remove(f_out)
    Sys.sleep(60)
  })
}


# ── Verify downloaded files ───────────────────────────────────────────────────

all_files <- list.files(download_dir, pattern = "\\.nc$", full.names = TRUE)
message("Files on disk : ", length(all_files))

# Identify any missing months
expected <- paste0(
  "ESACCI_SST_L4_",
  rep(years_download, each = length(months_download)), "_",
  sprintf("%02d", rep(months_download, length(years_download))),
  ".nc"
)
missing <- expected[!expected %in% basename(all_files)]
if (length(missing) > 0) {
  message("Missing months (", length(missing), ") :")
  print(head(missing, 20))
} else {
  message("All months present — download complete!")
}


