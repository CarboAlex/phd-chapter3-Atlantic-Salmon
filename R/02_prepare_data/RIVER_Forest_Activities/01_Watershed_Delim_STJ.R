

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: HEADER
# ══════════════════════════════════════════════════════════════════════════════
#
# Downloads a DEM and delineates the watershed of the Saint-Jean River using
# hydrological preprocessing (depression filling, D8 flow direction and
# accumulation, outlet snapping) via WhiteboxTools.
#
# Inputs:
#   - Internet connection (DEM fetched via elevatr; river line via osmdata)
#
# Outputs:
#   - R/01_derived_data/ForestActivities/dem_utm_STJ.tif
#       DEM reprojected to EPSG:32198
#   - R/01_derived_data/ForestActivities/dem_filled_STJ.tif
#       Depression-filled DEM
#   - R/01_derived_data/ForestActivities/d8_pointer_STJ.tif
#       D8 flow direction raster
#   - R/01_derived_data/ForestActivities/flow_accum_STJ.tif
#       D8 flow accumulation raster (cells)
#   - R/01_derived_data/ForestActivities/outlet_STJ.shp
#       Lowest point along the river line (pour point)
#   - R/01_derived_data/ForestActivities/snapped_outlet_STJ.shp
#       Pour point snapped to the highest-accumulation cell within 200 m
#   - R/01_derived_data/ForestActivities/raster_watershed_STJ.tif
#       Watershed raster delineated from the snapped pour point
#   - R/01_derived_data/ForestActivities/poly_watershed_STJ.shp
#       Watershed polygon (dissolved, valid geometry)
#
# ══════════════════════════════════════════════════════════════════════════════


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: LIBRARIES
# ══════════════════════════════════════════════════════════════════════════════

library(sf)
library(terra)
library(elevatr)
library(osmdata)
library(dplyr)
library(whitebox)


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
# SECTION 5: CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

out_dir   <- "R/01_derived_data/ForestActivities"  # Output directory
out_river <- "STJ"                                  # River identifier tag

# Bounding box covering the Saint-Jean watershed (WGS84)
bbx <- c(
  xmin = -65.683279,
  ymin =  48.5,
  xmax = -64.434943,
  ymax =  48.946834
)

snap_dist <- 200  # Snapping distance (m) to align pour point with flow accumulation


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7: MAIN SCRIPT
# ══════════════════════════════════════════════════════════════════════════════

# ── Build reference area polygon ──────────────────────────────────────────────

bbox_sf <- st_as_sfc(st_bbox(bbx, crs = 4326))
bbox_sf <- st_as_sf(bbox_sf)


# ── Download and reproject DEM ────────────────────────────────────────────────

# Fetch DEM at zoom level 12 clipped to the bounding box
dem_wgs <- get_elev_raster(
  locations = bbox_sf,
  z         = 12,
  clip      = "bbox"
)

plot(dem_wgs)

# Convert to terra SpatRaster and reproject to Quebec Lambert (EPSG:32198)
dem <- rast(dem_wgs)
dem <- project(dem, "EPSG:32198")


# ── Download river line from OpenStreetMap ────────────────────────────────────

river_osm <- opq(bbox = bbx, timeout = 200) |>
  add_osm_feature(key = "waterway", value = "river") |>
  osmdata_sf()

# OSM stores names with accents; iconv strips them for ASCII-safe comparison
river_line <- river_osm$osm_lines[
  which(iconv(river_osm$osm_lines$name, to = "ASCII//TRANSLIT") ==
          "Riviere Saint-Jean"), ]
river_line <- st_transform(river_line, crs(dem))

plot(dem)
plot(river_line, add = TRUE, col = "red", lwd = 2)


# ── Hydrological preprocessing ────────────────────────────────────────────────

# Write the projected DEM to disk — required as file input by WhiteboxTools
mnt_file <- file.path(out_dir, paste0("dem_utm_", out_river, ".tif"))
writeRaster(dem, mnt_file, overwrite = TRUE)

# Fill sinks so that all cells drain to the outlet without interruption
dem_filled <- file.path(out_dir, paste0("dem_filled_", out_river, ".tif"))
wbt_fill_depressions(dem = mnt_file, output = dem_filled)

# Compute D8 flow direction and accumulation from the filled DEM
d8_pointer <- file.path(out_dir, paste0("d8_pointer_", out_river, ".tif"))
flow_acc   <- file.path(out_dir, paste0("flow_accum_",  out_river, ".tif"))
wbt_d8_pointer(dem = dem_filled, output = d8_pointer)
wbt_d8_flow_accumulation(
  input    = dem_filled,
  output   = flow_acc,
  out_type = "cells"
)


# ── Identify the lowest point on the river (pour point) ───────────────────────

# Sample points along the river line at 50 m intervals
river_pts      <- st_line_sample(river_line, density = 1 / 50) |>
  st_cast("POINT")

# Extract DEM elevation at each point and select the lowest one as outlet
river_pts_vect <- vect(river_pts)
elev           <- terra::extract(dem, river_pts_vect)[, 2]
outlet         <- river_pts[which.min(elev), ]

plot(dem)
plot(river_line, add = TRUE, col = "blue", lwd = 2)
plot(outlet,     add = TRUE, col = "red",  pch = 19, cex = 1.5)

# Export outlet as shapefile for use with WhiteboxTools
outlet_file <- file.path(out_dir, paste0("outlet_", out_river, ".shp"))
st_write(outlet, outlet_file, delete_layer = TRUE)

# Snap the pour point to the nearest high-accumulation cell within snap_dist
snap_point_file <- file.path(out_dir, paste0("snapped_outlet_", out_river, ".shp"))
wbt_snap_pour_points(
  pour_pts   = outlet_file,
  flow_accum = flow_acc,
  output     = snap_point_file,
  snap_dist  = snap_dist
)


# ── Delineate watershed ───────────────────────────────────────────────────────

watershed_raster <- file.path(out_dir, paste0("raster_watershed_", out_river, ".tif"))
wbt_watershed(
  d8_pntr  = d8_pointer,
  pour_pts = snap_point_file,
  output   = watershed_raster
)

# Convert raster watershed to a dissolved, valid polygon
basin_r    <- rast(watershed_raster)
basin_poly <- as.polygons(basin_r, dissolve = TRUE)
basin_sf   <- st_as_sf(basin_poly)
basin_sf   <- st_make_valid(basin_sf)

plot(dem,      col = hcl.colors(50, "Terrain"), alpha = 0.7)
plot(basin_sf, add = TRUE,
     col    = adjustcolor("#e31a1c", alpha.f = 0.3),
     border = "#e31a1c",
     lwd    = 2.5)
plot(river_line, add = TRUE, col = "blue", lwd = 2.5)


# ── Export watershed polygon ──────────────────────────────────────────────────

st_write(basin_sf,
         file.path(out_dir, paste0("poly_watershed_", out_river, ".shp")))

