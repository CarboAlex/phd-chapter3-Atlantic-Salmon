

# ══════════════════════════════════════════════════════════════════════════════
# CREATE MIGRATION CORRIDOR SEGMENTS — St. Jean River
#
# Divides the St. Jean migration corridor polygon into evenly spaced perpendicular
# segments using coast and open-sea reference lines. The segmented corridor is
# used downstream to compute segment-level thermal habitat indices.
#
# Inputs  : R/00_raw_data/ReferenceMigrationLayers/StJeanMigrationCorridor.shp
# Outputs : R/00_raw_data/ReferenceMigrationLayers/StJeanMigrationCorridor_Parts.shp
# ══════════════════════════════════════════════════════════════════════════════


library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(nngeo)
library(units)
library(lwgeom)


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

# ── World basemap ─────────────────────────────────────────────────────────────
# SOURCE: rnaturalearth R package (Natural Earth, medium scale)
# Filtered to North America only to reduce object size.
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(name_en %in% c("United States of America", "Canada", "Greenland"))

# ── St. Jean migration corridor (unsegmented) ─────────────────────────────────
# SOURCE: Migration corridor polygon delineated manually in GIS
# STATUS: Included in this repository
stj_cor <- st_read(file.path(base_path, "R/00_raw_data/ReferenceMigrationLayers/StJeanMigrationCorridor.shp"))



# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

# Coordinate reference systems
wgs84             <- 4326   # WGS84 geographic
nad83_can_lambert <- 6623   # NAD83(CSRS) / Statistics Canada Lambert — equal-area
# conic projection suited for national-scale analyses in Canada

# Bounding box of the study area (WGS84 decimal degrees)
x_west_ns  <- -70
x_east_ns  <- -52
y_north_ns <-  54
y_south_ns <-  46

# Number of segments to divide the corridor into
n_segments <- 5



# ══════════════════════════════════════════════════════════════════════════════
# FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

# ── create_migration_corridor() ───────────────────────────────────────────────
# Purpose    : Divides a migration corridor polygon into evenly spaced perpendicular
#              segments using coast and open-sea reference lines as guides.
# Parameters :
#   corridor_polygon   sf POLYGON    Original corridor polygon to split
#   coast_lines        sf LINESTRING Coastward edge of the corridor
#   sea_lines          sf LINESTRING Seaward edge of the corridor
#   number_of_squares  integer       Number of segments desired
#   crs                integer       CRS for intermediate operations
#                                    (default: 32198, NAD83 / Quebec Lambert)
#   extension_distance numeric       Distance (m) to extend cutting lines beyond
#                                    corridor edges to ensure clean splits (default: 1e6)
# Returns    : sf object with one polygon per segment; column FID numbers them 1-n.

create_migration_corridor <- function(corridor_polygon, coast_lines, sea_lines,
                                      number_of_squares, crs = 32198,
                                      extension_distance = 1e6) {
  
  # Sample evenly spaced points along the coastline reference line
  desired_segment_length <- as.numeric(st_length(coast_lines) / number_of_squares)
  coast_line_length      <- as.numeric(st_length(coast_lines))
  positions              <- seq(0, coast_line_length, by = desired_segment_length)
  fractions              <- positions / coast_line_length
  coast_points           <- st_line_sample(coast_lines, sample = fractions)
  coast_points           <- st_sf(geometry = coast_points)
  
  # For each coast sample point, find its nearest counterpart on the sea line
  nearest_point_list <- list()
  coast_points_coord <- st_coordinates(coast_points)[, c(1, 2)]
  
  for (i in seq_along(coast_points_coord[, 1])) {
    temp_point <- st_sfc(st_point(c(coast_points_coord[i, 1], coast_points_coord[i, 2])),
                         crs = st_crs(coast_lines))
    temp_point <- st_transform(temp_point, st_crs(sea_lines))
    nearest_point_list[[i]] <- st_nearest_points(temp_point, sea_lines)
  }
  
  # Build connecting lines between each coast-sea point pair
  geom_list        <- lapply(nearest_point_list, function(x) st_geometry(x)[[1]])
  connecting_lines <- st_sf(geometry = st_sfc(geom_list, crs = st_crs(nearest_point_list[[1]])))
  
  # Cast connecting lines into individual endpoints for use as cutting line anchors
  points_ns <- st_cast(connecting_lines, "POINT")
  
  # Internal helper — extend a line segment beyond the corridor to ensure a clean split
  extend_line <- function(p1, p2, crs, extension_distance) {
    direction <- p2 - p1
    unit_dir  <- direction / sqrt(sum(direction^2))
    new_p1    <- p1 - unit_dir * extension_distance
    new_p2    <- p2 + unit_dir * extension_distance
    st_sfc(st_linestring(rbind(new_p1, new_p2)), crs = crs)
  }
  
  # Iteratively split the corridor using each cutting line
  sequence_split    <- seq(3, nrow(points_ns), 2)
  remaining_polygon <- corridor_polygon
  corridor_pieces   <- list()
  
  for (i in seq_along(sequence_split)) {
    p1       <- st_coordinates(points_ns[sequence_split[i],     ])
    p2       <- st_coordinates(points_ns[sequence_split[i] + 1, ])
    cut_line <- extend_line(p1, p2,
                            crs                = st_crs(remaining_polygon),
                            extension_distance = extension_distance)
    
    polygon_split <- st_split(remaining_polygon, cut_line)
    pieces        <- st_collection_extract(polygon_split, "POLYGON")
    pieces_sf     <- st_sf(geometry = pieces)
    
    if (i != length(sequence_split)) {
      corridor_pieces[[i]] <- pieces_sf[2, ]   # completed segment
      remaining_polygon    <- pieces_sf[1, ]   # portion still to be split
    } else {
      corridor_pieces[[i]] <- remaining_polygon   # final segment
    }
  }
  
  # Assemble all segments and number them sequentially
  corridor_final     <- do.call(rbind, corridor_pieces)
  corridor_final$FID <- 1:nrow(corridor_final)
  
  return(corridor_final)
}


# ── convert_bbox_to_nad83() ───────────────────────────────────────────────────
# Purpose    : Converts a WGS84 bounding box into a target projected CRS.
#              Used to set ggplot2 map extents in projected coordinates.
# Parameters :
#   x_west, y_north, x_east, y_south  numeric  Corner coordinates in WGS84
#   crs_target                         integer  Target CRS (default: nad83_can_lambert)
# Returns    : Named list with x_west, x_east, y_north, y_south in the target CRS.

convert_bbox_to_nad83 <- function(x_west, y_north, x_east, y_south,
                                  crs_target = nad83_can_lambert) {
  
  point_upper_left_wgs84  <- st_sfc(st_point(c(x_west, y_north)), crs = wgs84)
  point_upper_left_nad83  <- st_transform(point_upper_left_wgs84,  crs = crs_target)
  
  point_lower_right_wgs84 <- st_sfc(st_point(c(x_east, y_south)), crs = wgs84)
  point_lower_right_nad83 <- st_transform(point_lower_right_wgs84, crs = crs_target)
  
  list(
    x_west_nad83  = st_coordinates(point_upper_left_nad83)[1],    # X of upper-left corner
    x_east_nad83  = st_coordinates(point_lower_right_nad83)[1],   # X of lower-right corner
    y_north_nad83 = st_coordinates(point_upper_left_nad83)[2],    # Y of upper-left corner
    y_south_nad83 = st_coordinates(point_lower_right_nad83)[2]    # Y of lower-right corner
  )
}



# ══════════════════════════════════════════════════════════════════════════════
# MAIN SCRIPT
# ══════════════════════════════════════════════════════════════════════════════

# ── Convert study area bounding box to NAD83 ──────────────────────────────────
ns_bbox_nad83 <- convert_bbox_to_nad83(
  x_west  = x_west_ns - 1,   # slight westward offset to avoid clipping the western edge
  y_north = y_north_ns,
  x_east  = x_east_ns,
  y_south = y_south_ns
)
x_west_ns_nad83  <- ns_bbox_nad83$x_west_nad83
x_east_ns_nad83  <- ns_bbox_nad83$x_east_nad83
y_north_ns_nad83 <- ns_bbox_nad83$y_north_nad83
y_south_ns_nad83 <- ns_bbox_nad83$y_south_nad83


# ── Crop and reproject basemap ────────────────────────────────────────────────
# Crop to the study area extent to reduce processing time, then reproject to NAD83.
world_original <- world
bbox_coords    <- st_bbox(c(xmin = -75, xmax = -40, ymin = 40, ymax = 58), crs = wgs84)
bbox_poly      <- st_as_sfc(bbox_coords)
world          <- st_intersection(world, bbox_poly)

world_nad83          <- st_transform(world,          crs = nad83_can_lambert)
world_original_nad83 <- st_transform(world_original, crs = nad83_can_lambert)


# ── Split Canada into individual polygons ─────────────────────────────────────
# A buffer around the Quebec-Labrador polygon will later define the corridor
# width. Casting the multi-polygon into parts allows targeting that polygon
# specifically by index.
world_parts          <- st_cast(world_nad83, "POLYGON")
world_parts$ID       <- 1:nrow(world_parts)
world_parts$Location <- NA
world_parts$Location[10] <- "Quebec_Labrador"
world_parts$Location[12] <- "Gaspesie"
world_parts$Location[13] <- "New_Foundland"
world_parts$Location[14] <- "Anticosti"

# Visual check — confirm polygon assignments are correct
ggplot() +
  geom_sf(data = world_parts, aes(fill = Location))


# ── Extract coast and sea reference lines ─────────────────────────────────────
# Decompose the corridor polygon boundary into individual segments, then
# classify them as belonging to the coastward or seaward side of the corridor.

lines_ns <- st_segments(stj_cor)

# Remove the two short closing segments at the polygon ends; these connect the
# coast and sea sides and would interfere with perpendicular splitting.
lines_ns <- lines_ns[-c(8, 49), ]
print("Note that these numbers may have to change if the corridor changes.")

# Visual check — all remaining segments should outline the corridor sides
ggplot() +
  geom_sf(data = lines_ns, color = "red", size = 1, alpha = 0.4) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white", colour = "white", linetype = "solid"),
        panel.grid.major = element_blank())

coast_lines <- lines_ns[c(1:7, 48:49), ]
sea_lines   <- lines_ns[c(8:47),       ]

# Visual check — red = coast side, blue = open-sea side
ggplot() +
  geom_sf(data = coast_lines, color = "red",  size = 2) +
  geom_sf(data = sea_lines,   color = "blue", size = 2) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white", colour = "white", linetype = "solid"),
        panel.grid.major = element_blank())

# Merge individual segments into single continuous lines for each side,
# as the splitting algorithm requires a single linestring per reference
coast_lines <- st_line_merge(st_union(coast_lines))
sea_lines   <- st_line_merge(st_union(sea_lines))

# Ensure both lines run in the same direction, as required by the splitting algorithm
if (st_coordinates(coast_lines)[nrow(st_coordinates(coast_lines)), 2] == st_coordinates(sea_lines)[1, 2]) {
  sea_lines_coord      <- st_coordinates(sea_lines)
  sea_lines_coord[, 1] <- rev(sea_lines_coord[, 1])
  sea_lines_coord[, 2] <- rev(sea_lines_coord[, 2])
  sea_lines <- st_sfc(st_linestring(sea_lines_coord[, 1:2], dim = "XY"), crs = nad83_can_lambert)
}


# ── Segment the migration corridor ────────────────────────────────────────────
stj_cor_all <- create_migration_corridor(
  corridor_polygon   = stj_cor,
  coast_lines        = coast_lines,
  sea_lines          = sea_lines,
  number_of_squares  = n_segments,
  crs                = 32198,
  extension_distance = 1e6
)

# Visual check — each segment should be a distinct, contiguous polygon
ggplot() +
  geom_sf(data = world_nad83) +
  geom_sf(data = stj_cor_all, aes(fill = factor(FID)), color = "black", size = 1, alpha = 0.4) +
  coord_sf(xlim = c(x_west_ns_nad83, x_east_ns_nad83),
           ylim = c(y_south_ns_nad83, y_north_ns_nad83),
           expand = FALSE) +
  labs(x = "Longitude", y = "Latitude", fill = "Block") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white", colour = "white", linetype = "solid"),
        panel.grid.major = element_blank())


# ── Export ────────────────────────────────────────────────────────────────────
st_write(stj_cor_all, file.path(base_path, "R/00_raw_data/ReferenceMigrationLayers/StJeanMigrationCorridor_Parts.shp"), 
         delete_layer = TRUE)


