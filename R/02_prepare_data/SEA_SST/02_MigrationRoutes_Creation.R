

# =============================================================================
# Atlantic Salmon Marine Migration Corridors — Spatial Polygon Construction
# =============================================================================
# 
# Purpose : Delineate five spatial polygons representing the marine migration
#           route of Atlantic salmon post-smolts from river mouth to the
#           Labrador Sea. These polygons are used to extract satellite-derived
#           environmental variables (SST, zooplankton proxy) along the route.
#
# Polygons:
#   1. Trinité migration corridor (North Shore)
#   2. Saint-Jean migration corridor (Gulf of St. Lawrence)
#   3. Strait of Belle Isle (SOBI)
#   4. Labrador Sea entrance (post-Belle Isle zone)
#   5. Labrador Sea (offshore zone based on NAFO survey stations)
#
# Output  : Five shapefiles exported to R/00_raw_data/ReferenceMigrationLayers/
# =============================================================================


# =============================================================================
# Setup: working directory, libraries, base data
# =============================================================================

# ── Working directory ─────────────────────────────────────────────────────────
# Automatically detects the repository root based on the script location.
# No hardcoded path — works on any machine.
current_file <- rstudioapi::getActiveDocumentContext()$path
path_parts   <- unlist(strsplit(normalizePath(current_file, winslash = "/"), "/"))
target_index <- which(path_parts == "phd-chapter3-Atlantic-Salmon")
base_path    <- paste(path_parts[1:target_index], collapse = "/")
setwd(base_path)


# ── Libraries ─────────────────────────────────────────────────────────────────
library(tidyr)
library(dplyr)
library(ggplot2)
library(maps)
library(rnaturalearth)
library(sf)
library(nngeo)
library(smoothr)
library(ggpubr)
library(lwgeom)


# ── Coordinate reference systems ──────────────────────────────────────────────
wgs84           <- 4326
nad83_can_lambert <- 6623
# NAD83(CSRS) / Statistics Canada Lambert: conic equal-area projection designed
# for national-scale analyses in Canada. Preserves area accurately with low
# distortion across the entire country — well suited for ecological and spatial
# studies over large regions such as eastern Quebec and Labrador.


# ── Map extent limits (WGS84 decimal degrees) ────────────────────────────────
# Full study area (North Shore overview map)
x_west_ns  <- -70;  x_east_ns  <- -52
y_south_ns <- 46;   y_north_ns <- 54

# Inset: Trinité river mouth
x_west_ns_t  <- -68;  x_east_ns_t  <- -66
y_south_ns_t <- 49;   y_north_ns_t <- 50

# Inset: Saint-Jean river mouth
x_west_sj  <- -64.8;  x_east_sj  <- -64
y_south_sj <- 48.6;   y_north_sj <- 49


# ── Base map: land polygons ───────────────────────────────────────────────────
# Load and crop to the study region for efficiency
world_original <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(name_en %in% c("United States of America", "Canada", "Greenland"))

bbox_coords <- st_bbox(c(xmin = -75, xmax = -40, ymin = 40, ymax = 58), crs = wgs84)
bbox_poly   <- st_as_sfc(bbox_coords)
world       <- st_intersection(world_original, bbox_poly)

# Project to NAD83 Canada Lambert
world_nad83          <- st_transform(world,          crs = nad83_can_lambert)
world_original_nad83 <- st_transform(world_original, crs = nad83_can_lambert)

# Smoothed version: reduces artefacts from large estuaries / rivers in buffers
world_nad83_smooth <- smooth(world_nad83, method = "ksmooth", smoothness = 8)


# ── Convert map extent corners to NAD83 ──────────────────────────────────────
# Full study area
point_lim_upper_left_wgs84  <- st_sfc(st_point(c(x_west_ns - 1, y_north_ns)), crs = wgs84)
point_lim_upper_left_nad83  <- st_transform(point_lim_upper_left_wgs84,  crs = nad83_can_lambert)
point_lim_lower_right_wgs84 <- st_sfc(st_point(c(x_east_ns, y_south_ns)),    crs = wgs84)
point_lim_lower_right_nad83 <- st_transform(point_lim_lower_right_wgs84, crs = nad83_can_lambert)
x_west_ns_nad83  <- st_coordinates(point_lim_upper_left_nad83)[1]   # ≈ −162 692
x_east_ns_nad83  <- st_coordinates(point_lim_lower_right_nad83)[1]  # ≈ 1 267 072
y_north_ns_nad83 <- st_coordinates(point_lim_upper_left_nad83)[2]   # ≈ 1 118 953
y_south_ns_nad83 <- st_coordinates(point_lim_lower_right_nad83)[2]  # ≈ 367 087

# Inset: Trinité mouth point
point_lim_upper_left_wgs84_t  <- st_sfc(st_point(c(x_west_ns_t, y_north_ns_t)), crs = wgs84)
point_lim_upper_left_nad83_t  <- st_transform(point_lim_upper_left_wgs84_t,  crs = nad83_can_lambert)
point_lim_lower_right_wgs84_t <- st_sfc(st_point(c(x_east_ns_t, y_south_ns_t)), crs = wgs84)
point_lim_lower_right_nad83_t <- st_transform(point_lim_lower_right_wgs84_t, crs = nad83_can_lambert)
x_west_ns_nad83_t  <- st_coordinates(point_lim_upper_left_nad83_t)[1]
x_east_ns_nad83_t  <- st_coordinates(point_lim_lower_right_nad83_t)[1]
y_north_ns_nad83_t <- st_coordinates(point_lim_upper_left_nad83_t)[2]
y_south_ns_nad83_t <- st_coordinates(point_lim_lower_right_nad83_t)[2]

# Inset: Saint-Jean mouth point
point_lim_upper_left_wgs84_sj  <- st_sfc(st_point(c(x_west_sj, y_north_sj)), crs = wgs84)
point_lim_upper_left_nad83_sj  <- st_transform(point_lim_upper_left_wgs84_sj,  crs = nad83_can_lambert)
point_lim_lower_right_wgs84_sj <- st_sfc(st_point(c(x_east_sj, y_south_sj)), crs = wgs84)
point_lim_lower_right_nad83_sj <- st_transform(point_lim_lower_right_wgs84_sj, crs = nad83_can_lambert)
x_west_nad83_sj  <- st_coordinates(point_lim_upper_left_nad83_sj)[1]
x_east_nad83_sj  <- st_coordinates(point_lim_lower_right_nad83_sj)[1]
y_north_nad83_sj <- st_coordinates(point_lim_upper_left_nad83_sj)[2]
y_south_nad83_sj <- st_coordinates(point_lim_lower_right_nad83_sj)[2]


# ── Shared parameters ─────────────────────────────────────────────────────────
migration_route_width <- 100000  # Buffer width for migration corridors (m)

# Decompose world map into individual polygons (used throughout)
world_parts      <- st_cast(world_nad83, "POLYGON")
world_parts$ID   <- 1:nrow(world_parts)
# Key polygon IDs (verify after any map update):
#   10 = Quebec / Labrador mainland
#   12 = Gaspésie peninsula
#   13 = Newfoundland
#   14 = Anticosti Island

# Coastline of the Quebec/Labrador polygon as a line object (reused across sections)
world_parts_quebec_line <- st_line_merge(st_union(st_segments(world_parts[10, ])))

# Limit line between the two migration corridors and the Strait of Belle Isle
# (used to clip both the Trinité and Saint-Jean corridors at the same boundary)
limit_lab <- c(-57.146811, 51.423274)  # Quebec / Labrador side
limit_tn  <- c(-56.797615, 51.244339)  # Newfoundland side
limit_lab_tn_line <- st_transform(
  st_sfc(st_linestring(rbind(limit_lab, limit_tn)), crs = wgs84),
  crs = nad83_can_lambert
)

extend_line <- function(line_sf, extension_m = 10000) {
  # Extends a two-point line symmetrically by extension_m metres on each end.
  coords   <- st_coordinates(line_sf)
  if (nrow(coords) != 2) stop("The line must have exactly two points.")
  vec      <- coords[2, ] - coords[1, ]
  unit_vec <- vec / sqrt(sum(vec^2))
  new_coords <- rbind(coords[1, ] - extension_m * unit_vec,
                      coords[2, ] + extension_m * unit_vec)
  st_sfc(st_linestring(new_coords), crs = st_crs(line_sf))
}
limit_lab_tn_line <- extend_line(limit_lab_tn_line, extension_m = 100000)


# =============================================================================
# Trinité migration corridor (North Shore)
# =============================================================================
# Strategy: buffer the Quebec/Labrador coastline polygon outward by
# `migration_route_width`, then clip to retain only the nearshore band between
# the Trinité river mouth and the Strait of Belle Isle entrance.

# ── Key coordinates ───────────────────────────────────────────────────────────
lon_t <- -67.302472
lat_t <-  49.418111
distance_from_mouth_migration_route_start_t <- 5000  # Distance (m) from mouth where corridor starts

lon_ns_limit <- -56.830072  # North-eastern limit of corridor (Labrador coast)
lat_ns_limit <-  51.516840


# ── Step 1: Create initial 100 km coastal buffer ──────────────────────────────
area_ns <- st_buffer(world_parts[10, ], migration_route_width, endCapStyle = "SQUARE")
area_ns <- st_difference(area_ns, world_parts[10, ])  # Keep only the water-side band


# ── Step 2: Define SW and NE boundary polygons to clip the buffer ─────────────

# South-western boundary — near Trinité river mouth
river_mouth_t_wgs84   <- st_sfc(st_point(c(lon_t, lat_t)), crs = wgs84)
coords_mouth_t_nad83  <- st_transform(river_mouth_t_wgs84, crs = nad83_can_lambert)
coords_mouth_t_nad83  <- st_nearest_points(coords_mouth_t_nad83, world_parts_quebec_line)
mouth_on_shore_t      <- st_cast(coords_mouth_t_nad83, "POINT")[2]

# Find the corridor start point on the shore (5 km seaward from mouth)
circle_t                  <- st_buffer(mouth_on_shore_t, dist = distance_from_mouth_migration_route_start_t)
intersection_points_t     <- st_intersection(circle_t, world_parts_quebec_line)
start_corridor_ns_coord   <- st_coordinates(intersection_points_t)[
  which(st_coordinates(intersection_points_t)[, 2] == min(st_coordinates(intersection_points_t)[, 2])), 1:2
]
start_corridor_ns <- st_sfc(st_point(c(start_corridor_ns_coord[1], start_corridor_ns_coord[2])),
                            crs = nad83_can_lambert)

south_west_limit_ns_square <- matrix(c(
  start_corridor_ns_coord[1],                             start_corridor_ns_coord[2] - 5000,
  start_corridor_ns_coord[1] + migration_route_width + 100000, start_corridor_ns_coord[2] - 5000,
  start_corridor_ns_coord[1] + migration_route_width + 100000, start_corridor_ns_coord[2],
  start_corridor_ns_coord[1] - 10000,                    start_corridor_ns_coord[2],
  start_corridor_ns_coord[1],                             start_corridor_ns_coord[2] - 5000
), ncol = 2, byrow = TRUE)
south_west_limit_ns_square <- st_sfc(st_polygon(list(south_west_limit_ns_square)), crs = nad83_can_lambert)

# North-eastern boundary — near Labrador coast
lab_limit_wgs84   <- st_sfc(st_point(c(lon_ns_limit, lat_ns_limit)), crs = wgs84)
lab_limit_nad83   <- st_transform(lab_limit_wgs84, crs = nad83_can_lambert)
lab_limit_nad83   <- st_nearest_points(lab_limit_nad83, world_parts_quebec_line)
lab_limit_on_shore       <- st_cast(lab_limit_nad83, "POINT")[2]
lab_limit_on_shore_coord <- st_coordinates(lab_limit_on_shore)

north_east_limit_ns_square <- matrix(c(
  lab_limit_on_shore_coord[1],                             lab_limit_on_shore_coord[2] - 5000,
  lab_limit_on_shore_coord[1] + migration_route_width + 100000, lab_limit_on_shore_coord[2] - 5000,
  lab_limit_on_shore_coord[1] + migration_route_width + 100000, lab_limit_on_shore_coord[2],
  lab_limit_on_shore_coord[1] - 10000,                    lab_limit_on_shore_coord[2],
  lab_limit_on_shore_coord[1],                             lab_limit_on_shore_coord[2] - 5000
), ncol = 2, byrow = TRUE)
north_east_limit_ns_square <- st_sfc(st_polygon(list(north_east_limit_ns_square)), crs = nad83_can_lambert)


# ── Step 3: Build the boundary mask and clip the buffer ───────────────────────
temp_buffer_area_ns <- st_buffer(area_ns, 20000, endCapStyle = "SQUARE")
temp_buffer_area_ns <- st_difference(temp_buffer_area_ns, area_ns)

temp_buffer_area_ns <- st_union(rbind(
  st_as_sf(st_sf(geometry = st_geometry(north_east_limit_ns_square))),
  st_as_sf(st_sf(geometry = st_geometry(south_west_limit_ns_square))),
  st_as_sf(st_sf(geometry = st_geometry(temp_buffer_area_ns)))
))

area_ns <- st_difference(area_ns, temp_buffer_area_ns)


# ── Step 4: Isolate corridor polygon and remove land overlap ──────────────────
area_ns       <- st_sf(geometry = st_geometry(area_ns))
area_ns_parts <- st_cast(area_ns, "POLYGON")
area_ns_parts$ID <- 1:nrow(area_ns_parts)
area_ns_parts <- area_ns_parts[area_ns_parts$ID == 2, ]

area_ns_final <- st_difference(area_ns_parts, st_union(world_parts[c(10, 12, 13, 14), ]))


# ── Step 5: Separate coast and open-sea lines; crop south of Anticosti ────────
lines_ns <- st_segments(area_ns_final)

# Remove the two polygon end-cap segments (not useful for this step)
lines_ns <- lines_ns[-c(126, 186, 187, 188), ]
message("Note: line indices 126, 186-188 may need adjustment if the corridor geometry changes.")

lines_ns$Distance_From_Coast <- as.numeric(st_distance(lines_ns, world_parts[10, ]))
coast_lines <- lines_ns %>% filter(Distance_From_Coast <  500)
sea_lines   <- lines_ns %>% filter(Distance_From_Coast > 5000)

coast_lines <- st_line_merge(st_union(coast_lines))
sea_lines   <- st_line_merge(st_union(sea_lines))

# Polygon to crop the corridor south of Anticosti Island
sea_lines_coord   <- as.data.frame(st_coordinates(sea_lines))
pt_coord          <- sea_lines_coord[which(sea_lines_coord$X == min(sea_lines_coord$X)), ]
starting_point_crop_tr_cor <- st_as_sf(pt_coord, coords = c("X", "Y"), crs = nad83_can_lambert)

anticosti_west_point_google_map <- st_transform(
  st_sfc(st_point(c(-64.519987, 49.870548)), crs = wgs84), nad83_can_lambert
)
anticosti_west_point_nearest <- st_cast(
  st_nearest_points(anticosti_west_point_google_map, world_parts[14, ]), "POINT"
)[2]

reference_crop_tr_cor <- matrix(c(
  st_coordinates(starting_point_crop_tr_cor)[1],       st_coordinates(starting_point_crop_tr_cor)[2],
  st_coordinates(anticosti_west_point_nearest)[1],     st_coordinates(anticosti_west_point_nearest)[2],
  st_coordinates(st_transform(st_sfc(st_point(c(-62.014456, 49.235988)), crs = wgs84), nad83_can_lambert))[1],
  st_coordinates(st_transform(st_sfc(st_point(c(-62.014456, 49.235988)), crs = wgs84), nad83_can_lambert))[2],
  st_coordinates(st_transform(st_sfc(st_point(c(-65.078203, 48.695625)), crs = wgs84), nad83_can_lambert))[1],
  st_coordinates(st_transform(st_sfc(st_point(c(-65.078203, 48.695625)), crs = wgs84), nad83_can_lambert))[2],
  st_coordinates(starting_point_crop_tr_cor)[1],       st_coordinates(starting_point_crop_tr_cor)[2]
), ncol = 2, byrow = TRUE)
reference_crop_tr_cor <- st_sfc(st_polygon(list(reference_crop_tr_cor)), crs = nad83_can_lambert)

area_ns_final <- st_difference(area_ns_final, reference_crop_tr_cor)


# ── Step 6: Keep only the main corridor polygon; clip at Belle Isle line ───────
polygons_list     <- st_cast(area_ns_final, "POLYGON")
polygons_list$ID  <- 1:nrow(polygons_list)
tr_area_ns_final  <- polygons_list[2, ]

tr_area_split    <- st_split(tr_area_ns_final, limit_lab_tn_line)
tr_area_ns_final <- st_collection_extract(tr_area_split, "POLYGON")[1, ]


# ── Diagnostic plot ───────────────────────────────────────────────────────────
ggplot() +
  geom_sf(data = world_nad83) +
  geom_sf(data = tr_area_ns_final, fill = "blue", color = "blue", size = 1, alpha = 0.4) +
  coord_sf(xlim = c(-249419.3, x_east_ns_nad83),
           ylim = c(y_south_ns_nad83, y_north_ns_nad83), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_blank())


# =============================================================================
# Saint-Jean migration corridor (Gulf of St. Lawrence)
# =============================================================================
# Strategy: define a curved migration line from the Gaspésie coast tip, past
# Anticosti Island, to the Belle Isle limit. Buffer this line to form the
# corridor polygon, then remove land overlap.

# ── Key coordinates ───────────────────────────────────────────────────────────
lon_sj <- -64.381914
lat_sj <-  48.779793

lat_mouth_gaspe_bay       <- 48.696210;  lon_mouth_gaspe_bay       <- -64.157488
lat_large_anticosti_50km  <- 48.712857;  lon_large_anticosti_50km  <- -61.341346
lat_large_anticosti_25km  <- 48.877082;  lon_large_anticosti_25km  <- -61.584725


# ── Step 1: Build the migration line ─────────────────────────────────────────
mouth_sj               <- st_transform(st_sfc(st_point(c(lon_sj, lat_sj)),                       crs = wgs84), nad83_can_lambert)
via_mouth_gaspe_bay    <- st_transform(st_sfc(st_point(c(lon_mouth_gaspe_bay, lat_mouth_gaspe_bay)), crs = wgs84), nad83_can_lambert)
via_large_anticosti_50km <- st_transform(st_sfc(st_point(c(lon_large_anticosti_50km, lat_large_anticosti_50km)), crs = wgs84), nad83_can_lambert)
via_large_anticosti_25km <- st_transform(st_sfc(st_point(c(lon_large_anticosti_25km, lat_large_anticosti_25km)), crs = wgs84), nad83_can_lambert)

# Midpoint of the Belle Isle limit line (corridor endpoint)
end1_tn  <- st_point(st_coordinates(limit_lab_tn_line)[1, 1:2])
end2_lab <- st_point(st_coordinates(limit_lab_tn_line)[2, 1:2])
mid_end  <- st_sfc(st_point(c(mean(c(end1_tn[1], end2_lab[1])),
                              mean(c(end1_tn[2], end2_lab[2])))),
                   crs = nad83_can_lambert)

migration_line_sj_50km <- st_sf(
  geometry = st_sfc(st_linestring(rbind(
    st_point(st_coordinates(via_mouth_gaspe_bay)),
    st_point(st_coordinates(via_large_anticosti_50km)),
    st_point(st_coordinates(mid_end))
  )), crs = nad83_can_lambert)
)

migration_line_sj_25km <- st_sf(
  geometry = st_sfc(st_linestring(rbind(
    st_point(st_coordinates(via_mouth_gaspe_bay)),
    st_point(st_coordinates(via_large_anticosti_25km)),
    st_point(st_coordinates(mid_end))
  )), crs = nad83_can_lambert)
)


# ── Step 2: Buffer the migration lines ───────────────────────────────────────
migration_corridor_sj_50km <- st_buffer(migration_line_sj_50km, migration_route_width / 2, endCapStyle = "SQUARE")
migration_corridor_sj_25km <- st_buffer(migration_line_sj_25km, migration_route_width / 4, endCapStyle = "SQUARE")


# ── Step 3: Remove land overlap and clip at Belle Isle line ──────────────────
sj_area_final_50km  <- st_difference(migration_corridor_sj_50km, st_union(world_parts[c(10, 13), ]))
area_sj_final_25km  <- st_difference(migration_corridor_sj_25km, st_union(world_parts[c(10, 12, 13, 14), ]))

sj_area_split      <- st_split(sj_area_final_50km, limit_lab_tn_line)
sj_area_final_50km <- st_collection_extract(sj_area_split, "POLYGON")[1, ]


# ── Diagnostic plot ───────────────────────────────────────────────────────────
ggplot() +
  geom_sf(data = world_nad83) +
  geom_sf(data = sj_area_final_50km, fill = "blue", color = "blue", size = 1, alpha = 0.4) +
  coord_sf(xlim = c(x_west_ns_nad83, x_east_ns_nad83),
           ylim = c(y_south_ns_nad83, y_north_ns_nad83), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_blank())


# =============================================================================
# Strait of Belle Isle
# =============================================================================
# Strategy: define a quadrilateral polygon spanning the Strait between the
# southern limit line (shared with sections 1-2) and northern open-water
# anchor points on the Labrador and Newfoundland coasts.

# ── Key coordinates ───────────────────────────────────────────────────────────
wgs84_north_lab_point_bis <- c(-55.603023, 52.195050)
wgs84_north_tn_point_bis  <- c(-55.426483, 51.399200)


# ── Step 1: Snap corner points to coastlines ─────────────────────────────────
north_lab_point_bis <- st_transform(st_sfc(st_point(wgs84_north_lab_point_bis), crs = wgs84), nad83_can_lambert)
north_lab_point_bis <- st_cast(st_nearest_points(north_lab_point_bis, world_parts[10, ]), "POINT")[2]

north_tn_point_bis  <- st_transform(st_sfc(st_point(wgs84_north_tn_point_bis), crs = wgs84), nad83_can_lambert)
north_tn_point_bis  <- st_cast(st_nearest_points(north_tn_point_bis, world_parts[13, ]), "POINT")[2]

south_lab_point_bis <- end1_tn   # Southern corners = endpoints of limit_lab_tn_line
south_tn_point_bis  <- end2_lab


# ── Step 2: Build and clip the polygon ───────────────────────────────────────
bis_square <- matrix(c(
  st_coordinates(south_lab_point_bis)[1], st_coordinates(south_lab_point_bis)[2],
  st_coordinates(south_tn_point_bis)[1],  st_coordinates(south_tn_point_bis)[2],
  st_coordinates(north_tn_point_bis)[1],  st_coordinates(north_tn_point_bis)[2],
  st_coordinates(north_lab_point_bis)[1], st_coordinates(north_lab_point_bis)[2],
  st_coordinates(south_lab_point_bis)[1], st_coordinates(south_lab_point_bis)[2]
), ncol = 2, byrow = TRUE)
bis_square <- st_sfc(st_polygon(list(bis_square)), crs = nad83_can_lambert)

bis_area <- st_difference(bis_square, world_nad83)[2]
bis_area <- st_cast(bis_area, "POLYGON")[1]


# ── Diagnostic plot ───────────────────────────────────────────────────────────
ggplot(data = world_nad83) +
  geom_sf() +
  geom_sf(data = bis_area, fill = "blue", color = "blue", size = 1, alpha = 0.4) +
  coord_sf(xlim = c(x_west_ns_nad83, x_east_ns_nad83),
           ylim = c(y_south_ns_nad83, y_north_ns_nad83 + 500000), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_blank())


# =============================================================================
# Labrador Sea entrance (post-Strait zone)
# =============================================================================
# Strategy: construct a ~300 × 300 km square starting just north of the Strait
# of Belle Isle, representing the first open-ocean zone encountered by salmon
# after exiting the Strait.

# ── World parts: extended to full domain (Greenland + Quebec + NL) ────────────
world_parts_all      <- st_cast(world_original_nad83, "POLYGON")
world_parts_all$ID   <- 1:nrow(world_parts_all)
# Key IDs in the full-extent map:
#   128 = Greenland
#   148 = Quebec / Labrador
#   169 = Newfoundland

world_parts_all_lab_tn <- world_parts_all[c(128, 148, 169), ]

# Clip to east of longitude −70° to discard irrelevant western landmasses
limit_west   <- st_linestring(rbind(c(-70, 30), c(-70, 90))) |> st_sfc(crs = wgs84)
limit_west   <- st_transform(limit_west, st_crs(world_parts_all_lab_tn))
limit_west_x <- st_coordinates(limit_west)[1, "X"]
bbox         <- st_bbox(world_parts_all_lab_tn)
polygone_est <- st_sfc(st_polygon(list(rbind(
  c(limit_west_x, bbox["ymin"]),
  c(bbox["xmax"], bbox["ymin"]),
  c(bbox["xmax"], bbox["ymax"]),
  c(limit_west_x, bbox["ymax"]),
  c(limit_west_x, bbox["ymin"])
))), crs = st_crs(world_parts_all_lab_tn))
world_parts_all_lab_tn <- st_intersection(world_parts_all_lab_tn, polygone_est)


# ── Step 1: Build the square just north of the Strait ─────────────────────────
x_start_square_bis <- -55.774177 + 0.5
y_start_square_bis <-  52.779143

coords_square_bis_wgs84 <- st_sfc(st_point(c(x_start_square_bis, y_start_square_bis)), crs = 4326)
coords_square_bis_nad83 <- st_transform(coords_square_bis_wgs84, crs = nad83_can_lambert)
coords_square_bis_nad83 <- st_coordinates(coords_square_bis_nad83)
coords_square_bis_nad83[2] <- coords_square_bis_nad83[2] + 50000  # Shift 50 km northward

dimemsions_square_bis <- 300000
move                  <- 50000

coords_square_bis_nad83 <- matrix(c(
  coords_square_bis_nad83[1] - move,                         coords_square_bis_nad83[2],
  coords_square_bis_nad83[1] + dimemsions_square_bis,        coords_square_bis_nad83[2],
  coords_square_bis_nad83[1] + dimemsions_square_bis,        coords_square_bis_nad83[2] - dimemsions_square_bis,
  coords_square_bis_nad83[1],                                coords_square_bis_nad83[2] - dimemsions_square_bis,
  coords_square_bis_nad83[1] - move,                         coords_square_bis_nad83[2]
), ncol = 2, byrow = TRUE)
square_bis_nad83 <- st_sfc(st_polygon(list(coords_square_bis_nad83)), crs = nad83_can_lambert)


# ── Step 2: Remove overlapping land and Strait polygon ───────────────────────
square_bis_nad83 <- st_difference(square_bis_nad83, bis_area)
square_bis_nad83 <- st_difference(square_bis_nad83, st_union(world_parts_all_lab_tn))


# ── Diagnostic plot ───────────────────────────────────────────────────────────
ggplot() +
  geom_sf(data = world_parts_all_lab_tn) +
  geom_sf(data = square_bis_nad83, fill = "blue", color = "blue", size = 1, alpha = 0.4) +
  coord_sf(
    xlim = c(
      st_coordinates(st_transform(st_sfc(st_point(c(-60, 55)), crs = wgs84), nad83_can_lambert))[1],
      st_coordinates(st_transform(st_sfc(st_point(c(-50, 50)), crs = wgs84), nad83_can_lambert))[1]
    ),
    ylim = c(
      st_coordinates(st_transform(st_sfc(st_point(c(-55, 48)), crs = wgs84), nad83_can_lambert))[2],
      st_coordinates(st_transform(st_sfc(st_point(c(-60, 55)), crs = wgs84), nad83_can_lambert))[2]
    ),
    expand = FALSE
  ) +
  theme(panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_blank())


# =============================================================================
# Labrador Sea (offshore zone)
# =============================================================================
# Strategy: buffer NAFO sampling stations from Bradbury et al. (2021) to
# define the offshore Labrador Sea zone where Atlantic salmon overwinter.
# The Labrador coast buffer is constructed separately and the Strait polygon
# is excluded to avoid spatial overlap.
#
# Reference: Bradbury et al. (2021) ICES J. Mar. Sci., 78(4), 1434-1443.
# (Range-wide genetic assignment confirms long-distance oceanic migration)

# ── Step 1: Labrador coast buffer (boundary polygon) ─────────────────────────
# Custom polygon covering the relevant portion of the Labrador coast and
# eastern Newfoundland coast — used to extract the marine area of interest.
coords_lab_select <- matrix(c(
  -64,       61,
  -64,       59.822586,
  -64,       55.132846,
  -58.120901, 54.322472,
  -57.502102, 53.210424,
  wgs84_north_lab_point_bis, 
  -55.906735, 51.593859,
  -57.608837, 49.403272,
  -54,       47.906144,
  -54,       45,
  -40,       45,
  -40,       61,
  -64,       61
), ncol = 2, byrow = TRUE)

poly_lab_select <- st_sfc(st_polygon(list(coords_lab_select)), crs = wgs84)
poly_lab_select <- st_transform(poly_lab_select, crs = nad83_can_lambert)

lab_buffer_area <- st_difference(poly_lab_select, st_union(world_parts_all_lab_tn))
lab_buffer_area <- st_difference(lab_buffer_area, bis_area)


# ── Step 2: Offshore station buffers (Bradbury et al. 2021) ──────────────────
lab_sampling_coord <- data.frame(
  Station   = c("LBS_1F", "LBS_2G", "LBS_2H", "LBS_2J"),
  Latitude  = c(58.436, 56.511, 54.089, 50.804),
  Longitude = c(-55.58, -52.37, -50.67, -48.8)
)

# ── Diagnostic plot ───────────────────────────────────────────────────────────
ggplot() +
  geom_sf(data = world_original_nad83) +
  geom_sf(data = lab_buffer_area, fill = "blue", color = "blue", size = 1, alpha = 0.4) +
  coord_sf(xlim = c(x_west_ns_nad83, x_east_ns_nad83 + 1000000),
           ylim = c(y_south_ns_nad83 - 300000, y_north_ns_nad83 + 1100000),
           expand = FALSE) +
  theme_minimal()


# =============================================================================
# Export shapefiles
# =============================================================================
# All five polygons are dissolved to single-part geometries before export.
# Output directory: 01_RawData/ReferenceMigrationLayers/

tr_area_ns_final_unique  <- st_union(tr_area_ns_final)
sj_area_final_50km_unique <- st_union(sj_area_final_50km)

st_write(tr_area_ns_final_unique,   "R/00_raw_data/ReferenceMigrationLayers/TriniteMigrationCorridor.shp", delete_layer = TRUE)
st_write(sj_area_final_50km_unique, "R/00_raw_data/ReferenceMigrationLayers/StJeanMigrationCorridor.shp",  delete_layer = TRUE)
st_write(bis_area,                  "R/00_raw_data/ReferenceMigrationLayers/BelleIsleStrait.shp",          delete_layer = TRUE)
st_write(square_bis_nad83,          "R/00_raw_data/ReferenceMigrationLayers/LabradorSeaEntrance.shp",      delete_layer = TRUE)
st_write(lab_buffer_area,           "R/00_raw_data/ReferenceMigrationLayers/LabradorSea.shp",              delete_layer = TRUE)




