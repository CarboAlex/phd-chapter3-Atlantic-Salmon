

# |========================|
# |     LOAD LIBRARIES     |
# |========================|

library(dplyr)
library(lubridate)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(terra)
library(raster)
library(gdistance)
library(reshape2)
library(purrr)
library(fst)



# ══════════════════════════════════════════════════════════════════════════════
# DATA IMPORT
# ══════════════════════════════════════════════════════════════════════════════

# ── Telemetry data ────────────────────────────────────────────────────────────
# SOURCE: Lilly et al. (unpublished manuscript, under review).
# STATUS: Confidential — not included in this repository.
#         This file is a subset of the authors' full telemetry dataset.
#         An example file with identical structure and simulated data is
#         available at the same path with the suffix _EXAMPLE.csv.
#         To request access to the full dataset, contact the authors directly
#         once the manuscript has been published.
#
# To run with example data, replace the filename below with:
#   "master_dataset_smolt_manuscript_EXAMPLE.csv"


base_path <- "C:/Users/carbo/Documents/Alex/Ecole/Doctorat/Chapitre_3"
tel <- read.csv(
  file  = file.path(base_path, "01_RawData/Telemetry_Lilly_et_al/master_dataset_smolt_manuscript_dec_10_2025 QC.csv"),
  sep    = ";",
  dec    = ",",
  header = TRUE
)

tel <- read.csv("C:/Users/carbo/Downloads/master_dataset_smolt_manuscript_EXAMPLE.csv",
                sep    = ";",
                dec    = ",",
                header = TRUE)



# |========================|
# |     PLOT RECEPTORS     |
# |========================|

# Load coastlines
coastline <- ne_coastline(scale = "large", returnclass = "sf")

# Create a sf object with needed stations
stations <- tel %>%
  dplyr::distinct(Station_Name, Latitude, Longitude, Location) %>%
  dplyr::filter(Location %in% c("gulf_of_st_lawrence_1", "strait_of_belle_isle", "st_anthony", "southern_labrador", "newfoundland_1"))
stations_sf <- st_as_sf(
  stations,
  coords = c("Longitude", "Latitude"),
  crs = 4326   # WGS84
)

# Map limits
bbox <- st_bbox(c(
  xmin = -59, xmax = -54,
  ymin = 50,  ymax = 53
), crs = st_crs(4326))

# Plot map
map <- ggplot() +
  geom_sf(data = coastline, color = "grey40", linewidth = 0.3) +
  geom_sf(data = stations_sf, size = 2, aes(color = Location)) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"]),
           expand = FALSE) +
  scale_color_discrete(
    labels = c(
      gulf_of_st_lawrence_1 = "Gulf of St. Lawrence",
      strait_of_belle_isle = "Strait of Belle Isle (SOBI)",
      st_anthony = "St-Anthony (STAN)",
      southern_labrador = "Southern Labrador (SLAB)",
      newfoundland_1 = "Newfoundland 1"
    )
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 13)
map


# |=================================================|
# |     PREPARE DATA TO EXTRACT MIGRATION DATES     |
# |=================================================|

tel <- tel %>%
  # Rename column names
  dplyr::rename(id=serialno_tag, date_hit=Date_Hit, hit_location=Location, origin_river=river) %>%
  # Keep one observation per individual per day
  dplyr::group_by(id, date_hit) %>%
  dplyr::slice(1) %>% 
  dplyr::ungroup() %>%
  # Make sure date is in Date
  dplyr::mutate(date_hit = as.Date(date_hit))



# |======================================================================================================|
# |     ESTABLISH ARRIVAL DATE AT SOBI FOR TIME SPENT IN ST. LAWRENCE GULF AND DATE ENTERING IN SOBI     |
# |======================================================================================================|
# This date will be used to determine the time spent in the Gulf and for entering date in SOBI, so only the first date of detection at SOBI per individual is considered.
# The distribution created here will be used to determine the entry date into SOBI: 25th percentile of the distribution.

# ---------- KEEP THE FIRST OBSERVATION BY id × deploy_year ---------- 
# St. Jean
tel_first_obs_stj <- tel %>%
  # Keep only locations of interest and rivers of interest
  dplyr::filter(Latitude >= 50.5,
                hit_location %in% c("strait_of_belle_isle", "gulf_of_st_lawrence_1"),
                origin_river %in% c("Restigouche River", "Cascapedia River", "Saint Jean River")) %>%
  dplyr::mutate(date_hit = as.Date(date_hit)) %>%
  
  # Initial detection by id × origin_river × hit_location × deploy_year
  dplyr::group_by(id, deploy_year, origin_river, hit_location) %>%
  dplyr::slice_min(date_hit, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  
  # Apply priority for strait_of_belle_isle
  dplyr::group_by(id, deploy_year) %>%
  # If strait_of_belle_isle exists for this id × year, keep only this line
  dplyr::filter(!(any(hit_location == "strait_of_belle_isle") & hit_location == "gulf_of_st_lawrence_1")) %>%
  dplyr::ungroup() %>%
  
  # Keep only the desired columns
  dplyr::select(id, date_hit, origin_river, hit_location, deploy_year)
nrow(tel_first_obs_stj)

# Number of individual per origin river and per hit location
tel_first_obs_stj %>% dplyr::group_by(origin_river, hit_location) %>% dplyr::summarise(n = length(hit_location))

sobi_daily_detections_stj <- tel_first_obs_stj %>%
  dplyr::mutate(
    year = deploy_year,
    doy  = yday(date_hit)
  ) %>%
  count(year, doy, name = "n_detections")

# Trinite
tel_first_obs_tri <- tel %>%
  # Keep only locations of interest and rivers of interest
  dplyr::filter(Latitude >= 50.5,
                hit_location %in% c("strait_of_belle_isle", "gulf_of_st_lawrence_1"),
                origin_river %in% c("Trinite River", "Saguenay River")) %>%
  dplyr::mutate(date_hit = as.Date(date_hit)) %>%
  
  # Initial detection by id × origin_river × hit_location × deploy_year
  dplyr::group_by(id, deploy_year, origin_river, hit_location) %>%
  dplyr::slice_min(date_hit, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  
  # Apply priority for strait_of_belle_isle
  dplyr::group_by(id, deploy_year) %>%
  # If strait_of_belle_isle exists for this id × year, keep only this line
  dplyr::filter(!(any(hit_location == "strait_of_belle_isle") & hit_location == "gulf_of_st_lawrence_1")) %>%
  dplyr::ungroup() %>%
  
  # Keep only the desired columns
  dplyr::select(id, date_hit, origin_river, hit_location, deploy_year)
nrow(tel_first_obs_tri)

# Number of individual per origin river and per hit location
tel_first_obs_tri %>% dplyr::group_by(origin_river, hit_location) %>% dplyr::summarise(n = length(hit_location))

sobi_daily_detections_tri <- tel_first_obs_tri %>%
  dplyr::mutate(
    year = deploy_year,
    doy  = yday(date_hit)
  ) %>%
  count(year, doy, name = "n_detections")


# ---------- CALCULATE MEDIAN AND STANDARD DISTRIBUTION ----------
sobi_arrival_stj_median_sd <- data.frame(median_detection = median(sobi_daily_detections_stj$doy),
                                         sd_detection     = sd(sobi_daily_detections_stj$doy))
sobi_arrival_tri_median_sd <- data.frame(median_detection = median(sobi_daily_detections_tri$doy),
                                         sd_detection     = sd(sobi_daily_detections_tri$doy))


# |============================================================================================|
# |     ESTABLISH ARRIVAL DATE AT LABRADOR SEA FOR QUITTING SOBI AND ENTERING LABRADOR SEA     |
# |============================================================================================|
# This mean date will be used to determine the moment individuals quit SOBI, so only the last date of detection at Southern Labrador and at NFL per individual is considered.
# The distribution created here will be used to determine the SOBI quitting date for each river: 75th percentile of the distribution for St-Jean and 95th percentile for Trinite, because the sample sizes are small for the latter and the distributions overlap.

# ---------- KEEP THE LAST OBSERVATION BY id × deploy_year ----------
# St. Jean
tel_last_obs_stj <- tel %>%
  # Keep only locations of interest and rivers of interest
  dplyr::filter(Latitude >= 50.5,
                hit_location %in% c("southern_labrador", "newfoundland_1", "st_anthony"),
                origin_river %in% c("Restigouche River", "Cascapedia River", "Saint Jean River")) %>%
  dplyr::mutate(date_hit = as.Date(date_hit)) %>%
  
  # Last detection by id × origin_river × hit_location × deploy_year
  dplyr::group_by(id, deploy_year, origin_river, hit_location) %>%
  dplyr::slice_max(date_hit, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  
  # Priority order of sites
  dplyr::mutate(
    priority = dplyr::case_when(
      hit_location == "southern_labrador" ~ 1,
      hit_location == "newfoundland_1"    ~ 2,
      hit_location == "st_anthony"        ~ 3
    )
  ) %>%
  
  # Keep only ONE line per individual × year according to priority
  dplyr::group_by(id, deploy_year) %>%
  dplyr::slice_min(priority, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  
  dplyr::select(id, date_hit, origin_river, hit_location, deploy_year)
nrow(tel_last_obs_stj)

stan_slab_daily_detections_stj <- tel_last_obs_stj %>%
  dplyr::mutate(
    year = deploy_year,
    doy  = yday(date_hit)
  ) %>%
  count(year, doy, name = "n_detections")

# Trinite
tel_last_obs_tri <- tel %>%
  # Keep only locations of interest and rivers of interest
  dplyr::filter(Latitude >= 50.5,
                hit_location %in% c("southern_labrador", "newfoundland_1", "st_anthony"),
                origin_river %in% c("Trinite River", "Saguenay River")) %>%
  dplyr::mutate(date_hit = as.Date(date_hit)) %>%
  
  # Last detection by id × origin_river × hit_location × deploy_year
  dplyr::group_by(id, deploy_year, origin_river, hit_location) %>%
  dplyr::slice_max(date_hit, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  
  # Priority order of sites
  dplyr::mutate(
    priority = dplyr::case_when(
      hit_location == "southern_labrador" ~ 1,
      hit_location == "newfoundland_1"    ~ 2,
      hit_location == "st_anthony"        ~ 3
    )
  ) %>%
  
  # Keep only ONE line per individual × year according to priority
  dplyr::group_by(id, deploy_year) %>%
  dplyr::slice_min(priority, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  
  dplyr::select(id, date_hit, origin_river, hit_location, deploy_year)
nrow(tel_last_obs_tri)

stan_slab_daily_detections_tri <- tel_last_obs_tri %>%
  dplyr::mutate(
    year = deploy_year,
    doy  = yday(date_hit)
  ) %>%
  count(year, doy, name = "n_detections")



# ---------- CALCULATE MEAN AND STANDARD DISTRIBUTION ----------
stan_slab_daily_detections_both <- rbind(stan_slab_daily_detections_stj, stan_slab_daily_detections_tri)
stan_slab_arrival_stj_median_sd <- data.frame(median_detection = median(stan_slab_daily_detections_both$doy),
                                              sd_detection       = sd(stan_slab_daily_detections_both$doy))
stan_slab_arrival_tri_median_sd <- data.frame(median_detection = median(stan_slab_daily_detections_both$doy),
                                              sd_detection       = sd(stan_slab_daily_detections_both$doy))


# ---------- CREATE FINAL DATAFRAME OF DATES ----------
dates_final <- data.frame(river = c("stj", "stj", "tri", "tri"),
           arriving_to = c("sobi", "lab", "sobi", "lab"),
           date_doy = c(round(sobi_arrival_stj_median_sd$median_detection),
                        round(stan_slab_arrival_stj_median_sd$median_detection),
                        round(sobi_arrival_tri_median_sd$median_detection),
                        round(stan_slab_arrival_tri_median_sd$median_detection)),
           sd_date = c(round(sobi_arrival_stj_median_sd$sd_detection),
                       round(stan_slab_arrival_stj_median_sd$sd_detection),
                       round(sobi_arrival_tri_median_sd$sd_detection),
                       round(stan_slab_arrival_tri_median_sd$sd_detection))) %>%
  dplyr::mutate(date = format(as.Date(date_doy - 1, origin = "2001-01-01"), "%m-%d"))


# |===================================|
# |     EXPORT PLOT DISTRIBUTIONS     |
# |===================================|

pdf("02_PrepareData_Scripts/SST/MigrationDates_Distributions_FINAL.pdf", width = 12)

map

ggplot(dates_final, aes(x = date_doy, y = river, color = factor(arriving_to, levels = c("sobi", "lab")))) +
  geom_point(position = position_dodge(width = 0.4), size = 3) +
  geom_errorbarh(aes(xmin = date_doy - sd_date, xmax = date_doy + sd_date),
                 position = position_dodge(width = 0.4), height = 0.15, linewidth = 1) +
  geom_text(aes(label = format(as.Date(date_doy - 1, origin = "2001-01-01"), "%m-%d")),
            position = position_dodge(width = 0.4), vjust = -1.2, size = 3.5, show.legend = FALSE) +
  geom_text(aes(label = paste0("±", sd_date, " days")),
            position = position_dodge(width = 0.4), vjust = 1.8, size = 3, show.legend = FALSE) +
  scale_y_discrete(labels = c(stj = "St. Jean", tri = "Trinite")) +
  scale_x_continuous(
    name = "Date",
    limits = c(170, 240),
    breaks = seq(170, 240, by = 10),
    labels = function(x) format(as.Date(x - 1, origin = "2001-01-01"), "%m-%d")
  ) +
  labs(y = NULL, color = "Arriving at:") +
  scale_color_discrete(labels = c(sobi = "SOBI", lab  = "Labrador\nSea")) +
  theme_minimal(base_size = 18)

dev.off()


# |=============================================================================|
# |     ESTABLISH POST-SMOLT SWIM SPEED FROM LITTERATURE AND TELEMETRY DATA     |
# |=============================================================================|

# ---------- LITTERATURE SWIMING DISTANCES AND SPEED ----------

# | Source               | Distance / Time              | Average speed            |
# | -------------------- | ---------------------------- | ------------------------ |
# | Lefevre et al. 2012  | ~630 - 800 km in 44 d        | ~14.4 - 18.2 km/d        |
# | Chaput et al. 2019   | ~800 km in ~36-48 d          | ~17 – 22 km/d            |
# | Shelton et al. 1997  | 713 – 874 km in 38 – 51 d    | 14–21 km/d               |
# | Holst et al. 1996    | ~2000 km                     | ~22 km/d                 |
# | Reddin & Short, 1991 |                              | ~22 km/d                 |

# Lefevre, M. A., Stokesbury, M. J., Whoriskey, F. G., & Dadswell, M. J. (2012). Atlantic salmon post-smolt migration routes in the Gulf of St. Lawrence. ICES Journal of Marine Science, 69(6), 981-990.
# Chaput, G., Carr, J., Daniels, J., Tinker, S., Jonsen, I., & Whoriskey, F. (2019). Atlantic salmon (Salmo salar) smolt and early post-smolt migration and survival inferred from multi-year and multi-stock acoustic telemetry studies in the Gulf of St. Lawrence, northwest Atlantic. ICES Journal of Marine Science, 76(4), 1107-1121.
# Shelton, R. G. J., Turrell, W. R., Macdonald, A., McLaren, I. S., & Nicoll, N. T. (1997). Records of post-smolt Atlantic salmon, Salmo salar L., in the Faroe-Shetland Channel in June 1996. Fisheries Research, 31(1-2), 159-162.
# Holst, J. C., Hansen, L. P., & Holm, M. (1996). Observations of abundance, stock composition, body size and food of postsmolts of Atlantic salmon in the NE Atlantic during summer.
# Reddin, D. G., & Short, P. B. (1991). Postsmolt Atlantic salmon (Salmo salar) in the Labrador Sea. Canadian Journal of Fisheries and Aquatic Sciences, 48(1), 2-6.

# ---------- CREATE A DATASET OF INDIVIDUALS FOR WHICH THERE ARE DETECTIONS AT SOBI ENTRY AND EXIT ----------
# tel_entry_sobi <- tel %>%
#   dplyr::filter(hit_location %in% c("gulf_of_st_lawrence_1", "strait_of_belle_isle")) %>%
#   dplyr::rename(hit_location_entry = hit_location, date_hit_entry = date_hit, SerialNo_Rec_entry = SerialNo_Rec, Station_Name_entry = Station_Name, Latitude_entry = Latitude, Longitude_entry = Longitude) %>%
#   dplyr::select(id, origin_river, hit_location_entry, date_hit_entry, SerialNo_Rec_entry, Station_Name_entry, Latitude_entry, Longitude_entry, fl_mm)
# tel_exit_sobi <- tel %>%
#   dplyr::filter(hit_location %in% c("newfoundland_1", "southern_labrador", "st_anthony")) %>%
#   dplyr::rename(hit_location_exit = hit_location, date_hit_exit = date_hit, SerialNo_Rec_exit = SerialNo_Rec, Station_Name_exit = Station_Name, Latitude_exit = Latitude, Longitude_exit = Longitude) %>%
#   dplyr::select(id, origin_river, hit_location_exit, date_hit_exit, SerialNo_Rec_exit, Station_Name_exit, Latitude_exit, Longitude_exit, fl_mm)
# tel_distances_sobi <- merge(tel_entry_sobi, tel_exit_sobi, by = c("id", "origin_river"))


# ---------- CALCULATE LEAST COST DISTANCE BETWEEN RECEPTOR STATIONS ----------

# Extract stations
station_distances <- tel %>%
  dplyr::select(SerialNo_Rec, Station_Name, Latitude, Longitude, hit_location)

# Prepare the points (stations)
stations_sf <- station_distances %>%
  dplyr::filter(Latitude >= 50.5,
                hit_location %in% c("gulf_of_st_lawrence_1", "strait_of_belle_isle", "newfoundland_1", "southern_labrador", "st_anthony")) %>%
  distinct(SerialNo_Rec, Latitude, Longitude, hit_location) %>%   # remove duplicates
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

# Load land data
land <- ne_countries(scale = "large", returnclass = "sf")

ggplot() +
  geom_sf(data = land) +
  geom_sf(data = stations_sf, aes(color = hit_location)) +
  coord_sf(xlim = c(bbox["xmin"]-10, bbox["xmax"]+2),
           ylim = c(bbox["ymin"]-5, bbox["ymax"]),
           expand = FALSE)


# Create a cost raster grid (water vs. land)
ext <- ext(stations_sf) + 0.1 # Extend of area to consider
r <- rast(extent = ext, resolution = 0.001, crs = "EPSG:4326") # Create the raster with 0.01 degree resolution
values(r) <- 1 # Initialize: water = 1
r_land <- rasterize(vect(land), r, field = 1) # Rasterize the earth (infinite cost)
r[r_land == 1] <- NA # Terre = NA (impraticable)
plot(r)
plot(st_geometry(stations_sf), add = TRUE, pch = 16, col = "red", cex = 0.7)

# Complete distance matrix
stations_sf <- st_transform(stations_sf, 32620) # Convert in a meters projection: UTM zone 20N
land <- st_transform(land, 32620) # Convert in a meters projection: UTM zone 20N
r_raster <- raster(r) # Convert to "RasterLayer" since transition() needs this format
r_raster <- projectRaster(r_raster, crs = 32620, method = "bilinear") # Convert in a meters projection: UTM zone 20N
tr <- transition(r_raster, function(x) 1/mean(x), 8) # Build the transition surface (possible paths)
tr <- geoCorrection(tr, type="c") # Adjusts the transition values so that the distances calculated by costDistance better correspond to the actual distance on the field
dist_mat_m <- costDistance(tr, st_coordinates(stations_sf), st_coordinates(stations_sf)) # Calculate least cost distances (meters)


# Reshape so it's a long format
ids <- stations_sf$SerialNo_Rec
rownames(dist_mat_m) <- ids
colnames(dist_mat_m) <- ids
dist_long <- reshape2::melt(dist_mat_m, 
                            varnames = c("SerialNo_Rec_entry", "SerialNo_Rec_exit"), 
                            value.name = "distance_m")

# ---------- ADD DISTANCE TO INDIVIDUALS ----------
distances_id <- merge(tel_distances_sobi, dist_long, by = c("SerialNo_Rec_entry", "SerialNo_Rec_exit")) %>%
  dplyr::mutate(distance_km = distance_m / 1000) %>%
  dplyr::select(-distance_m)


# ---------- NUMBER OF DAYS BETWEEN 2 POINTS ----------
distances_id <- distances_id %>%
  # Keep the columns of interest and calculates travel time
  dplyr::select(id, origin_river, SerialNo_Rec_entry, SerialNo_Rec_exit, 
                hit_location_entry, hit_location_exit, date_hit_entry, 
                date_hit_exit, distance_km, fl_mm.y) %>%
  dplyr::mutate(
    date_hit_entry = yday(date_hit_entry),
    date_hit_exit  = yday(date_hit_exit),
    travel_time    = date_hit_exit - date_hit_entry
  ) %>%
  # Priority: strait_of_belle_isle > gulf_of_st_lawrence_1
  group_by(id, SerialNo_Rec_exit) %>%
  slice_min(order_by = if_else(hit_location_entry == "strait_of_belle_isle", 0, 1), n = 1) %>%
  ungroup() %>%
  filter(travel_time > 0) %>%
  dplyr::mutate(swim_speed = distance_km / travel_time,
                region = ifelse(origin_river %in% c("Restigouche River", "Cascapedia River", "Saint Jean River"), "Gaspesie", "North Shore"))

mean_swim_speed <- round(mean(distances_id$swim_speed)) # km / day
round(300 / mean_swim_speed)

# Published estimates of Atlantic salmon post-smolt swimming speeds generally fall within a range of approximately 14 to 22 km/d. Studies conducted 
# in the Gulf of St. Lawrence during the early marine phase, shortly after river exit, report mean speeds of ~14–18 km/d (Lefevre et al. 2012) and 
# ~17–22 km/d (Chaput et al. 2019). Similar values have been reported elsewhere in the North Atlantic, including the Faroe–Shetland Channel, where 
# post-smolts covered 713–874 km in 38–51 days (14–21 km/d; Shelton et al. 1997). Observations from more offshore or later migration stages indicate 
# comparable or slightly higher values, with mean speeds around 22 km/d reported in the Northeast Atlantic (Holst et al. 1996) and in the Labrador 
# Sea (Reddin & Short 1991).
# 
# It is important to note that many of these estimates are derived from fish in the first phase of their marine migration, often shortly after river 
# exit, when individuals are smaller and may exhibit lower sustained swimming speeds. In contrast, our estimate of a mean swimming speed of approximately 
# 23 km/d is based on fish originating from the same river systems and observed farther along their migration, in the Strait of Belle Isle region. At 
# this stage, individuals are expected to be slightly larger and more developmentally advanced, which is expected to translate into higher swimming 
# capacity. Consequently, although our estimate lies at the upper end of the range reported in the literature, it remains biologically plausible and 
# consistent with existing studies when differences in fish size and migration stage are taken into account.
# 
# Thus, considering that the polygon at the entrance to the Labrador Sea is 300 km wide, we estimate that it will take (300 km) / (23 km/day) = 13 days 
# to cross this polygon.


# ---------- PLOT MIGRATION ROUTES ----------

# Prepare stations with projected coordinates
stations_proj <- station_distances %>%
  distinct(SerialNo_Rec, Station_Name, Latitude, Longitude) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(32620)

# Single table of routes to be calculated
pairs <- distances_id %>%
  distinct(id, SerialNo_Rec_entry, SerialNo_Rec_exit, origin_river)

# Function that calculates a path
get_path <- function(entry, exit, fish_id, river){
  
  start_pt <- stations_proj %>% filter(SerialNo_Rec == entry)
  end_pt   <- stations_proj %>% filter(SerialNo_Rec == exit)
  
  if(nrow(start_pt) == 0 | nrow(end_pt) == 0) return(NULL)
  
  path_sp <- shortestPath(
    tr,
    st_coordinates(start_pt),
    st_coordinates(end_pt),
    output = "SpatialLines"
  )
  
  proj4string(path_sp) <- CRS("+proj=utm +zone=20 +datum=WGS84 +units=m +no_defs")
  
  path_sf <- st_as_sf(path_sp)
  
  path_sf$id <- fish_id
  path_sf$river <- river
  path_sf$entry <- entry
  path_sf$exit  <- exit
  
  return(path_sf)
}

# Calculate all paths
paths_list <- pmap(
  list(pairs$SerialNo_Rec_entry,
       pairs$SerialNo_Rec_exit,
       pairs$id,
       pairs$origin_river),
  get_path
)
paths_sf <- do.call(rbind, paths_list)

# Map
bbox_ll <- st_bbox(c(
  xmin = -59, xmax = -54,
  ymin = 50,  ymax = 53
), crs = st_crs(4326))

bbox_utm <- st_transform(st_as_sfc(bbox_ll), 32620) %>% st_bbox()

land_utm      <- st_transform(land, 32620)

ggplot() +
  geom_sf(data = land_utm, fill = "white", color = "grey50", linewidth = 0.4) +
  geom_sf(data = paths_sf, linewidth = 0.8, alpha = 0.6) +
  geom_sf(data = stations_proj, color = "RED", size = 1.5) +
  coord_sf(xlim = c(bbox_utm["xmin"], bbox_utm["xmax"]),
           ylim = c(bbox_utm["ymin"], bbox_utm["ymax"]),
           expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 13) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA))


# |======================================================================|
# |     ESTABLISH DATES FOR LABRADOR SEA EXIT (FOR BOTH 1SW AND 2SW)     |
# |======================================================================|

# According to Imlay et al. (2024):
### The summer period for post-smolts spanned August to November and during the 
### second year, the summer period spanned May to November; the winter period 
### spanned December to April in both years.

# Thus, May 1st is considered as the exit of Labrador Sea polygon for both 1SW and 2SW
lab_sea_1sw_exit <- yday(as.Date("2001-05-01"))
lab_sea_2sw_exit <- yday(as.Date("2001-05-01"))


# |=======================================================|
# |     CREATE FINAL DATASET WITH DATES FOR EACH YEAR     |
# |=======================================================|

all_dates_sea <- bind_rows(
  data.frame(
    river                  = "stj",
    cohort                 = 1980:2025,
    sobi_entry             = as.numeric(round(sobi_arrival_stj_median_sd$median_detection)),
    sobi_exit              = round(stan_slab_arrival_stj_median_sd$median_detection),
    entering_lab_sea_entry = as.numeric(round(stan_slab_arrival_stj_median_sd$median_detection)),
    entering_lab_sea_exit  = as.numeric(round(stan_slab_arrival_stj_median_sd$median_detection)) + round(300 / as.numeric(mean_swim_speed)),
    lab_sea_entry          = as.numeric(round(stan_slab_arrival_stj_median_sd$median_detection)),
    lab_sea_exit_1sw       = lab_sea_1sw_exit,
    lab_sea_exit_2sw       = lab_sea_2sw_exit    #,
    # gulf_entry_1sw         = NA, # Will have to be determined based on swimming speeds (number of days preceding the date of entry into the river).
    # gulf_exit_1sw          = NA, # Will have to be established based on the dates of entry into the river.
    # gulf_entry_2sw         = NA, # Will have to be determined based on swimming speeds (number of days preceding the date of entry into the river).
    # gulf_exit_2sw          = NA  # Will have to be established based on the dates of entry into the river.
  ),
  
  data.frame(
    river                  = "tri",
    cohort                 = 1980:2025,
    sobi_entry             = as.numeric(round(sobi_arrival_tri_median_sd$median_detection)),
    sobi_exit              = as.numeric(round(stan_slab_arrival_tri_median_sd$median_detection)),
    entering_lab_sea_entry = as.numeric(round(stan_slab_arrival_tri_median_sd$median_detection)),
    entering_lab_sea_exit  = as.numeric(round(stan_slab_arrival_tri_median_sd$median_detection)) + round(300 / as.numeric(mean_swim_speed)),
    lab_sea_entry          = as.numeric(round(stan_slab_arrival_tri_median_sd$median_detection)),
    lab_sea_exit_1sw       = lab_sea_1sw_exit,
    lab_sea_exit_2sw       = lab_sea_2sw_exit    #,
    # gulf_entry_1sw         = NA, # Will have to be determined based on swimming speeds (number of days preceding the date of entry into the river).
    # gulf_exit_1sw          = NA, # Will have to be established based on the dates of entry into the river.
    # gulf_entry_2sw         = NA, # Will have to be determined based on swimming speeds (number of days preceding the date of entry into the river).
    # gulf_exit_2sw          = NA  # Will have to be established based on the dates of entry into the river.
  )
)

# Import smolting dates
all_dates_smolt <- bind_rows(read_fst("02_PrepareData_Scripts/Intermediate_Cleaned_Data/Year_Median_Date_StJean.fst"),
                             read_fst("02_PrepareData_Scripts/Intermediate_Cleaned_Data/Year_Median_Date_Trinite.fst")) %>%
  dplyr::select(river, cohort, median_date_doy,quantile_25_date, quantile_75_date) %>%
  dplyr::rename(smolting_doy = median_date_doy, smolt_quant_25_doy = quantile_25_date, smolt_quant_75_doy = quantile_75_date) %>%
  dplyr::mutate(smolt_quant_25_doy = yday(smolt_quant_25_doy),
                smolt_quant_75_doy = yday(smolt_quant_75_doy))

# Combine sea migration dates and smolting dates 
all_dates_migration <- merge(all_dates_smolt, all_dates_sea, by = c("river", "cohort"))

# Export dates
write_fst(all_dates_migration, path = "02_PrepareData_Scripts/Intermediate_Cleaned_Data/Year_Migration_Dates_River_Sea.fst")





########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################







# gr_stj <- read_excel("~/Alex/Ecole/Doctorat/Chapitre_3/01_RawData/MELCCFP_Data/RawData_MELCCFP_AllYears/donnees-2001-2024/donnees-2001-2024/2014/Saumon/BD/ST-JEAN_SAISIE_2014.xlsx", 
#                      sheet = "GRACIATIONS") %>%
#   dplyr::select(DATE, CATEGORIE) %>%
#   dplyr::mutate(DATE = as.Date(DATE)) %>%
#   dplyr::group_by(DATE, CATEGORIE) %>%
#   dplyr::summarise(n = length(CATEGORIE)) %>%
#   as.data.frame() %>%
#   mutate(type = "gracie")
# capt_stj <- read_excel("~/Alex/Ecole/Doctorat/Chapitre_3/01_RawData/MELCCFP_Data/RawData_MELCCFP_AllYears/donnees-2001-2024/donnees-2001-2024/2014/Saumon/BD/ST-JEAN_SAISIE_2014.xlsx", 
#                      sheet = "Captures St-Jean", skip = 1) %>%
#   dplyr::select(DATE, CATEGORIE) %>%
#   dplyr::mutate(DATE = as.Date(DATE)) %>%
#   dplyr::group_by(DATE, CATEGORIE) %>%
#   dplyr::summarise(n = length(CATEGORIE)) %>%
#   as.data.frame() %>%
#   mutate(type = "captured")
# 
# test <- rbind(gr_stj, capt_stj)
# 
# 
# 
# date_limits <- range(capt_stj$DATE, na.rm = TRUE)
# 
# ggarrange(
# ggplot(capt_stj %>% filter(CATEGORIE == "M"), aes(x = DATE, y = n, fill = factor(type))) +
#   geom_col(alpha = 0.8) +
#   scale_x_date(limits = date_limits) +
#   theme_minimal(base_size = 14) +
#   labs(y = "Count", x = NULL, title = "1SW"),
# ggplot(capt_stj %>% filter(CATEGORIE == "R"), aes(x = DATE, y = n, fill = factor(type))) +
#   geom_col(alpha = 0.8) +
#   scale_x_date(limits = date_limits) +
#   theme_minimal(base_size = 14) +
#   labs(y = "Count", x = NULL, title = "2SW"),
# ncol = 1, common.legend = TRUE)



