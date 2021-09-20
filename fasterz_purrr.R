# ==============================================
# Outside of loop
# ==============================================
# define continent
continent_vector <- as.character(c("Antarctica"))

# load world map polygons
world <- st_read(dsn = "0_raw data/ne_10m_admin_0_countries/ne_10m_admin_0_countries",
                 layer = "ne_10m_admin_0_countries",
                 crs = 4326)
st_crs(world)
crs(world)

# adjust map
world <- world[world$CONTINENT %!in% continent_vector,]

# save the map's projection
crs_longlat <- st_crs(world)
crs_longlat <- crs(world, asText = TRUE)

# mercator projection
crs_merc <- "+proj=merc +ellps=GRS80 +units=m +no_defs"

crs_4326 <- CRS("+init=epsg:4326")
new <- st_crs(4326, asText = T)

# ==============================================
# Load cities
# ==============================================
cities_master <- load_cities_func()

cities_master <- cities_master %>% 
  mutate_all(na_if,"") %>% 
  filter(!is.na(latitude) & !is.na(longitude) & !is.na(city)) %>%
  mutate(city = factor(city),
         country = factor(country)) %>%
  arrange(country, city, year)

cities_master <- cities_master[cities_master$CONTINENT %!in% continent_vector,]

r <- raster(ncol = 360, nrow = 180)
extent(r) <- extent(world)
crs(r) <- crs_longlat

sf_world <- rio::import("2_data in stages/11_non_spatial/sf_world.Rdata")

# ==============================================
# Remove contemporary borders borders
# ==============================================
world <- st_union(world)
world <- st_transform(world, crs = crs_merc)

ww <- as(world, "Spatial")
worldOwin <- as.owin(ww)
class(worldOwin)

# ===========================================
# SPLIT UP DATA
# ===========================================
year_sec <- seq(0, 1500, 10)
#year_sec <- c(400, 1000, 1600)
length(year_sec)

cities_purrr <- cities_master %>% 
  filter(year %in% year_sec)

table(cities_purrr$year)

cit_list <- split(cities_purrr,
                  f = cities_purrr$year)

# ===========================================
# (1) DEFINE VORONOI FUNCTION
# ===========================================
fun_v1 <- function(cities_raw, cities_spatial, cities_voronoi, voronoi_sf) {
  
  cities_spatial <- SpatialPointsDataFrame(coords = data.frame(cities_raw$longitude,
                                                               cities_raw$latitude),
                                           data = cities_raw, proj4string = CRS(crs_longlat))
  st_crs(cities_spatial)
  cities_spatial <- spTransform(cities_spatial, crs_merc)
  
  cities_voronoi <- as.ppp.SpatialPoints(cities_spatial)
  
  Window(cities_voronoi) <- worldOwin
  
  cities_voronoi <- as(dirichlet(cities_voronoi), "SpatialPolygons")
  
  proj4string(cities_voronoi) <- crs_merc
  cities_voronoi <- spTransform(cities_voronoi, crs_longlat)
  
  cities_voronoi <- as(cities_voronoi, "SpatialPolygonsDataFrame")
  
  cities_spatial <- spTransform(cities_spatial, crs_longlat)
  cities_voronoi@data <- over(cities_voronoi, cities_spatial)
  
  voronoi_sf <- st_as_sf(cities_voronoi)
  voronoi_sf$polygon_area_km2 <- unclass(st_area(voronoi_sf))/1000000
  
  return(voronoi_sf)
  
}

# not in function (): voronoi_crs <- crs(voronoi_sf, asText = TRUE)

# ===========================================
# (2) DEFINE RASTERIZATION FUNCTION
# ===========================================
fun_v2 <- function(voronoi_sf, voronoi_df, rp_ratified, rat, sf_dat, polygon_area) {
  
  polygon_area <- voronoi_sf %>%
    as.data.frame() %>%
    dplyr::select(city_yr, polygon_area_km2)
  
  voronoi_df <- data.frame(ID = 1:length(unique(voronoi_sf$city_yr)),
                           voronoi = unique(voronoi_sf$city_yr))
  voronoi_sf$ID <- voronoi_df$ID[match(voronoi_sf$city_yr, voronoi_df$voronoi)]
  
  #extent(r) <- extent(voronoi_sf)
  
  rp_ratified <- fasterize(voronoi_sf,
                        r,
                        field = "ID",
                        fun = "last",
                        background = NA_real_,
                        by = NULL)
  
  rp_ratified <- ratify(rp_ratified)

  rat <- levels(rp_ratified)[[1]]
  rat$voronoi <- voronoi_df$voronoi[match(rat$ID, voronoi_df$ID)]
  rat$IDs <- voronoi_df$ID[match(rat$voronoi, voronoi_df$voronoi)]
  levels(rp_ratified) <- rat
  
  sf_dat <- as(rp_ratified, Class = "SpatialPolygonsDataFrame") %>% 
    st_as_sf() %>% 
    mutate(ID = layer) %>%
    left_join(., rat, by = "ID") %>%
    mutate(city_yr = as.character(voronoi)) %>%
    dplyr::select(c(ID, city_yr)) %>% # geometry?
    left_join(., cities_master, by = "city_yr") %>% ### was: cities_raw
    mutate(cell_id = row_number(),
           cell_year_id = paste(cell_id,year, sep = "_")) %>% 
    left_join(., sf_world, by = "cell_id") %>%
    left_join(., polygon_area, by = "city_yr") %>%
    as.data.frame() %>%
    dplyr::select(cell_id,
                  year,
                  GEOUNIT_cell,
                  city,
                  city_id,
                  GEOUNIT,
                  pol_control_i,
                  city_pop_i,
                  polygon_area_km2,
                  latitude,
                  longitude,
                  source,
                  iso3,
                  CONTINENT_cell,
                  REGION_WB_cell,
                  REGION_UN_cell,
                  SUBREGION_cell,
                  iso3_cell,
                  lon_cell,
                  lat_cell,
                  cell_area_km2,
                  cell_year_id,
                  city_pop,
                  intpl_length,
                  CONTINENT,
                  REGION_WB,
                  REGION_UN,
                  SUBREGION,
                  geometry)
  
  return(sf_dat)
  
}

# ===========================================
# RUN FUNCTION (1)
# ===========================================
library(progressr)

future::plan(multiprocess)
tic()
with_progress({
  
  p <- progressor(steps = length(cit_list))
  
  step_1x <- future_map(cit_list, fun_v1, .progress = T)
  
})
toc()
future::plan(sequential)

future::plan(multiprocess)  # set workers = ... + gc = TRUE?
tic()
step_1 <- future_map(cit_list, fun_v1, .progress = T)
toc()
future::plan(sequential)

panel_polygons <- bind_rows(step_1, .id = "year_id") %>%
  dplyr::select(c(city_id, year, city_yr, geometry))

st_write(panel_polygons, dsn = "2_data in stages/10_polygons_df/10_polygons.shp",
         delete_dsn = T)

# ===========================================
# RUN FUNCTION (2)
# ===========================================
future::plan(multiprocess)
tic()
step_2 <- future_map(step_1, fun_v2, .progress = T)
toc()
future::plan(sequential)

# ===========================================
# UNLIST DATA
# ===========================================
panel_grid <- bind_rows(step_2, .id = "year_id") %>% 
  dplyr::select(-c(year_id, geometry)) %>% 
  arrange(year, cell_id) %>% 
  mutate(cell_city = paste(cell_id, city_id, sep = "_"))

nrow(panel_grid)/length(year_sec)
nrow(panel_grid)/length(year_sec) - nrow(sf_world)

# ===========================================
# CALCULATE CELL-TO-CITY DISTANCES
# ===========================================
# === Save computing power by calculating distance for unique 'cell-city' pairs ===
uniques <- panel_grid %>% 
  distinct(., cell_city, .keep_all = TRUE)

# === City geometries ===
city_coords <- uniques %>% 
  as.data.frame() %>% 
  dplyr::select(cell_city,latitude,longitude) %>% 
  st_as_sf(., coords = c("longitude", "latitude"),
           crs = 4326)

city_cell_matching <- city_coords$cell_city

city_coords <- city_coords %>% 
  st_geometry()

# === Cell geometries ===
cell_coords <- uniques %>% 
  as.data.frame() %>% 
  dplyr::select(cell_city, lat_cell, lon_cell) %>% 
  st_as_sf(., coords = c("lon_cell", "lat_cell"),
           crs = 4326) %>%   
  st_geometry()

# === Distance vector ===
dist_vector <- st_distance(cell_coords,
                           city_coords,
                           by_element = TRUE)

dist <- dist_vector %>%
  unclass() %>% 
  as.data.frame() %>% 
  mutate(cell_city = city_cell_matching)

data.table::setnames(dist,
                     old = c('.'),
                     new = c('dist_m'),
                     skip_absent = TRUE)

any(panel_grid$cell_city %!in% dist$cell_city)
any(dist$cell_city %!in% panel_grid$cell_city)

panel_grid <- panel_grid %>% 
  left_join(., dist, by = "cell_city") %>% 
  dplyr::select(-c(cell_city))

any(is.na(panel_grid$city))
any(is.na(panel_grid$year))
any(is.na(panel_grid$dist_m))

# ===========================================
# EXPORT DATA
# ===========================================
save(panel_grid, file = "2_data in stages/10_grid_panel.Rdata")

haven::write_dta(data = panel_grid,
                 path = "2_data in stages/10_grid_panel.dta",
                 version = 14)


