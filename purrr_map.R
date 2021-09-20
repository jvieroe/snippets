# ==============================================
# Outside of loop
# ==============================================
# define continent
continent_vector <- as.character(c("Antarctica"))

# load world map polygons
world <- rgdal::readOGR(dsn = "0_raw data/ne_10m_admin_0_countries/ne_10m_admin_0_countries",
                        layer = "ne_10m_admin_0_countries")

# adjust map
world <- world[world@data$CONTINENT %!in% continent_vector,]

# save the map's projection
crs_longlat <- proj4string(world)

# mercator projection
crs_merc <- "+proj=merc +ellps=GRS80 +units=m +no_defs"

# ==============================================
# Remove contemporary borders borders
# ==============================================
world <- gUnaryUnion(world)
worldOwin <- spTransform(world, crs_merc)
worldOwin <- as.owin(worldOwin)

class(worldOwin)

# ==============================================
# Load cities
# ==============================================
cities_master <- load_cities_func()

cities_master <- cities_master %>% 
  mutate_all(na_if,"") %>% 
  filter(!is.na(latitude) & !is.na(longitude) & !is.na(city)) %>%
  mutate(city = factor(city),
         country = factor(country)) %>%
  mutate(city_id_year = paste(as.character(city_id), year, sep = "_")) %>% 
  arrange(country, city, year)

cities_master <- cities_master[cities_master$CONTINENT %!in% continent_vector,]

# === Remove geographical duplicates within years ===
dups <- cities_master %>%
  mutate(lat = round(latitude, digits = 3),
         lon = round(longitude, digits = 3)) %>%
  group_by(lat, lon, year) %>%
  mutate(dupe = n()>1)
dups <- dups %>%
  dplyr::select(city,country,latitude,longitude,dupe,lat,lon,pol_control,year,city_id) %>%
  filter(dupe == T) %>%
  arrange(year, latitude, longitude)

nn <- nrow(dups)/2

dups <- dups %>%
  mutate(toDelete = seq(1, nrow(dups), nn)) %>%
  filter(toDelete == 1)
cities_master <- cities_master[cities_master$city %!in% dups$city,]
rm(dups, nn)

r <- raster(ncol = 360, nrow = 180) # 720, 360
sf_world <- rio::import("2_data in stages/11_non_spatial/sf_world.Rdata")

# ===========================================
# SPLIT UP VERSION
# ===========================================
year_sec <- seq(100, 200, 10)
#year_sec <- c(400, 1000, 1600)
length(year_sec)

cities_purrr <- cities_master %>% 
  filter(year %in% year_sec)

table(cities_purrr$year)

cit_list <- split(cities_purrr,
                  f = cities_purrr$year)

# cit1 <- cit_list[[1]]
# cit_sp <- SpatialPointsDataFrame(coords = data.frame(cit1$longitude,
#                                                      cit1$latitude),
#                                  data = cit1, proj4string = CRS(crs_longlat))
# cit_sp <- spTransform(cit_sp, crs_merc)
# cit_vor <- as(cit_sp, "ppp")
# cit_vor2 <- as.ppp(cit_sp)
# class(cit_vor)
# class(cit_vor2)

fun_v1 <- function(cities_raw, cities_spatial, cities_voronoi, voronoi_sf) {
  
  cities_spatial <- SpatialPointsDataFrame(coords = data.frame(cities_raw$longitude,
                                                               cities_raw$latitude),
                                           data = cities_raw, proj4string = CRS(crs_longlat))
  
  cities_spatial <- spTransform(cities_spatial, crs_merc)
  
  #cities_voronoi <- as(cities_spatial, "ppp")
  
  #cities_voronoi <- as.ppp.ppp(cities_spatial, W = worldOwin)
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

# x <- as.ppp.ppp()
# as.ppp.SpatialPointsDataFrame(X)

fun_v2 <- function(voronoi_sf, voronoi_df, voronoi_crs, rp_ratified, rat, sf_dat, polygon_area) {
  
  polygon_area <- voronoi_sf %>%
    as.data.frame() %>%
    dplyr::select(city_id_year, polygon_area_km2)
  
  # Rasterize
  voronoi_df <- data.frame(ID = 1:length(unique(voronoi_sf$city_id_year)),
                           voronoi = unique(voronoi_sf$city_id_year))
  voronoi_sf$ID <- voronoi_df$ID[match(voronoi_sf$city_id_year, voronoi_df$voronoi)]
  
  voronoi_crs <- crs(voronoi_sf, asText = TRUE)
  
  extent(r) <- extent(voronoi_sf)
  
  rp_ratified <- ratify(rasterize(voronoi_sf, r, field = "ID"))
  
  #rp <- projectRaster(rp, crs = voronoi_crs, method = "ngb")
  
  rat <- levels(rp_ratified)[[1]]
  rat$voronoi <- voronoi_df$voronoi[match(rat$ID, voronoi_df$ID)]
  rat$IDs <- voronoi_df$ID[match(rat$voronoi, voronoi_df$voronoi)]
  levels(rp_ratified) <- rat
  
  sf_dat <- as(rp_ratified, Class = "SpatialPolygonsDataFrame")
  
  sf_dat <- st_as_sf(sf_dat)
  
  sf_dat <- sf_dat %>%
    mutate(ID = layer) %>%
    left_join(., rat, by = "ID") %>%
    mutate(city_id_year = as.character(voronoi)) %>%
    dplyr::select(c(ID, city_id_year)) %>%
    left_join(., cities_master, by = "city_id_year") %>% ### was: cities_raw
    mutate(cell_id = 1:nrow(sf_dat),
           cell_year_id = paste(cell_id,year, sep = "_")) %>%
    left_join(., sf_world, by = "cell_id") %>%
    left_join(., polygon_area, by = "city_id_year") %>%
    as.data.frame() %>%
    dplyr::select(cell_id,
                  year,
                  GEOUNIT_cell,
                  cell_area_km2,
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

# library(future)
# library(furrr)

future::plan(multiprocess)  # set workers = ... + gc = TRUE?
tic()
step_1 <- future_map(cit_list, fun_v1, .progress = T)
toc()
future::plan(sequential)

df1 <- step_1[[2]]
any(is.na(df1$city))
df2 <- df1 %>% 
  mutate(id = row_number()) %>% 
  filter(is.na(city)) %>% 
  dplyr::select(c(city, country, id))

df3 <- cities_purrr %>% 
  filter(year == 1940 & country == "russia") %>% 
  filter(city %!in% df1$city)

df4 <- cities_purrr %>% 
  filter(year == 1950 & country == "russia") %>%
  arrange(longitude, latitude)


rm(df1, df2)

(19249.52/60)/71

(471.48/60)/11

future::plan(multiprocess)
tic()
step_2 <- future_map(step_1, fun_v2, .progress = T)
toc()
future::plan(sequential)

df <- bind_rows(step_2, .id = "year_id")
class(df)

df <- df %>% 
  dplyr::select(-c(geometry))

nrow(df)/length(year_sec)

any(is.na(df$year))

table(df$year) # 140 / 150
table(cities_purrr$year)

# df2 <- df2 %>% 
#   mutate(year_id = as.numeric(year_id)) %>% 
#   mutate(year = year_id)


save(df, file = "2_data in stages/10_grid_panel.Rdata")

haven::write_dta(data = df,
                 path = "2_data in stages/10_grid_panel.dta",
                 version = 14)

# ==============================================
# Add cell-to-city distances
# ==============================================
rm(list=ls())

df <- rio::import("2_data in stages/10_grid_panel.Rdata")
table(df$year)

spatial_merge <- st_read(dsn = "2_data in stages/11_spatial_dimensions/11_spatial.shp",
                         layer = "11_spatial")

merge_df <- spatial_merge %>% 
  right_join(df, by = "cell_id") %>% 
  arrange(year, cell_id) %>% 
  mutate(cell_city = paste(cell_id, city_id, sep = "_"))
rm(df)

any(is.na(merge_df$cell_id))
any(duplicated(merge_df$cell_id))
tt <- duplicated(merge_df$cell_id) %>% 
  as.data.frame() %>% 
  mutate(cell_id = 1:nrow(tt)) %>% 
  filter(. == "TRUE")

table(merge_df$year)

head(merge_df)

dist_data <- merge_df %>% 
  distinct(., cell_city, .keep_all = TRUE) %>% 
  dplyr::select(c(cell_city,
                  latitude, longitude, geometry))

centroids <- st_centroid(dist_data) %>% 
  st_transform(crs = 4326) %>% 
  dplyr::select(cell_city,geometry)

city_coords <- dist_data %>% 
  as.data.frame() %>% 
  dplyr::select(latitude, longitude, cell_city)
city_coords <- st_as_sf(city_coords,
                        coords = c("longitude", "latitude"),
                        crs = 4326)

#city_coords <- st_geometry(city_coords)
#head(city_coords)

dist_vector <- st_distance(centroids, city_coords, by_element = TRUE)

dist <- dist_vector %>%
  unclass() %>% 
  as.data.frame()
dist <- dist %>% 
  mutate(match_id = 1:nrow(dist))
data.table::setnames(dist,
                     old = c('.'),
                     new = c('dist_m'))
rm(city_coords, dist_vector, centroids)

dist_data <- dist_data %>% 
  mutate(match_id = 1:nrow(dist)) %>% 
  left_join(., dist, by = "match_id") %>% 
  as.data.frame() %>% 
  dplyr::select(-c(latitude, longitude, match_id, geometry))

merge_df <- merge_df %>% 
  left_join(., dist_data, by = "cell_city")


avg_df <- merge_df %>% 
  as.data.frame() %>% 
  #dplyr::select(-geometry) %>% 
  group_by(cell_id) %>% 
  dplyr::summarize(avg_dist = mean(dist_m))
class(avg_df)

plot_dat <- merge_df %>% 
  filter(year == 1200) %>% 
  left_join(., avg_df, by = "cell_id")

# ggplot(data = plot_dat) +
#   geom_sf(aes(fill = avg_dist), color = NA)
  
