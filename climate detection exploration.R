library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(camtrapR)
library(ggfortify)
library(ggrepel)
library(forcats)
library(lubridate)
library(unmarked)
library(camtrapR)
library(overlap)
library(ggmap)
library(reshape2)
library(corrplot)
library(terra)
library(mapview)
library(sf)
library(exactextractr)
library(xml2)
library(geodata)
library(viridis)
library(scales)
library(colorspace)
library(spOccupancy)
library(ClimateNAr)
library(ubms)
library(data.table)
library(purrr)
library(ISOweek)
library(ggspatial)
library(maptiles)

#Read in camera detections
all_detections <-read.csv("~/chilcotin_bobcats/data/kwasi.bobcats.etc.csv")

#Read in environmental covariates
site_covars <- read.csv("~/chilcotin_bobcats/data/env.vars.2021.csv")

#Format date/time

all_detections$date.time <- parse_date_time(
  all_detections$date.time,
  orders = c("mdy HMS", "mdy HM", "mdy"),
  tz = "Canada/Pacific"
)


#Join LatLong and Elevation to detections

##Rename sites to station in covars

site_covars <- site_covars %>%
  rename(station = site)

all_detections <- all_detections %>%
  left_join(site_covars, by = "station")

#Visualize detections per year

species_year_counts <- all_detections %>%
  group_by(species, year) %>% 
  summarise(n_detections = n(), .groups = "drop")

species_year_counts <- species_year_counts %>%
  mutate(year = factor(year, levels = 2017:2025))%>%
  filter(!is.na(year))

ggplot(species_year_counts, aes(x = year, y = n_detections, fill = species)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ species, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Detections per Species per Year",
    x = "Year",
    y = "Number of Detections"
  )


#Overlap between stations

stations_sf <- site_covars %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Summaries per species × station
sp_station <- all_detections %>%
  group_by(species, station) %>%
  summarise(n = n(), .groups = "drop") %>%
  left_join(stations_sf, by = "station")

sp_pair <- sp_station %>% filter(species %in% c("bobcat", "lynx"))

ggplot(sp_pair) +
  geom_sf(aes(color = species, size = n), alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(title = "Spatial Overlap: Bobcat & Lynx",
       color = "Species", size = "Detections")

#Map detection of species per site

##Convert to sf objects
sites_sf <- site_covars %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)


det_summary <- all_detections %>%
  group_by(station, species) %>%
  summarise(n = n(), .groups = "drop")

sites_det <- sites_sf %>%
  left_join(det_summary, by = "station")

species_to_plot <- "bobcat"

ggplot() +
  geom_sf(data = sites_sf, color = "grey70") +
  geom_sf(data = filter(sites_det, species == species_to_plot),
          aes(size = n, color = n)) +
  scale_color_viridis_c() +
  labs(title = paste("Detections of", species_to_plot),
       size = "Number of detections",
       color = "Number of detections")


# get bounding box of your camera sites
bb <- st_bbox(sites_sf)

# download OSM tiles
osm_tiles <- get_tiles(bb, provider = "OpenStreetMap", zoom = 10)

ggplot() +
  layer_spatial(osm_tiles) +
  geom_sf(data = sites_sf, size = 2, color = "black") +
  geom_sf(data = filter(sites_det, species == "bobcat"),
          aes(color = n, size = n)) +
  scale_color_viridis_c() +
  coord_sf() +
  labs(title = "Bobcat detections with OSM basemap")

ggplot() +
  layer_spatial(osm_tiles) +
  geom_sf(data = sites_sf, size = 2, color = "black") +
  geom_sf(data = filter(sites_det, species == "canada lynx"),
          aes(color = n, size = n)) +
  scale_color_viridis_c() +
  coord_sf() +
  labs(title = "Lynx detections with OSM basemap")

ggplot() +
  layer_spatial(osm_tiles) +
  geom_sf(data = sites_sf, size = 2, color = "black") +
  geom_sf(data = filter(sites_det, species == "snowshoe hare"),
          aes(color = n, size = n)) +
  scale_color_viridis_c() +
  coord_sf() +
  labs(title = "Snowshoe hare detections with OSM basemap")

##Lynx across years

det_year <- all_detections %>%
  mutate(year = year(date.time))

lynx_yearly <- det_year %>%
  filter(species == "canada lynx") %>%
  group_by(station, year) %>%
  summarise(n = n(), .groups = "drop")

lynx_yearly_sf <- sites_sf %>%
  left_join(
    lynx_yearly %>% st_drop_geometry(),
    by = "station"
  ) %>%
  mutate(n = replace_na(n, 0))

ggplot() +
  layer_spatial(osm_tiles) +
  geom_sf(data = sites_sf, size = 1, color = "grey60") +
  geom_sf(data = filter(lynx_yearly_sf, year == 2020),
          aes(color = n, size = n)) +
  scale_color_viridis_c() +
  coord_sf() +
  labs(title = "Canada Lynx Detections — 2020")

ggplot(lynx_yearly_sf) +
  layer_spatial(osm_tiles) +
  geom_sf(aes(color = n, size = n)) +
  scale_color_viridis_c() +
  coord_sf() +
  facet_wrap(~ year) +
  labs(title = "Canada Lynx Detections by Year")

#Snowshoe hare across years

hare_yearly <- det_year %>%
  filter(species == "snowshoe hare") %>%
  group_by(station, year) %>%
  summarise(n = n(), .groups = "drop")

hare_yearly_sf <- sites_sf %>%
  left_join(
    hare_yearly %>% st_drop_geometry(),
    by = "station"
  ) %>%
  mutate(n = replace_na(n, 0))

ggplot() +
  layer_spatial(osm_tiles) +
  geom_sf(data = sites_sf, size = 1, color = "grey60") +
  geom_sf(data = filter(hare_yearly_sf, year == 2020),
          aes(color = n, size = n)) +
  scale_color_viridis_c() +
  coord_sf() +
  labs(title = "Snowshoe hare Detections — 2020")

ggplot(hare_yearly_sf) +
  layer_spatial(osm_tiles) +
  geom_sf(aes(color = n, size = n)) +
  scale_color_viridis_c() +
  coord_sf() +
  facet_wrap(~ year) +
  labs(title = "Snowshoe hare Detections by Year")


#Overlap between bobcat and lynx

overlap_sf <- sites_sf %>%                     # KEEP geometry here
  left_join(
    sites_det %>%
      st_drop_geometry() %>%                   # drop geometry only from detections
      filter(species %in% c("bobcat", "canada lynx")) %>%
      select(station, species, n) %>%
      tidyr::pivot_wider(
        names_from = species,
        values_from = n,
        values_fill = 0
      ),
    by = "station"
  ) %>%
  mutate(
    bobcat = replace_na(bobcat, 0),
    `canada lynx` = replace_na(`canada lynx`, 0),
    overlap = case_when(
      bobcat > 0 & `canada lynx` > 0 ~ "Both",
      bobcat > 0 ~ "Bobcat only",
      `canada lynx` > 0 ~ "Lynx only",
      TRUE ~ "None"
    )
  )

ggplot() +
  layer_spatial(osm_tiles) +
  geom_sf(data = sites_sf, size = 1, color = "grey60") +
  geom_sf(data = overlap_sf, aes(color = overlap), size = 3) +
  scale_color_manual(
    values = c(
      "Bobcat only" = "#1f78b4",
      "Lynx only"   = "#33a02c",
      "Both"        = "#e31a1c",
      "None"        = "grey80"
    )
  ) +
  coord_sf() +
  labs(title = "Bobcat–Lynx Overlap at Camera Sites")



ggplot() +
  layer_spatial(osm_tiles) +
  geom_sf(data = sites_sf, size = 1, color = "grey60") +
  geom_sf(data = overlap_sf, aes(color = overlap), size = 3) +
  scale_color_manual(
    values = c(
      "Bobcat only" = "#1f78b4",
      "Lynx only"   = "#33a02c",
      "Both"        = "#e31a1c"
    )
  ) +
  coord_sf() +
  labs(title = "Bobcat–Lynx Overlap at Camera Sites")

#Load in climate data
climate_files <- list.files("~/chilcotin_bobcats/data/climate_data", pattern = "csv$", full.names = TRUE)

climate_covars <- lapply(climate_files, read.csv)

names(climate_covars) <- c("climatecovars_2018",
                           "climatecovars_2019",
                           "climatecovars_2020",
                           "climatecovars_2021",
                           "climatecovars_2022",
                           "climatecovars_2023",
                           "climatecovars_2024",
                           "climatecovars_2025")


##Format Climate into long format with column for year

climate_panel <- climate_covars %>%
  imap(~ .x %>%
         mutate(
           year = as.integer(str_extract(.y, "\\d{4}")),  # extract 2018, 2019, etc.
           station = ID1                                  # rename ID1 to station
         )) %>%
  bind_rows()

snow_plot <- ggplot(climate_panel, aes(x = year, y = PAS, color = station))+
  geom_point()
  

snow_plot

##Compare climatic covariates vs site covariates
climate_site_covars <- climate_panel %>%
  left_join(site_covars, by = "station")

##Elevation x Snow
ggplot(climate_site_covars, aes(x = year, y = PAS, color = elev))+
  geom_point()+
  scale_color_viridis_c()

##Elevation x Winter Temp
ggplot(climate_site_covars, aes(x = year, y = MCMT, color = elev))+
  geom_point()+
  scale_color_viridis_c()

##Elevation x Summer Temp
ggplot(climate_site_covars, aes(x = year, y = MWMT, color = elev))+
  geom_point()+
  scale_color_viridis_c()

##Format detections into weekly
det_weekly <- all_detections %>%
  mutate(week = isoweek(date.time)) %>%
  group_by(station, species, year, week) %>%
  summarise(detections = n(), .groups = "drop")

##Join climate data

det_weekly <- det_weekly %>%
  left_join(climate_panel, by = c("station", "year"))

##Elevation x detection
ggplot(det_weekly, aes(x = el, y = detections)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE, color = "black") +
  facet_grid(year ~ species, scales = "free_y") +
  theme_bw()

elev_q <- det_weekly %>%
  group_by(species, year) %>%
  summarise(
    q10 = quantile(el, 0.10, na.rm = TRUE),
    q50 = quantile(el, 0.50, na.rm = TRUE),
    q90 = quantile(el, 0.90, na.rm = TRUE)
  )

ggplot(elev_q, aes(x = year, y = q50, color = species)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = q10, ymax = q90, fill = species), alpha = 0.2, color = NA) +
  theme_bw()

#Bobcat climate effects
bobcat <- det_weekly %>%
  filter(species == "bobcat")

ggplot(bobcat, aes(x = PAS, y = detections, color = factor(year))) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_bw()

##Snowshoe hare over time
det_weekly <- det_weekly %>%
  mutate(
    iso_week = sprintf("%d-W%02d-1", year, week),   # Monday of each ISO week
    date = ISOweek2date(iso_week)
  )


ggplot(det_weekly, aes(x = date, y = detections))+
  geom_point()+
  geom_line()

##Annual mean detection rate with snowpack

bobcat_year <- bobcat %>%
  group_by(year) %>%
  summarise(
    mean_det = mean(detections),
    snow = mean(PAS),
    temp = mean(MCMT)
  )

ggplot(bobcat_year, aes(x = snow, y = mean_det)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw()

ggplot(bobcat_year, aes(x = temp, y = mean_det)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw()


bobcat_year$snow_bin <- cut(bobcat_year$snow, breaks = 3)

ggplot(bobcat, aes(x = week, y = detections, color = factor(year))) +
  geom_line(alpha = 0.5) +
  facet_wrap(~ snow_bin) +
  theme_bw()

lynx<-det_weekly %>%
  filter(species == "canada lynx")

ggplot(lynx, aes(x = PAS, y = detections, color = factor(year))) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_bw()

lynx_year <- lynx %>%
  group_by(year) %>%
  summarise(
    mean_det = mean(detections),
    snow = mean(PAS),
    temp = mean(MCMT)
  )

ggplot(lynx_year, aes(x = snow, y = mean_det)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw()
