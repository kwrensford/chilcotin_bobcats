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


#Read in camera detections
all_detections <-read.csv("~/chilcotin_bobcats/data/kwasi.bobcats.etc.csv")

#Read in environmental covariates
site_covars <- read.csv("~/chilcotin_bobcats/data/env.vars.2021.csv")

#Format date/time
all_detections$date.time <- parse_date_time(all_detections$date.time, "ymd HMS", tz = "Canada/Pacific")

all_detections <- all_detections %>%
  mutate(
    date.time = parse_date_time(date.time, orders = c("ymd HMS", "ymd HM", "ymd H"))
  )

all_detections <- all_detections %>% filter(!is.na(date.time))


#Join LatLong and Elevation to detections

##Rename sites to station in covars

site_covars <- site_covars %>%
  rename(station = site)

all_detections <- all_detections %>%
  left_join(site_covars, by = "station")

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
    snow = mean(PAS)
  )

ggplot(bobcat_year, aes(x = snow, y = mean_det)) +
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
    snow = mean(PAS_wt)
  )

ggplot(lynx_year, aes(x = snow, y = mean_det)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw()
