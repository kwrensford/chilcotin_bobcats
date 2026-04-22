library(dplyr)
library(tidyr)
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

bobcats <- read.csv("~/chilcotin_bobcats/data/bobcats.2018.2025.csv")

bobcats$date.time <- parse_date_time(bobcats$date.time, "mdy HM", tz = "Canada/Pacific")

bobcats$year <- year(bobcats$date.time)

bobcat_counts_annual <- bobcats %>%
  count(year)

bobcat_counts_station <- bobcats %>%
  count(station)

ggplot(bobcat_counts_annual, aes(x = year, y = n)) +
  geom_col(fill = "steelblue") +
  labs(x = "Year",
       y = "Number of detections",
       title = "Detections per year") +
  theme_minimal()



bobcat_year_station <- bobcats %>%
  count(year, station)

ggplot(bobcat_year_station, aes(x = year, y = n, fill = station)) +
  geom_col(position = "stack") +
  labs(x = "Year",
       y = "Number of detections",
       fill = "Station",
       title = "Detections per year by station") +
  theme_minimal()
