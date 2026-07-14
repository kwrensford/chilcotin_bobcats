#Exploring independent detections

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
library(data.table)

#Read in camera detections
all_detections <-read.csv("~/chilcotin_bobcats/data/kwasi.bobcats.etc.csv")

#Format date/time

all_detections$date.time <- parse_date_time(
  all_detections$date.time,
  orders = c("mdy HMS", "mdy HM", "mdy"),
  tz = "Canada/Pacific"
)

#Generate delta time (intervals between detection)
det_intervals <- all_detections %>%
  arrange(species, station, date.time) %>%
  group_by(species, station) %>%
  mutate(delta_time = as.numeric(difftime(date.time, lag(date.time), units = "mins"))) %>%
  ungroup()

#Visualize detection intervals for each species

ggplot(det_intervals, aes(x = delta_time)) + 
  geom_histogram() + 
  xlim(NA, 100)
  facet_wrap(~ species)
