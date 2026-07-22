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

files <- list.files("~/chilcotin_bobcats/data/deployment_data", pattern = "csv$", full.names = TRUE)

deploy_list <- lapply(files, read.csv)

names(deploy_list) <- gsub("\\D", "", basename(files))

#Format date/time

all_detections$date.time <- parse_date_time(
  all_detections$date.time,
  orders = c("mdy HMS", "mdy HM", "mdy"),
  tz = "Canada/Pacific"
)

#Clean Deployment CSV's
read_deploy <- function(file) {
  
  # Read raw CSV with fill=TRUE so rows align
  raw <- read.csv(file, header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
  
  # Row 1: date names (for date columns)
  # Row 2: metadata names (for metadata columns)
  
  # Identify date columns: row 1 matches YYYY-MM-DD
  is_date <- grepl("^\\d{4}-\\d{2}-\\d{2}$", raw[1, ])
  
  # Build final column names:
  # - date columns get row 1 names
  # - metadata columns get row 2 names
  final_names <- ifelse(is_date, raw[1, ], raw[2, ])
  
  # Drop first two rows
  df <- raw[-c(1, 2), ]
  names(df) <- final_names
  
  df
}

deploy_list <- lapply(files, read_deploy)
names(deploy_list) <- gsub("\\D", "", basename(files))   # extract year



#Compute  start/end for each camera
get_start_end <- function(df, year) {
  
  # Identify date columns
  date_cols <- grep("^\\d{4}-\\d{2}-\\d{2}$", names(df), value = TRUE)
  
  # Convert to Date objects
  dates <- as.Date(date_cols)
  
  df %>%
    rowwise() %>%
    mutate(
      deploy_vec = list(c_across(all_of(date_cols))),
      start_date = if (any(unlist(deploy_vec) == 1)) dates[which(unlist(deploy_vec) == 1)[1]] else NA,
      end_date   = if (any(unlist(deploy_vec) == 1)) dates[tail(which(unlist(deploy_vec) == 1), 1)] else NA,
      year = year
    ) %>%
    ungroup() %>%
    select(Station, year, start_date, end_date)
}

results <- imap_dfr(deploy_list, ~ get_start_end(.x, .y))

library(tidyverse)

plot_deploy <- function(df) {
  
  # Identify date columns
  date_cols <- grep("^\\d{4}-\\d{2}-\\d{2}$", names(df), value = TRUE)
  
  df_long <- df %>%
    pivot_longer(
      cols = all_of(date_cols),
      names_to = "date",
      values_to = "active"
    ) %>%
    mutate(date = as.Date(date))
  
  ggplot(df_long, aes(x = date, y = Station, fill = active)) +
    geom_tile() +
    scale_fill_manual(values = c("0" = "grey90", "1" = "steelblue")) +
    theme_minimal() +
    labs(
      title = "Camera Deployment Activity",
      x = "Date",
      y = "Station",
      fill = "Active"
    )
}

ggplot(results, aes(y = Station)) +
  geom_segment(
    aes(x = start_date, xend = end_date, color = factor(year)),
    size = 3
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  labs(
    title = "Multi‑Year Camera Deployment Timeline",
    x = "Date",
    y = "Station",
    color = "Year"
  )


plot_deploy(deploy_list[["2018"]])
