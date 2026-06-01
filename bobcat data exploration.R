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

#Read in camera detections
all_detections <-read.csv("~/chilcotin_bobcats/data/kwasi.bobcats.etc.csv")

#Read in environmental covariates
site_covars <- read.csv("~/chilcotin_bobcats/data/env.vars.2021.csv")

#Read in camera deployment data

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(camtrapR)

build_multiseason_camop_with_problems <- function(files) {
  
  process_one <- function(path) {
    
    hdr <- read_csv(path, n_max = 2, col_names = FALSE, show_col_types = FALSE)
    row1 <- as.character(hdr[1, ])
    row2 <- as.character(hdr[2, ])
    clean_names <- ifelse(is.na(row1) | row1 == "", row2, row1)
    
    dat <- read_csv(path, skip = 2, col_names = clean_names, show_col_types = FALSE)
    
    date_cols <- grep("^\\d{4}-\\d{2}-\\d{2}$", names(dat), value = TRUE)
    year_val <- str_extract(basename(path), "\\d{4}") %>% as.integer()
    
    deploy_long <- dat %>%
      pivot_longer(cols = all_of(date_cols),
                   names_to = "date",
                   values_to = "operational") %>%
      mutate(date = ymd(date),
             Year = year_val)
    
    # Remove camera-year groups with no operational data
    deploy_long <- deploy_long %>%
      group_by(Camera, Year) %>%
      filter(!all(is.na(operational))) %>%
      ungroup()
    
    gaps <- deploy_long %>%
      arrange(Camera, date) %>%
      group_by(Camera, Year) %>%
      mutate(is_gap = operational == 0,
             change = is_gap != lag(is_gap, default = FALSE),
             block = cumsum(change)) %>%
      filter(is_gap) %>%
      group_by(Camera, Year, block) %>%
      summarise(Problem_from = min(date),
                Problem_to   = max(date),
                .groups = "drop")
    
    if (nrow(gaps) > 0) {
      gaps <- gaps %>%
        group_by(Camera, Year) %>%
        mutate(idx = row_number()) %>%
        pivot_wider(id_cols = c(Camera, Year),
                    names_from = idx,
                    values_from = c(Problem_from, Problem_to),
                    names_glue = "Problem{idx}_{.value}")
    }
    
    deploy_intervals <- deploy_long %>%
      group_by(Camera, Year) %>%
      summarise(SetupDate = min(date),
                RetrievalDate = max(date),
                .groups = "drop")
    
    if (nrow(gaps) > 0) {
      deploy_intervals <- left_join(deploy_intervals, gaps,
                                    by = c("Camera", "Year"))
    }
    
    return(deploy_intervals)
  }
  
  all_intervals <- bind_rows(lapply(files, process_one))
  
  # ---- DIAGNOSTIC OUTPUT ----
  message("Columns in CTtable BEFORE renaming:")
  print(names(all_intervals))
  message("First few rows of CTtable:")
  print(head(all_intervals))
  message("---- END DIAGNOSTIC ----")
  
  # ---- Try to normalize names ----
  all_intervals <- all_intervals %>%
    rename(
      Camera = any_of(c("Camera", "camera", "Camera.x", "Station")),
      Year   = any_of(c("Year", "year", "Session"))
    )
  
  message("Columns in CTtable AFTER renaming:")
  print(names(all_intervals))
  message("First few rows AFTER renaming:")
  print(head(all_intervals))
  message("---- END DIAGNOSTIC ----")
  
  # ---- Now call cameraOperation ----
  camop <- cameraOperation(
    CTtable      = all_intervals,
    stationCol   = "Camera",
    cameraCol    = "Camera",        # <-- REQUIRED
    setupCol     = "SetupDate",
    retrievalCol = "RetrievalDate",
    sessionCol   = "Year",
    byCamera     = TRUE,
    hasProblems  = TRUE,
    dateFormat   = "%Y-%m-%d",
    writecsv     = FALSE
  )
  
  return(list(
    CTtable = all_intervals,
    camop   = camop
  ))
}


files <- list.files("~/chilcotin_bobcats/data/deployment_data", pattern = "csv$", full.names = TRUE)
multi <- build_multiseason_camop_with_problems(files)


deploy_list <- lapply(files, read_deployment_csv)

all_intervals <- bind_rows(lapply(deploy_list, `[[`, "intervals"))

#Format deployment tables into long format

#Format date/time
all_detections$date.time <- parse_date_time(all_detections$date.time, "ymd HMS", tz = "Canada/Pacific")

#Join LatLong and Elevation to detections
all_detections$site <- all_detections$station

all_detections <- all_detections %>%
  left_join(site_covars, by = "site")

#Explore naive occupancy per species

##Bobcats
bobcat_site_occ <- all_detections %>%
  group_by(location, longitude, latitude)%>%
  summarise(coyotedetection = as.integer(any(species_common_name == "Coyote")))

naive_occ_coyote <- mean(coyote_site_occ$coyotedetection)

##Relationship between detections and covariates.

#Camera Operation Matrices
camops <- lapply(deploy_list, function(x) {
  cameraOperation(
    CTtable = deploy_intervals,
    stationCol = "Camera",
    setupCol = "SetupDate",
    retrievalCol = "RetrievalDate",
    sessionCol = "Year",     # <-- NEW
    cameraCol = "Camera",    # <-- REQUIRED when multiple rows per camera
    dateFormat = "%Y-%m-%d"
  )
})


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
