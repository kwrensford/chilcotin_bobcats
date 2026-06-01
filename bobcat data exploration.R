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

#Read in camera detections
all_detections <-read.csv("~/chilcotin_bobcats/data/kwasi.bobcats.etc.csv")

#Read in environmental covariates
site_covars <- read.csv("~/chilcotin_bobcats/data/env.vars.2021.csv")

#Read in camera deployment data

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
  
  # ---- Normalize names ----
  all_intervals <- all_intervals %>%
    rename(
      Camera = any_of(c("Camera", "camera", "Camera.x", "Station")),
      Year   = any_of(c("Year", "year", "Session"))
    )
  
  # ---- Fix problem column names ----
  all_intervals <- all_intervals %>%
    rename_with(~ gsub("Problem([0-9]+)_Problem_from", "Problem\\1_from", .x)) %>%
    rename_with(~ gsub("Problem([0-9]+)_Problem_to",   "Problem\\1_to",   .x))
  
  # ---- Build multiseason camOp ----
  camop <- cameraOperation(
    CTtable      = all_intervals,
    stationCol   = "Camera",
    cameraCol    = "Camera",
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

#Format date/time
all_detections$date.time <- parse_date_time(all_detections$date.time, "ymd HMS", tz = "Canada/Pacific")

all_detections <- all_detections %>%
  mutate(
    date.time = parse_date_time(date.time, orders = c("ymd HMS", "ymd HM", "ymd H"))
  )

all_detections <- all_detections %>% filter(!is.na(date.time))


#Join LatLong and Elevation to detections
all_detections$site <- all_detections$station

all_detections <- all_detections %>%
  left_join(site_covars, by = "site")

##Detection histories

build_multiseason_detection_history_camtrapR <- function(detections, camop) {
  
  # ---- Normalize station names to match camOp rownames ----
  normalize_station <- function(x) {
    x <- sub("_2$", "", x)      # remove _2 suffix
    x <- sub("off$", "", x)     # remove off suffix
    x <- sub("err$", "", x)     # remove err suffix
    x
  }
  
  detections <- detections %>%
    mutate(
      date.time = as.POSIXct(date.time),
      station_clean = normalize_station(station)
    )
  
  species_list <- sort(unique(detections$species))
  seasons <- names(camop)
  
  out <- list()
  
  for (sp in species_list) {
    
    sp_list <- list()
    
    for (ss in seasons) {
      
      camop_ss <- camop[[ss]]
      
      # Filter detections for this species × season
      det_ss <- detections %>%
        filter(species == sp, year == as.integer(ss)) %>%
        filter(!is.na(date.time))   # <-- NEW: drop NA timestamps
      
      
      # Normalize station names in this subset
      det_ss$station_clean <- normalize_station(det_ss$station)
      
      # ---- Check for station mismatches ----
      missing <- setdiff(det_ss$station_clean, rownames(camop_ss))
      if (length(missing) > 0) {
        message("WARNING: Season ", ss, " | Species ", sp)
        message("Stations in detections but NOT in camOp: ", paste(missing, collapse=", "))
        
        # Drop unmatched stations (camtrapR cannot use them)
        det_ss <- det_ss %>% filter(!station_clean %in% missing)
      }
      
      # If no detections remain, return an all-zero matrix
      if (nrow(det_ss) == 0) {
        sp_list[[ss]] <- matrix(
          0,
          nrow = nrow(camop_ss),
          ncol = ncol(camop_ss),
          dimnames = dimnames(camop_ss)
        )
        next
      }
      
      # ---- Build detection history using camtrapR ----
      dh <- detectionHistory(
        recordTable = det_ss,
        species = sp,                      # REQUIRED in your version
        camOp = camop_ss,
        stationCol = "station_clean",
        speciesCol = "species",
        recordDateTimeCol = "date.time",
        occasionLength = 1,
        day1 = "station",
        includeEffort = FALSE
      )
      
      sp_list[[ss]] <- dh$detection_history
    }
    
    out[[sp]] <- sp_list
  }
  
  return(out)
}



##SP Occupancy


convert_to_spoccupancy_array <- function(det_list) {
  
  seasons <- names(det_list)
  n_years <- length(seasons)
  
  # ---- 1. Determine full station set ----
  all_stations <- Reduce(union, lapply(det_list, rownames))
  
  # ---- 2. Determine full occasion set (days) ----
  all_occasions <- Reduce(union, lapply(det_list, colnames))
  
  # ---- 3. Pad each season matrix to full dimensions ----
  padded <- lapply(det_list, function(mat) {
    
    # Create full matrix
    full <- matrix(
      0,
      nrow = length(all_stations),
      ncol = length(all_occasions),
      dimnames = list(all_stations, all_occasions)
    )
    
    # Insert existing data
    full[rownames(mat), colnames(mat)] <- mat
    
    return(full)
  })
  
  # ---- 4. Build 3D array ----
  y <- array(
    NA,
    dim = c(length(all_stations), length(all_occasions), n_years),
    dimnames = list(all_stations, all_occasions, seasons)
  )
  
  for (i in seq_along(seasons)) {
    y[,,i] <- padded[[ seasons[i] ]]
  }
  
  return(y)
}

convert_to_spoccupancy_array <- function(det_list) {
  
  seasons <- names(det_list)
  n_years <- length(seasons)
  
  # ---- 1. Determine full station set ----
  all_stations <- Reduce(union, lapply(det_list, rownames))
  
  # ---- 2. Determine full occasion set (days) ----
  all_occasions <- Reduce(union, lapply(det_list, colnames))
  
  # ---- 3. Pad each season matrix to full dimensions ----
  padded <- lapply(det_list, function(mat) {
    
    # Create full matrix
    full <- matrix(
      0,
      nrow = length(all_stations),
      ncol = length(all_occasions),
      dimnames = list(all_stations, all_occasions)
    )
    
    # Insert existing data
    full[rownames(mat), colnames(mat)] <- mat
    
    return(full)
  })
  
  # ---- 4. Build 3D array ----
  y <- array(
    NA,
    dim = c(length(all_stations), length(all_occasions), n_years),
    dimnames = list(all_stations, all_occasions, seasons)
  )
  
  for (i in seq_along(seasons)) {
    y[,,i] <- padded[[ seasons[i] ]]
  }
  
  return(y)
}

all_arrays <- build_all_detection_arrays(detHist)

##SPOccupancyt Detection Array

convert_to_spoccupancy_array <- function(det_list) {
  
  seasons <- names(det_list)
  n_years <- length(seasons)
  
  # ---- 1. Determine full station set ----
  all_stations <- Reduce(union, lapply(det_list, rownames))
  
  # ---- 2. Determine full occasion set (days) ----
  all_occasions <- Reduce(union, lapply(det_list, colnames))
  
  # ---- 3. Pad each season matrix to full dimensions ----
  padded <- lapply(det_list, function(mat) {
    
    # Create full matrix
    full <- matrix(
      0,
      nrow = length(all_stations),
      ncol = length(all_occasions),
      dimnames = list(all_stations, all_occasions)
    )
    
    # Insert existing data
    full[rownames(mat), colnames(mat)] <- mat
    
    return(full)
  })
  
  # ---- 4. Build 3D array ----
  y <- array(
    NA,
    dim = c(length(all_stations), length(all_occasions), n_years),
    dimnames = list(all_stations, all_occasions, seasons)
  )
  
  for (i in seq_along(seasons)) {
    y[,,i] <- padded[[ seasons[i] ]]
  }
  
  return(y)
}

#Align covariates with detection array
normalize_station <- function(x) {
  x <- sub("_2$", "", x)
  x <- sub("off$", "", x)
  x <- sub("err$", "", x)
  x
}

site_covs <- site_covariates_raw %>%
  mutate(station = normalize_station(station)) %>%   # same normalizer as detections
  filter(station %in% stations) %>%
  distinct(station, .keep_all = TRUE) %>%
  right_join(
    tibble(station = stations),
    by = "station"
  ) %>%
  arrange(match(station, stations))


site_covs <- site_covs %>% column_to_rownames("station")

##Diel activity

activityDensity(
  recordTable = all_detections,
  species = "red squirrel",
  speciesCol = "species",
  recordDateTimeCol = "date.time",
  recordDateTimeFormat = "ymd HMS",
  writePNG = FALSE
)
