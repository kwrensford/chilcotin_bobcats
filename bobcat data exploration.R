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

##Splitting camop tables into stacked array
split_camop_by_year <- function(camop_matrix) {
  
  # Extract year from rownames
  row_years <- sub(".*__SESS_([0-9]{4})__.*", "\\1", rownames(camop_matrix))
  
  # Extract year from column names (YYYY-MM-DD)
  col_years <- substr(colnames(camop_matrix), 1, 4)
  
  years <- sort(unique(row_years))
  
  out <- list()
  
  for (yr in years) {
    
    row_idx <- which(row_years == yr)
    col_idx <- which(col_years == yr)
    
    m <- camop_matrix[row_idx, col_idx, drop = FALSE]
    m <- as.matrix(m)
    
    # Clean rownames
    station <- sub("__.*$", "", rownames(m))
    rownames(m) <- station
    
    out[[yr]] <- m
  }
  
  return(out)
}

clean_camop <- split_camop_by_year(multi$camop)

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
        occasionLength = 7,
        day1 = "station",
        includeEffort = FALSE,
        timeZone = "America/Vancouver"
      )
      
      sp_list[[ss]] <- dh$detection_history
    }
    
    out[[sp]] <- sp_list
  }
  
  return(out)
}

detHist <- build_multiseason_detection_history_camtrapR(
  detections = all_detections,
  camop = clean_camop
)

#Save detection matrices as csv's

iwalk(detHist, function(year_list, species) {
  iwalk(year_list, function(mat, yr) {
    
    # build path: data/species_year.csv
    out_path <- file.path("data", paste0(species, "_", yr, ".csv"))
    
    write.csv(mat, out_path, row.names = FALSE)
  })
})


##Diel activity

activityDensity(
  recordTable = all_detections,
  species = "snowshoe hare",
  speciesCol = "species",
  recordDateTimeCol = "date.time",
  recordDateTimeFormat = "ymd HMS",
  writePNG = FALSE
)

activityOverlap(recordTable = all_detections,
                speciesCol = "species",
                speciesA = "canada lynx",
                speciesB = "bobcat",
                recordDateTimeCol = "date.time",
                recordDateTimeFormat = "ymd HMS",
                writePNG = FALSE)

#Climate NA data
climate_years <- c("Year_2018.ann",
                   "Year_2019.ann", 
                   "Year_2020.ann", 
                   "Year_2021.ann", 
                   "Year_2022.ann",
                   "Year_2023.ann",
                   "Year_2024.ann",
                   "Year_2025.ann")

##To downloard: Apply across all years (Important, Climate NA only allows 5 calls per hour, so may need to break up into batches)
climate_list <- lapply(climate_years, function(yr){
  ClimateNA_API2(
    ClimateBC_NA = "BC",
    inputFile = "~/chilcotin_bobcats/data/climatena_coords.csv",
    period = yr,
    MSY = "SY"
  )
})
#When downloaded, move into a single folder within your data folder called "climate_data"

#If already downloaded, load in from data folder
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
