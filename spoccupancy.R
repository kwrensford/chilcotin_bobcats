#SP Occupancy Analyses

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
library(stars)
library(concaveman)

#Read in detection data

#Read in covariate data

#Single Species single season model

##Read in detection data for 2025 bobcats
y <- detHist$bobcat$`2025`
dim(y)

##Site covariates
det_sites <- rownames(y)
cov_sites <- site_covars$station
common_sites <- intersect(det_sites, cov_sites)

y <- y[common_sites, ]

site_covs_raw <- site_covars 

site_covs <- site_covs_raw[site_covs_raw$station %in% common_sites, ]
site_covs <- site_covs[match(common_sites, site_covs$station), ]

rownames(site_covs) <- site_covs$station

#Annual climate covariate
annual_covs_2025 <- climate_covars$climatecovars_2025

annual_covs_2025 <- annual_covs_2025 %>%
  rename(station = ID1)

annual_covs_2025 <- annual_covs_2025[annual_covs_2025$station %in% common_sites, ]
annual_covs_2025 <- annual_covs_2025[match(common_sites, annual_covs_2025$station)]

rownames(annual_covs_2025) <- annual_covs_2025$station

#Coordinates
site_coords <- site_covs %>%
  select(station, lat, lon)

coords <- as.matrix(site_coords[, c("lat", "lon")])

identical(rownames(y), rownames(site_covs))
identical(rownames(y), rownames(annual_covs_2025))
identical(rownames(y), rownames(coords))

#Build spOccupancy data list
site_covs$PAS <- annual_covs_2025$PAS

data_list <- list(
  y = y,
  occ.covs = site_covs,
  coords = coords
)

bobcat_model <- spPGOcc(
  occ.formula = ~ elev + PAS,
  det.formula = ~ 1,
  data = data_list,
  n.batch = 3000,
  batch.length = 25,
  n.burn = 3000,
  n.chains = 3
)


summary(bobcat_model)

#Predict and plot
coords_sf <- st_as_sf(site_coords, coords = c("lon", "lat"), crs = 4326)
coords_sf <- st_transform(coords_sf, 3005)

study_area<-concaveman(coords_sf)
st_crs(study_area) <- 3005 

study_area <- st_convex_hull(st_union(coords_sf))

plot(study_area)

grid <- st_make_grid(
  study_area,
  cellsize = 500,
  what = "centers"
) %>%
  st_as_sf() %>%
  st_intersection(study_area)
#Join covariates
coords_sf <- coords_sf %>% 
  left_join(
  site_covs,
  by = "station"
)
#Convert spatial sf objects to spat vectors
coords_v <- vect(coords_sf)
grid_v <- vect(grid)

#Compute nearest-neighbor covariates
##Compute distance matrix
dmat <- distance(grid_v, coords_v)

#Identify nearest camera for each grid cell
nearest_idx <- apply(dmat, 1, which.min)

#Assign covariates from nearest camera
grid$elev <- coords_sf$elev[nearest_idx]
grid$PAS  <- coords_sf$PAS[nearest_idx]

#Standardize using model means
elev_mean <- mean(site_covs$elev)
elev_sd   <- sd(site_covs$elev)

PAS_mean  <- mean(site_covs$PAS)
PAS_sd    <- sd(site_covs$PAS)

grid$elev_std <- (grid$elev - elev_mean) / elev_sd
grid$PAS_std  <- (grid$PAS  - PAS_mean)  / PAS_sd

occ_covs_pred <- data.frame(
  elev = grid$elev_std,
  PAS  = grid$PAS_std
)


#Predict Occupancy
coords_pred <- st_coordinates(grid)

X.0 <- cbind(
  1,
  grid$elev_std,
  grid$PAS_std
)


pred <- predict(
  bobcat_model,
  X.0 = X.0,
  coords.0 = coords_pred,   # prediction coordinates
  type = "occupancy"
)


X0 <- cbind()
bobcat_model_pred <- predict(bobcat_model, X0)


plot_data_bobcat <- data.frame(x = coords_sf$lon,
                               y = coords_sf$lat,
                               mean_psi = apply(pred$psi.0.samples, 2, mean),
                               sd.psi = apply(pred$psi.0.samples, 2, sd))

dat_stars <- st_as_stars(plot_data_bobcat, dims = c('x' , 'y'))

ggplot() + 
  geom_stars(data = dat_stars, aes(x = x, y = y, fill = mean_psi)) + 
  scale_fill_viridis_c(na.value = "transparent") + 
  labs(x = "Longitude", y = "Latitude", fill = "",
       title = "Mean bobcat occurence probability") + 
  theme_bw()

##SPOccupancy Detection Array

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

site_covs <- site_covars %>%
  mutate(station = normalize_station(station)) %>%   # same normalizer as detections
  filter(station %in% stations) %>%
  distinct(station, .keep_all = TRUE) %>%
  right_join(
    tibble(station = stations),
    by = "station"
  ) %>%
  arrange(match(station, stations))


site_covs <- site_covs %>% column_to_rownames("station")

#Model syntax (Bobcat)

bobcat_occ_formula <- ~ scale(year) + scale(elev) + scale(PAS)
bobcat_det_formula <- ~ 

stPGOcc()

#Model syntax (Lynx)

lynx_occ_formula <- ~ scale(year) + scale(elev) + scale(PAS)
lynx_det_formula <- ~ scale()

lynx_model <- stPGOcc(occ.formula = ,
                      det.formula = lynx_det_formula,
                      data = lynx_data,
                      inits = lynx_sp_inits,
                      priors = lynx_sp_priors,
                      cov.model = 