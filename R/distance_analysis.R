#' @title Viewshed Distance Analysis
#' @description Since visibility changes with increasing distance, a threshold for the study area can be evaluated beyond which an observer can't see any terrain.
#' The \code{distance_analysis} function calculates the proportion of visible area for each distance value.
#'
#' @param observer object of class \code{sf} with POINT geometries; Observer location(s) from where the distance analysis should be computed. It's best to provide a sample of points within the complete study area, to evaluate the Visibility threshold of the study area.
#' @param max_distance numeric; Maximum buffer distance
#' @param dsm_rast object of class \code{\link[terra]{rast}}; \code{\link[terra]{rast}} of the DSM
#' @param dtm_rast object of class \code{\link[terra]{rast}}; \code{\link[terra]{rast}} of the DTM
#' @param observer_height numeric > 0; Height of the observer (e.g. 1.7 meters)
#' @param raster_res optional; NULL or numeric > 0; Resolution that the viewshed raster should be aggregated to. Must be a multible of the dsm_rast resolution
#' @param plot optional; Plot (Visibility ~ Distance)
#' @param progress logical; Show progress bar and computation time?
#' @param cores numeric; The number of cores to use
#'
#' @return object of class \code{\link[tibble]{tibble}}
#' @export
#' 
#' @importFrom magrittr %>%
#' @importFrom sf st_buffer
#' @importFrom sf st_coordinates
#' @importFrom sf st_crs
#' @importFrom sf st_geometry_type
#' @importFrom sf st_make_grid
#' @importFrom terra crs
#' @importFrom terra extract
#' @importFrom terra res
#' @importFrom terra crop
#' @importFrom terra mask
#' @importFrom terra vect
#' @importFrom terra aggregate
#' @importFrom terra rowFromY
#' @importFrom terra colFromX
#' @importFrom terra values
#' @importFrom terra ncol
#' @importFrom terra boundaries
#' @importFrom terra xyFromCell
#' @importFrom terra rast
#' @importFrom raster raster
#' @importFrom dplyr as_tibble
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_smooth
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_light
#' @importFrom scales percent_format

distance_analysis <- function(observer, dsm_rast, dtm_rast, 
                     max_distance = 800, observer_height = 1.7, 
                     raster_res = NULL, plot = FALSE,
                     progress = FALSE, cores = 1) {
  #### 1. Check input ####
  # observer
  if (!is(observer, "sf")) {
    stop("observer must be a sf object")
  } else if (sf::st_crs(observer)$units != "m") {
    stop("observer CRS unit needs to be metric")
  } else if (!as.character(sf::st_geometry_type(observer, by_geometry = FALSE)) %in% c("POINT", "MULTIPOINT")) {
    stop("observer must be POINT")
  }
  
  # dsm_rast
  if (!is(dsm_rast, "SpatRaster")) {
    stop("dsm_rast must be a SpatRaster object")
  } else if (sf::st_crs(terra::crs(dsm_rast))$epsg != sf::st_crs(observer)$epsg) {
    stop("dsm_rast must have the same CRS as observer")
  }
  
  # dtm_rast
  if (!is(dtm_rast, "SpatRaster")) {
    stop("dtm_rast must be a SpatRaster object")
  } else if (sf::st_crs(terra::crs(dtm_rast))$epsg != sf::st_crs(observer)$epsg) {
    stop("dtm_rast must have the same CRS as observer")
  }
  
  # max_distance
  max_distance <- round(max_distance, digits = 0)
  
  # raster_res
  dsm_res <- min(terra::res(dsm_rast))
  if (is.null(raster_res)) {
    raster_res = dsm_res
  } else if (raster_res < min(terra::res(dsm_rast))) {
    stop("raster_res must be higher than the resolution of dsm_rast")
  } else if ((raster_res %% dsm_res) != 0) {
    stop(paste0("raster_res must be a multible of the dsm_rast resolution. Try raster_res = ", raster_res - (raster_res %% dsm_res)))
  }
  rm(dsm_res)
  
  
  #### 2. Prepare Data for viewshed analysis ####
  # Max AOI
  max_aoi <- observer %>% 
    sf::st_bbox() %>% 
    sf::st_as_sfc() %>% 
    sf::st_buffer(max_distance) %>% 
    sf::st_as_sf()
  
  max_aoi <- terra::vect(max_aoi) %>% 
    terra::crop(dsm_rast)
  
  #### 3. Check RAM size ####
  grid_fact <- rast_fits_vect_fact(max_aoi = max_aoi, max_distance = max_distance, raster_res = raster_res)
  
  aoi_grid <- sf::st_make_grid(max_aoi, n = ifelse(grid_fact == 1, 1, grid_fact*4))
  
  
  #### 4. Main loop ####
  distance_tbl <- tibble(
    Distance = as.integer(),
    Visible_perc = as.double()
  )
  
  if(progress) {
    message("Progress:")
    
    if(length(aoi_grid) > 1) {
      pb = txtProgressBar(min = 0, max = length(aoi_grid), initial = 0, style = 3)
    }
    start_time <- Sys.time()
  }
  for (i in seq_along(aoi_grid)) {
    if ( progress && (length(aoi_grid) > 1) ) setTxtProgressBar(pb, i)
    
    this_aoi <- aoi_grid[i]
    this_observer <- observer[this_aoi,]
    
    #### 1. Crop DSM to max AOI and change resolution ####
    this_rast <- dsm_rast %>% 
      terra::crop(terra::vect(sf::st_buffer(this_aoi, max_distance)))
    
    if(raster_res != min(raster::res(this_rast))) {
      terra::terraOptions(progress = 0)
      this_rast <- terra::aggregate(this_rast, fact = raster_res/terra::res(this_rast))
      terra::terraOptions(progress = 3)
    }
    
    # Convert to vector
    dsm_vec <- terra::values(this_rast, mat = FALSE)
    dsm_cpp_rast <- this_rast %>% terra::rast() %>% raster::raster()
    
    # Coordinates of start point
    x0 <- sf::st_coordinates(this_observer)[,1]
    y0 <- sf::st_coordinates(this_observer)[,2]
    
    # Observer heights
    height_0_vec <- unlist(terra::extract(dtm_rast, cbind(x0, y0)), use.names = F) + observer_height
    
    #### 2. Remove points outside the DSM or DTM ####
    invalid_points <- unique(c(
      which(is.na(terra::extract(this_rast, cbind(x0, y0)))), # points outside the DSM
      which(is.na(height_0_vec)) # points outside the DTM
    ))
    
    # Remove invalid points
    if (length(invalid_points) > 0) {
      observer <- observer[-invalid_points, ]
      x0 <- x0[-invalid_points]
      y0 <- y0[-invalid_points]
      height_0_vec <- height_0_vec[-invalid_points]
    }
    
    # Skip if all values are NaN
    if(length(x0) == 0) next
    
    #### 3. Compute viewshed ####
    # Start row/col
    r0 <- terra::rowFromY(this_rast, y0)
    c0 <- terra::colFromX(this_rast, x0)
    
    this_distance_tbl <- viewshed_distance_analysis_cpp(dsm_cpp_rast, dsm_vec,
                                                        c0, r0, max_distance, height_0_vec,
                                                        ifelse(length(x0) > cores, cores, length(x0)),
                                                        ifelse(length(aoi_grid) > 1, FALSE, progress))
    colnames(this_distance_tbl) <- c("V1", "V2", "V3")
    
    this_distance_tbl <- this_distance_tbl %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(Visible_perc = V3 / V2) %>%
      dplyr::rename(Distance = V1) %>%
      dplyr::select(Distance, Visible_perc)
    
    distance_tbl <- dplyr::add_row(distance_tbl, this_distance_tbl)
    rm(dsm_vec); invisible(gc())
  }
  
  if (progress) {
    time_dif <- round(cores * ((as.numeric(difftime(Sys.time(), start_time, units = "s"))*1000) / nrow(observer)), 2)
    cat("\n")
    
    time_total <- round(as.numeric(difftime(Sys.time(), start_time, units = "m")))
    if(time_total < 1){
      time_total <- round(as.numeric(difftime(Sys.time(), start_time, units = "s")))
      
      if(time_total < 1){
        time_total <- round(as.numeric(difftime(Sys.time(), start_time, units = "s")))*1000
        message(paste("Total runtime:", time_total, " milliseconds"))
      } else {
        message(paste("Total runtime:", time_total, " seconds"))
      }
    } else {
      message(paste("Total runtime:", time_total, " minutes"))
    }
    
    message(paste("Average time for a single point:", time_dif, "milliseconds"))
  }
  return(distance_tbl)
}