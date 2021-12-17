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
  if(progress) {
    message("Preprocessing:")
    pb = txtProgressBar(min = 0, max = 4, initial = 0, style = 2)
  }
  # Max AOI
  max_aoi <- observer %>% 
    sf::st_bbox() %>% 
    sf::st_as_sfc() %>% 
    sf::st_buffer(max_distance)
  
  max_aoi <- terra::vect(max_aoi) %>% 
    terra::crop(dsm_rast)
  
  if (progress) setTxtProgressBar(pb, 1)
  
  # Crop DSM to max AOI and change resolution
  dsm_rast <- terra::crop(dsm_rast, max_aoi)
  
  
  if(raster_res != min(raster::res(dsm_rast))) {
    terra::terraOptions(progress = 0)
    dsm_rast <- terra::aggregate(dsm_rast, fact = raster_res/terra::res(dsm_rast))
    terra::terraOptions(progress = 3)
  }
  
  if (progress) setTxtProgressBar(pb, 2)
  
  dsm_vec <- terra::values(dsm_rast, mat = FALSE)
  dsm_cpp_rast <- dsm_rast %>% terra::rast() %>% raster::raster()
  
  if (progress) setTxtProgressBar(pb, 3)
  
  # Coordinates of start point
  x0 <- sf::st_coordinates(observer)[,1]
  y0 <- sf::st_coordinates(observer)[,2]
  
  # Observer heights
  height_0_vec <- unlist(terra::extract(dtm_rast, cbind(x0, y0)), use.names = F) + observer_height
  
  
  #### 3. Remove points outside the DSM or DTM ####
  invalid_points <- unique(c(
    which(is.na(terra::extract(dsm_rast, cbind(x0, y0)))), # points outside the DSM
    which(is.na(height_0_vec)) # points outside the DTM
  ))
  
  # Remove invalid points
  if (length(invalid_points) > 0) {
    observer <- observer[-invalid_points, ]
    x0 <- x0[-invalid_points]
    y0 <- y0[-invalid_points]
    height_0_vec <- height_0_vec[-invalid_points]
  }
  
  if (progress) {
    setTxtProgressBar(pb, 4)
    cat("\n")
  }
  if (length(invalid_points) == 1) {
    message("1 point has been removed, because it was outside of the DSM or DTM")
  } else if (length(invalid_points) > 1) {
    message(paste(length(invalid_points), "points have been removed, because they were outside of the DSM or DTM"))
  }
  
  #### 4. Compute viewshed ####
  # Start row/col
  r0 <- terra::rowFromY(dsm_rast, y0)
  c0 <- terra::colFromX(dsm_rast, x0)
  
  # Apply viewshed (C++) function
  if (progress) {
    cat("\n")
    message(paste0("Computing Distance Analysis for ", nrow(observer), ifelse(nrow(observer)>1, " points:", " point:")))
    cat("\n")
    start_time <- Sys.time()
  }
  distance_tbl <- viewshed_distance_analysis_cpp(dsm_cpp_rast, dsm_vec,
                                                 c0, r0, max_distance, height_0_vec,
                                                 cores, progress)
  colnames(distance_tbl) <- c("V1", "V2", "V3")
  
  distance_tbl <- distance_tbl %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(Visible_perc = V3 / V2) %>%
    dplyr::rename(Distance = V1) %>%
    dplyr::select(Distance, Visible_perc)
  
  
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
  
  #### 5. Compare DSM with Visibility ####
  if (plot) {
    cat("\n")
    message("Building Plot...")
    
    p <- distance_tbl %>% 
      ggplot2::ggplot(ggplot2::aes(x = Distance, y = Visible_perc)) + 
      ggplot2::geom_smooth(formula = y ~ s(x, bs = "cs"), method = "gam") +
      #ggplot2::geom_line() +
      ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      ggplot2::labs(x = "Distance [m]", y = "Visibility [%]") +
      ggplot2::theme_light()    
    print(p)
  }
  
  rm(dsm_vec); invisible(gc())
  return(distance_tbl)
}