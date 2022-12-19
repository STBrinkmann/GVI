test_that("VVI calculation works", {
  library(terra)
  library(raster)
  library(sf)
  # Download DEM
  DEM_tmp <- tempfile(fileext = ".tif")
  download.file(url = "https://github.com/STBrinkmann/data/raw/main/GVI_Data/GVI_DEM.tif",
                destfile = DEM_tmp, mode="wb")
  
  # Download DSM
  DSM_tmp <- tempfile(fileext = ".tif")
  download.file(url = "https://github.com/STBrinkmann/data/raw/main/GVI_Data/GVI_DSM.tif",
                destfile = DSM_tmp, mode="wb")
  
  DEM <- terra::rast(DEM_tmp)
  DSM <- terra::rast(DSM_tmp)
  
  observer <- sf::st_sfc(
    sf::st_point(c(492243.3, 5454231.4)),
    crs = sf::st_crs(26910)) 
  observer <- sf::st_as_sf(observer)
  
  vvi <- vvi_from_sf(
    observer = observer,
    dsm_rast = DSM,
    dtm_rast = DEM,
    max_distance = 200,
    observer_height = 1.7,
    raster_res = NULL,
    cores = 1)
  
  vp <- visibility_proportion(
    observer = observer,
    dsm_rast = DSM,
    dtm_rast = DEM,
    max_distance = 200,
    observer_height = 1.7,
    raster_res = NULL
  )

  expect_equal(
    round(vvi$VVI, 3),
    round(vp, 3))
  
  vvi_poly25 <- vvi_from_sf(
    observer = observer %>% st_buffer(25),
    dsm_rast = DSM,
    dtm_rast = DEM,
    max_distance = 200,
    observer_height = 1.7,
    raster_res = NULL,
    cores = 1)
  
  expect_s3_class(
    vvi_poly25,
    "sf")

})
