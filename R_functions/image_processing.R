image_processing <- function(input_dsn, output_dsn, band_num = 1, boundary, crs_rdnew, crs_utm31, ref_raster = 0) {
  
  ## This function is created for pre-processing input UAV images. Pre-processing consists of cropping the input image,
  ## aggregating the cropped image to lower the resolution, transform the raster to rdnew coordinates, and lastly
  ## to mask the image according to the field boundary. The function creates a masked uav image with a spatial resolution
  ## of 1m and stores it in the given output folder. The coordinates and extents are according to RDNew coordinates.
  
  ## The arguments of this function are:
  ## input_dsn: the data source name (directory + file name) of the input raster file.
  ## output_dsn: the data source name (directory + file name) of the pre-processed raster file.
  ## band_num: the raster's band number in the geotiff file (if only one band is present, band_num = 1).
  ## boundary: the field boundary in spatialPolygons or spatialPolygonsDataFrame format.
  ## crs_rdnew: the project string representing RDNew coordinates.
  ## crs_utm31: the project string representing UTM31 coordinates.
  ## ref_raster: a reference raster with extent and resolution to which the input raster should be projected.
  
  #######################################################################################################################
  
  # define coordinates (in utm31 format) to create a rectangular extent for cropping the input raster
  coords = matrix(c(651196, 5687021, 651196, 5687656, 651662, 5687656, 651662, 5687021),
                  ncol = 2, byrow = TRUE)
  crop_poly <- Polygon(coords)
  crop_extent <- SpatialPolygons(list(Polygons(list(crop_poly), ID = "1")), proj4string = CRS(crs_utm31))
  
  # load the raster including band number and crop it according to the rectangular extent
  uav_raster <- raster(input_dsn, band_num)
  uav_crop <- crop(x = uav_raster, y = crop_extent)
  
  # aggregate the cropped image, project the new raster and mask it according to the field boundary
  # in case of a vegetation index, project it a second time according to the reference raster to match its parameters
  uav_agg <- aggregate(uav_crop, fact = 2/res(uav_crop)[1], fun = mean, expand = TRUE)
  uav_rdnew <- projectRaster(from = uav_agg, res = c(2, 2), crs = crs_rdnew, method = "bilinear")
  uav_image <- mask(x = uav_rdnew, mask = boundary)
  if (TRUE %in% grepl("15_08_21__pix4d", input_dsn)) {
    uav_image <- projectRaster(from = uav_image, to = ref_raster, method = "bilinear")
  }
  
  # write the raster as a new TIFF file to the output folder
  if (require(rgdal)) {
    writeRaster(x = uav_image, filename = output_dsn, format = "GTiff", overwrite = TRUE)
  }
}