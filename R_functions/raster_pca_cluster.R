raster_pca_cluster <- function(input_points, ref_raster, output_crs, output_dir_ras) {
  
  ## This function is built to create raster files out of derived principal componentns and delineated clusters.
  ## First of all, an point data frame is used to transform into a spatial points data frame with RDNew coordinates,
  ## and next, this spatial points DF is transformed into a raster. Next, the rasters are projected to match the same extents
  ## and resolution as a reference raster. In case of datasets representing clusters, the clusters are smoothed by means of
  ## a focal filter of 35x35 and modality function, in order to get rid of isolated pixels and sharp boundary edges.
  
  ## The arguments of this function are:
  ## input_points: the input point csv that is read and pre-processed in this function.
  ## ref_raster: a reference raster with extent and resolution to which the input raster should be projected.
  ## output_crs: the project string representing the output coordinate reference system.
  ## output_dir_ras: the output directory where the processed rasterfiles are stored.
  
  #############################################################################################################################################
  
  # set variable name
  var_name <- names(input_points)[3]
  
  # transform data into a SpatialPointsDataFrame and set correct project string
  formula <- as.formula(paste("~", names(input_points)[1], "+", names(input_points)[2]))
  coordinates(input_points) <- formula
  crs(input_points) <- output_crs
  
  # transform spatial points into a raster and set correct project string again
  var_raster <- rasterFromXYZ(input_points)
  crs(var_raster) <- output_crs
  
  # project the raster to match extents and resolution of other raster datasets
  var_raster <- projectRaster(from = var_raster, to = ref_raster, method = "bilinear")
  
  # in case of clusters, use the focal function to smooth the rasters
  if (TRUE %in% grepl("clus", names(input_points))) {
    smooth_raster <- focal(x = var_raster, w = matrix(1, nc = 35, nr = 35), fun = modal, na.rm = TRUE)
    var_raster <- mask(x = smooth_raster, mask = main_bound_rdnew)
  }
  
  # write the raster to the given output folder
  writeRaster(x = var_raster, filename = paste0(output_dir_ras, "/", var_name, "_raster.tif"),
              format = "GTiff", overwrite = TRUE)
}