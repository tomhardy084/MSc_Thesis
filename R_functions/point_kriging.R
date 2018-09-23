point_kriging <- function(input_dir, input_points, output_dir_ras, output_dir_graphs, output_crs,
                          boundary, new_predict, ref_raster) {
  
  ## This function first estimates a variogram from an input variable with its corresponding coordinates,
  ## which is used as input for the kriging function. Based on this variogram, the function constructs
  ## a kriged raster map containing the predictions on new locations for the input variable.
  ## These predictions are determined by a linear model with the input variable as response.
  ## The function returns the fitted variogram parameters as a dataframe, and stores a kriged raster map
  ## in the given output folder. The coordinates and extents are according to RDNew coordinates.
  
  ## The arguments of this function are:
  ## input_dir: the input directory where the point shapefiles are stored.
  ## input_points: the name of an input point file that is used in this function.
  ## output_dir_ras: the output directory where the processed rasterfiles are stored.
  ## output_dir_graphs: the output directory where the histograms and Moran Plots are stored.
  ## output_crs: the project string representing the output coordinate reference system.
  ## boundary: the field boundary in spatialPolygons or spatialPolygonsDataFrame format.
  ## new_predict: a shapefile in spatialPoints or spatialPointsDataFrame with locations of new predictions.
  ## ref_raster: a reference raster with extent and resolution to which the input raster should be projected.
  
  ###########################################################################################################
  
  # read point shapefile from input directory, and set correct project string (ignore warning)
  input_points_name <- gsub(".shp", "", input_points)
  points_rdnew <- readOGR(dsn = input_dir, layer = input_points_name)
  proj4string(points_rdnew) <- CRS(output_crs)
  
  # apply autofitVariogram function to estimate and fit a semivariogram based on linear model
  formula <- as.formula(paste(gsub(" points", "", input_points_name), "~ 1"))
  points_vgm <- autofitVariogram(formula = formula,
                                 input_data = points_rdnew,
                                 model = c("Sph", "Exp"),
                                 miscFitOptions = list(merge.small.bins = TRUE))
  
  # get experimental variogram parameters
  vgm <- points_vgm$var_model
  model <- as.character(vgm[2, "model"])
  nugget <- round(sum(vgm[1, "psill"]), digits = 3)
  sill <- round(sum(vgm[, "psill"]), digits = 3)
  range <- round(sum(vgm[, "range"]), digits = 3)
  
  # store parameters in a dataframe
  variogram_df <- data.frame(model, nugget, sill, range)
  rownames(variogram_df) <- input_points_name
  
  # create plot experimental and fitted variogram
  jpeg(filename = paste0(output_dir_graphs, "/Variogram ", input_points_name, ".jpg"), width = 320, height = 320)
  plot(points_vgm$exp_var, points_vgm$var_model, col = "blue", pch = 1,
       main = paste("Experimental and fitted variogram\n", input_points_name))
  dev.off()
  
  # perform kriging for new prediction locations, based on the autofitted variogram and point locations
  block_size <- as.integer(sqrt(length(ref_raster[!(is.na(ref_raster))]) / length(points_rdnew)))
  kriged_points <- krige(formula = formula, nmin = 3, nmax = 30,
                         block = c(block_size, block_size),
                         locations = points_rdnew,
                         newdata = new_predict,
                         model = points_vgm$var_model)
  
  # transform kriged map into a raster map and set the proper project string for RDNew
  kriged_raster <- rasterFromXYZ(kriged_points)
  crs(kriged_raster) <- prj_string_rdnew
  kriged_raster <- projectRaster(from = kriged_raster, to = ref_raster, method = "bilinear")
  
  # write the raster as a new TIFF file to the output folder
  if (require(rgdal)) {
    kriged_raster_dsn <- paste0(output_dir_ras, "/", gsub(" points", "", input_points_name), ".tif")
    writeRaster(x = kriged_raster, filename = kriged_raster_dsn, format = "GTiff", overwrite = TRUE)
  }
  
  return(variogram_df)
}