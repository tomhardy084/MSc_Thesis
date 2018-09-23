spatial_vars <- function(input_dir, input_raster_names, boundary) {
  
  ## This function takes raster files as input and masks them according to the field boundary.
  ## Next, it transforms the pixels of each raster file into spatial points dataframe and stores
  ## them as variables including x, y coordinates in a general dataframe, which is returned at
  ## the end of this function. The coordinate reference system of all files are according to RDNew.
  
  ## The arguments of this function are:
  ## input_dir: the input directory where the point shapefiles are stored.
  ## input_raster_names: the names of the input raster files that are used in this function.
  ## boundary: the field boundary in spatialPolygons or spatialPolygonsDataFrame format.
  
  ###########################################################################################################
  
  # initiate a loop to store all variables in one DF including x,y coordinates
  for (i in 1:length(input_raster_names)) {
    
    # set input data source name and load the raster as input
    input_dsn <- paste0(input_dir, "/", input_raster_names[i])
    input_ras <- raster(input_dsn)
    
    # mask raster according to field boundary (to get same extents again for all raster datasets)
    mask_ras <- mask(x = input_ras, mask = boundary)
    
    # change raster cells to spatial points, and points to dataframe
    spatial_points <- rasterToPoints(mask_ras)
    spatial_points <- as.data.frame(spatial_points)
    
    # create and extent dataframe to store variables of all raster datasets
    if (i == 1) {
      variable_df <- spatial_points
    } else {
      variable_df <- cbind(variable_df, spatial_points)
    }
  }
  
  return(variable_df)
}