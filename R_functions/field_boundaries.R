field_boundaries <- function(input_dir, input_name, output_dir_shp, output_name, output_crs) {
  
  ## This function reads input boundary shapefiles, and then sets the project string or transforms
  ## the file to RDNew coordinates. After that it optionally selects features, and in the end
  ## the transformed boundary shapefiles are stored in the given output folder.
  
  ## The arguments of this function are:
  ## input_dir: the input directory where the point shapefiles are stored.
  ## input_name: the name of an input shapefile that is used in this function.
  ## output_dir_shp: the output directory where the processed shapefiles are stored.
  ## output_name: the name of an output shapefile that is stored in the output directory.
  ## output_crs: the project string representing the output coordinate reference system.
  
  ################################################################################################
  
  # read field or fertilizer boundary
  layer = gsub(".shp", "", input_name)
  boundary <- readOGR(dsn = input_dir, layer = layer)
  
  # set correct project string for field boundary to RDNew (ignore warning)
  # or transform coordinates system for fertilizer boundaries to RDNew
  if (TRUE %in% c(grepl("Main_Field", input_name))) {
    proj4string(boundary) <- CRS(output_crs)
  } else {
    boundary <- spTransform(boundary, CRSobj = CRS(output_crs))
  }
  
  # In case of vertical fertilizer boundaries, select the one with no additional fertilizer
  if (TRUE %in% c(grepl("Vertical", input_name))) {
    boundary <- boundary[boundary@data$OBJECTID == 1, ]
  }
  
  # Store RDNew boundary shapefiles in output folder
  writeOGR(obj = boundary, dsn = output_dir_shp, layer = output_name,
           driver="ESRI Shapefile", overwrite_layer = TRUE)
}