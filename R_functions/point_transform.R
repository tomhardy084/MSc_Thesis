point_transform <- function(input_dir, input_crs, output_crs, coords, vars, var_names, boundary) {
  
  ## This function is created for transforming the reference system of input point data to RDNew coordinates.
  ## The function returns dataframes with coordinates according to this reference system and their corresponding variables.

  ## The arguments of this function are:
  ## input_dir: the input directory where a csv file with spatial point data is stored.
  ## input_crs: the project string representing the input reference system.
  ## output_crs: the project string representing the output reference system.
  ## coords: vector with two integers indicating the columns that represent the spatial coordinates.
  ## vars: vector with number of integers indicating the colums that represent the variable(s).
  ## var_names: vector with the names of the coordinates and variable(s).
  ## boundary: the field boundary in spatialPolygons or spatialPolygonsDataFrame format.
  
  ########################################################################################################################
  
  # read csv file, and select coordinates and variables
  point_file <- list.files(path = input_dir, pattern = glob2rx("*.csv"))
  points_csv <- read.csv(file = paste0(input_dir, point_file), header = TRUE, sep = ",", dec = ".")
  points_df <- points_csv[, c(coords, vars)]
  
  # if necessary, rearange coordinates to match the sequence (x, y) or (lon, lat)
  if (points_df[1, 1] > points_df[1, 2]) {
    points_df[, 1:2] <- points_df[, 2:1]
  }
  
  # assign variable names to dataframe
  names(points_df) <- var_names
  
  # transform data into a SpatialPointsDataFrame, and transform coordinates to RDNew
  formula <- as.formula(paste("~", names(points_df)[1], "+", names(points_df)[2]))
  coordinates(points_df) <- formula
  proj4string(points_df) <- CRS(input_crs)
  points_rdnew <- spTransform(points_df, CRS(output_crs))
  
  # select points inside the field boundary
  points_rdnew <- points_rdnew[boundary, ]
  
  # transform data back to dataframe and rearrange columns
  points_rdnew <- as.data.frame(points_rdnew)
  points_rdnew <- points_rdnew[, c(length(var_names) - 1, length(var_names), 1:(length(var_names) - 2))]
  
  return(points_rdnew)
}