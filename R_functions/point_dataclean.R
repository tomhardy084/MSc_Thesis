point_dataclean <- function(input_points, output_crs, output_dir_shp, output_dir_graphs, eucl_dist = 0) {
  
  ## This function is created for pre-processing input point datasets. Pre-processing consists of the detection and
  ## removal of outliers or inliers. Outliers are based on values lower or higher than three times the standard deviation,
  ## and inliers are based on local Moran's Index. The function stores pre-processed point shapefiles, boxplots and
  ## Moran Plots in the given output folders. The shapefile coordinates and extents are according to RDNew coordinates.

  ## The arguments of this function are:
  ## input_points: the input point csv that is read and pre-processed in this function.
  ## output_crs: the project string representing the output coordinate reference system.
  ## output_dir_shp: the output directory where the processed shapefiles are stored.
  ## output_dir_graphs: the output directory where the histograms and Moran Plots are stored.
  ## eucl_dist (only for inlier detection): the Euclidean distance [m] used to identify the number of neighbouring points.
  
  #########################################################################################################################
  
  # select values from input dataset that are not NA
  input_points_select <- which(!is.na(input_points[, 3]))
  input_points <- input_points[input_points_select, ]
  
  # create boxplot of point variable and store it in output folder
  jpeg(filename = paste0(output_dir_graphs, "/Boxplot ", names(input_points)[3], ".jpg"), width = 320, height = 320)
  boxplot(input_points[, 3], col = "gray", main = paste("Boxplot", names(input_points)[3]), ylab = names(input_points)[3])
  dev.off()
  
  if (TRUE %in% c(grepl("elevation", names(input_points)), grepl("soil", names(input_points)))) {
    # outlier detection and removal for elevation, soil depth and soil OM
    # calculate mean, standard deviation, and lower and upper bounds
    points_mean <- mean(input_points[, 3])
    points_std <- sd(input_points[, 3])
    points_lower <- points_mean - 3 * points_std
    points_upper <- points_mean + 3 * points_std
    
    # remove outliers based on lower and upper bounds
    points_select <- which(input_points[, 3] >= points_lower & input_points[, 3] <= points_upper)
    points_reduced <- input_points[points_select, ]
    
  } else {
    # inlier detection and removal for soil ec and crop yield
    # identify neighbouring points by Eucledian distance, calculate spatial weights and store them in a list
    coords <- coordinates(input_points[, 1:2])
    neighbours <- dnearneigh(coords, 0, eucl_dist)
    neighbours_list <- nb2listw(neighbours, style = "S", zero.policy = TRUE)
    
    # create Moran Plot
    jpeg(filename = paste0(output_dir_graphs, "/Moran Plot ", names(input_points)[3], ".jpg"), width = 320, height = 320)
    moran_plot <- moran.plot(x = input_points[, 3], listw = neighbours_list, quiet = TRUE, labels = FALSE,
                             zero.policy = TRUE, main = "Moran Plot with potential inliers", col = "cyan", pch = 20,
                             xlab = names(input_points)[3], ylab = paste(names(input_points)[3], "spatially lagged"))
    dev.off()
    
    # define Local Moran Index and add columns with this information to points DF
    moran_index <- localmoran(x = input_points[, 3], listw = neighbours_list, zero.policy = TRUE,
                              p.adjust.method = "none", alternative ="less")
    input_points <- data.frame(input_points, moran_index)
    
    # selection of cases that have a positive Local Moran Index, are statistically non-significant (p > 0.05)
    # and have no inliers anymore identified with Moran Plot
    points_select <- which(input_points$Ii > 0 & input_points$Pr.z...0. > 0.05)
    points_reduced <- input_points[points_select, 1:3]
  }
  
  # transform data into a SpatialPointsDataFrame again
  formula <- as.formula(paste("~", names(points_reduced)[1], "+", names(points_reduced)[2]))
  coordinates(points_reduced) <- formula
  proj4string(points_reduced) <- CRS(output_crs)
  
  # store the point shapefile in the output folder
  writeOGR(obj = points_reduced, dsn = output_dir_shp,
           layer = paste(names(input_points)[3], "points"),
           driver = "ESRI Shapefile", overwrite_layer = TRUE)
}