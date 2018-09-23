####################################
###  Project: MGI Thesis         ###
###  Title: data pre-processing  ###
###  Author: Tom Hardy           ###
###  Date: August 30, 2018       ###
####################################


### INITIALIZATION -----------------------------------------------------------------------------------------

# empty work environment
rm(list = ls())

# read names of directories from file and create the directories
directories <- read.delim("directories.txt", header = FALSE)
directories <- as.character(directories[, 1])
lapply(directories, FUN = function(X) {if (!(file.exists(X)))
{dir.create(file.path(X))}})

# read names of packages from file and install and load the packages
packages <- read.delim("packages.txt", header = FALSE)
packages <- as.character(packages[, 1])
lapply(packages, FUN = function(X) {if(!(X %in% rownames(installed.packages())))
  {install.packages(X)}})
lapply(packages, require, character.only = TRUE)

# read names of CRS objects from file and define project strings from those objects
crs_obj <- read.delim("crs objects.txt", header = FALSE)
crs_obj <- as.character(crs_obj[, 1])
prj_string_wgs84 <- crs_obj[1]
prj_string_utm31 <- crs_obj[2]
prj_string_rdnew <- crs_obj[3]

# set output directories to store shapefiles, rasterfiles, graphs and tables
output_dir_shp <- "R_output/R_shapefiles"
output_dir_ras <- "R_output/R_rasterfiles"
output_dir_graphs <- "R_output/R_graphs_tables"

# store function names in a list and source the functions
functions <- list.files("R_functions")
lapply(functions, FUN = function(X) {source(paste0("R_functions/", X))})



### 1. FIELD AND FERTILIZATION BOUNDARIES ------------------------------------------------------------------

# set input directory, and input and output names for all boundary shapefiles
input_dir <- "Data/Field boundaries/field and fert boundaries"
input_shp <- list.files(path = input_dir, pattern = glob2rx("*.shp"))
output_shp <- c("fert_bound_grid_rdnew", "fert_bound_horizontal_rdnew",
                "fert_bound_vertical_rdnew", "main_boundary_rdnew")

# create loop to call boundary processing function to transform all field and
# fertilizer boundaries to RDNew coordinates, and store them in output folder.
for (i in 1:length(input_shp)) {
  field_boundaries(input_dir = input_dir, input_name = input_shp[i],
                   output_dir_shp = output_dir_shp, output_name = output_shp[i],
                   output_crs = prj_string_rdnew)
}

# load main field boundary in R environment and set correct project string field boundary (ignore warning)
main_bound_rdnew <- readOGR(dsn = output_dir_shp, layer = "main_boundary_rdnew")
proj4string(main_bound_rdnew) <- CRS(prj_string_rdnew)



### 2. SOIL AND CROP IMAGES FROM UAVs ----------------------------------------------------------------------

## 2.1 UAV images RGB bare soil April 2015 -----------------------------------------------------------------

# set input and output directories and file names and store output raster names in a list
input_dir <- "Data/eBee images/20150414_Tegenover_de_Schuur_Orthomosaic_RGB/"
input_ras <- list.files(path = input_dir, pattern = glob2rx("*.tif"))
output_ras <- c("bare_soil_red.tif", "bare_soil_green.tif", "bare_soil_blue.tif")

# create loop to set input and output dsn and call uav_processing function
# to aggregate, transform and mask bare soil RGB images from April 2015
for (i in 1:length(output_ras)) {
  input_dsn <- paste0(input_dir, input_ras)
  output_dsn <- paste0(output_dir_ras, "/", output_ras[i])
  image_processing(input_dsn = input_dsn, output_dsn = output_dsn,
                   band_num = i, boundary = main_bound_rdnew,
                   crs_rdnew = prj_string_rdnew,
                   crs_utm31 = prj_string_utm31)
}

# load the pre-processed RGB raster files
bare_soil_red <- raster(paste0(output_dir_ras, "/bare_soil_red.tif"))
bare_soil_green <- raster(paste0(output_dir_ras, "/bare_soil_green.tif"))
bare_soil_blue <- raster(paste0(output_dir_ras, "/bare_soil_blue.tif"))

# calculate NDRG and SUMVIS indices and create True Color image
bare_soil_ndrg <- (bare_soil_red - bare_soil_green) / (bare_soil_red + bare_soil_green)
bare_soil_sumvis <- bare_soil_red + bare_soil_green + bare_soil_blue
bare_soil_true_color <- brick(bare_soil_red, bare_soil_green, bare_soil_blue)

# store the NDRG, SUMVIS and True Color images in output folder
writeRaster(x = bare_soil_ndrg, filename = paste0(output_dir_ras, "/bare_soil_ndrg.tif"),
            format = "GTiff", overwrite = TRUE)
writeRaster(x = bare_soil_sumvis, filename = paste0(output_dir_ras, "/bare_soil_sumvis.tif"),
            format = "GTiff", overwrite = TRUE)
writeRaster(x = bare_soil_true_color, filename = paste0(output_dir_ras, "/bare_soil_true_color.tif"),
            format = "GTiff", overwrite = TRUE)


## 2.2 UAV images NDVI and WDVI August 2015 ----------------------------------------------------------------

# set input and output directories and file names
input_dir <- "Data/eBee images/15_08_21__pix4d/"
input_ras <- list.files(path = input_dir, pattern = glob2rx("*.tif"))
output_ras <- gsub("15_08_21__pix4d_", "crop", input_ras)

# create loop to set input and output dsn and call uav_processing function
# to aggregate, transform and mask crop images from August 2015
for (i in 1:length(output_ras)) {
  input_dsn <- paste0(input_dir, input_ras[i])
  output_dsn <- paste0(output_dir_ras, "/", output_ras[i])
  image_processing(input_dsn = input_dsn, output_dsn = output_dsn,
                   boundary = main_bound_rdnew,
                   crs_rdnew = prj_string_rdnew,
                   crs_utm31 = prj_string_utm31,
                   ref_raster = bare_soil_red)
}

# load the pre-processed NIR, Red, Green raster files
crop_NIR <- raster(paste0(output_dir_ras, "/crop_NIR.tif"))
crop_RED <- raster(paste0(output_dir_ras, "/crop_RED.tif"))
crop_GREEN <- raster(paste0(output_dir_ras, "/crop_GREEN.tif"))

# create False Color image
crop_false_color <- brick(crop_NIR, crop_RED, crop_GREEN)

# store the False Color image in output folder
writeRaster(x = crop_false_color, filename = paste0(output_dir_ras, "/crop_false_color.tif"),
            format = "GTiff", overwrite = TRUE)



### 3. SOIL AND CROP MEASUREMENTS PROXIMAL SENSING AND SAMPLING --------------------------------------------

## 3.1 Point data processing soil electric conductivity ----------------------------------------------------

# set input and directory and files name, read csv file, and select data
input_dir <- "Data/Electric conductivity/"
coords <- 1:2
vars <- 4:7
var_names <- c("x", "y", "ECa_0_5m", "ECa_1_0m", "ECa_1_5m", "ECa_3_0m")

# call point transformation function for transforming EC datasets to RDNew coordinates
points_rdnew <- point_transform(input_dir = input_dir, input_crs = prj_string_wgs84,
                                output_crs = prj_string_rdnew, coords = coords, vars = vars,
                                var_names = var_names, boundary = main_bound_rdnew)

# create loop to call point data cleaning function for outlier and inlier detection and removal
for (i in 3:length(var_names)) {
  point_dataclean(input_points = points_rdnew[c(1:2, i)], output_crs = prj_string_rdnew, eucl_dist = 2.5,
                  output_dir_shp = output_dir_shp, output_dir_graphs = output_dir_graphs)
}


## 3.2 Point data processing elevation and A-horizon depth -------------------------------------------------

# set input and directory and files name, read csv file, and select data
input_dir <- "Data/Elevation and soil depth/"
coords <- 2:3
vars <- 4:5
var_names <- c("x", "y", "elevation", "soil_depth")

# call point transformation function for transforming elevation and soil depth datasets to RDNew coordinates
points_rdnew <- point_transform(input_dir = input_dir, input_crs = prj_string_rdnew,
                                output_crs = prj_string_rdnew, coords = coords, vars = vars,
                                var_names = var_names, boundary = main_bound_rdnew)

# create loop to call point data cleaning function for outlier and inlier detection and removal
for (i in 3:length(var_names)) {
  point_dataclean(input_points = points_rdnew[c(1:2, i)], output_crs = prj_string_rdnew,
                  output_dir_shp = output_dir_shp, output_dir_graphs = output_dir_graphs)
}


## 3.3 Point data processing soil organic matter -----------------------------------------------------------

# set input and directory and files name, read csv file, and select data
input_dir <- "Data/Organic matter/"
coords <- 4:5
vars <- 3
var_names <- c("x", "y", "soil_OM")

# call point transformation function for transforming soil OM dataset to RDNew coordinates
points_rdnew <- point_transform(input_dir = input_dir, input_crs = prj_string_wgs84,
                                output_crs = prj_string_rdnew, coords = coords, vars = vars,
                                var_names = var_names, boundary = main_bound_rdnew)

# call point data cleaning function for outlier and inlier detection and removal
point_dataclean(input_points = points_rdnew, output_crs = prj_string_rdnew,
                output_dir_shp = output_dir_shp, output_dir_graphs = output_dir_graphs)


## 3.4 Point data processing crop yield --------------------------------------------------------------------

# set input and directory and files name, read csv file, and select data
input_dir <- "Data/Crop yield/"
coords <- 1:2
vars <- 20
var_names <- c("x", "y", "crop_yield")

# call point transformation function for transforming crop yield dataset to RDNew coordinates
points_rdnew <- point_transform(input_dir = input_dir, input_crs = prj_string_wgs84,
                                output_crs = prj_string_rdnew, coords = coords, vars = vars,
                                var_names = var_names, boundary = main_bound_rdnew)

# call point data cleaning function for outlier and inlier detection and removal
point_dataclean(input_points = points_rdnew, output_crs = prj_string_rdnew, eucl_dist = 2,
                output_dir_shp = output_dir_shp, output_dir_graphs = output_dir_graphs)


## 3.5 Interpolation by kriging all point datasets ---------------------------------------------------------

# set input data source name and file names
input_dir <- "R_output/R_shapefiles"
input_shp <- list.files(path = input_dir, pattern = glob2rx("*points.shp"))

# use reference raster as input for prediction locations and change cells to points
krig_raster <- bare_soil_red
krig_points <- rasterToPoints(krig_raster)

# change points to dataframe, and dataframe to spatialpoints dataframe
krig_points <- as.data.frame(krig_points)
coordinates(krig_points) <- ~ x + y
crs(krig_points) <- prj_string_rdnew

# create loop to call point kriging function for spatial interpolation of all point datasets
for (i in 1:length(input_shp)) {
  point_krig <- point_kriging(input_dir = input_dir, input_points = input_shp[i],
                              output_dir_ras = output_dir_ras, output_dir_graphs = output_dir_graphs,
                              output_crs = prj_string_rdnew, boundary = main_bound_rdnew,
                              new_predict = krig_points, ref_raster = bare_soil_red)
  
  # create and extent dataframe to store variogram parameters of all point datasets
  if (i == 1) {
    variogram_parameters <- point_krig
  } else {
    variogram_parameters <- rbind(variogram_parameters, point_krig)
  }
}

# write variogram parameters table in given output folder
write.table(variogram_parameters, file = paste0(output_dir_graphs, "/variogram parameters.txt"),
            col.names = NA, quote = FALSE, sep = ",")


## 3.6 calculate total soil organic matter -----------------------------------------------------------------

# load the kriged soil organic matter and soil depth raster files
soil_om <- raster(paste0(output_dir_ras, "/soil_OM.tif"))
soil_depth <- raster(paste0(output_dir_ras, "/soil_depth.tif"))

# calculate total soil organic matter and store the raster file in output folder
soil_om_total <- soil_om / 100 * soil_depth * 1.13 * 10000
writeRaster(x = soil_om_total, filename = paste0(output_dir_ras, "/soil_om_total.tif"),
            format = "GTiff", overwrite = TRUE)


#### end ###################################################################################################