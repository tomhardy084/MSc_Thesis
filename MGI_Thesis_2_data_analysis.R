####################################
###  Project: MGI Thesis         ###
###  Title: data analysis        ###
###  Author: Tom Hardy           ###
###  Date: August 30, 2018       ###
####################################


### INITIALIZATION -----------------------------------------------------------------------------------------

# empty work environment
rm(list = ls())

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

# set output directories to store shapefiles, rasterfiles, graphs, tables and scatterplots
output_dir_shp <- "R_output/R_shapefiles"
output_dir_ras <- "R_output/R_rasterfiles"
output_dir_graphs <- "R_output/R_graphs_tables"
output_dir_scatter <- "R_output/R_scatterplots"

# store function names in a list and source the functions
functions <- list.files("R_functions")
lapply(functions, FUN = function(X) {source(paste0("R_functions/", X))})

# load main field boundary, set correct project string field boundary (ignore warning)
# and load bare_soil_red as reference raster
main_bound_rdnew <- readOGR(dsn = output_dir_shp, layer = "main_boundary_rdnew")
proj4string(main_bound_rdnew) <- CRS(prj_string_rdnew)
bare_soil_red <- raster(paste0(output_dir_ras, "/bare_soil_red.tif"))



### 1. CREATE DATAFRAME TO STORE ALL VARIABLES INCLUDING COORDINATES ---------------------------------------

# store input raster file names in a list
input_ras_names <- list.files(path = output_dir_ras, pattern = glob2rx("*.tif"))
input_ras_names <- input_ras_names[which(!(grepl("color", input_ras_names)))]

# call spatial variables function to store all raster datasets as variables into one dataframe
variable_df <- spatial_vars(input_dir = output_dir_ras,
                            input_raster_names = input_ras_names,
                            boundary = main_bound_rdnew)
  
# select and rearrange columns of variable data frame
variable_df <- variable_df[c("x", "y", "bare_soil_red", "bare_soil_green", "bare_soil_blue", "bare_soil_ndrg",
                             "bare_soil_sumvis", "ECa_0_5m", "ECa_1_0m", "ECa_1_5m", "ECa_3_0m", "elevation",
                             "soil_depth", "soil_OM", "soil_om_total", "crop_yield", "crop_NDVI", "crop_WDVI")]

# write data frame as csv file in output folder (to save memory in R)
write.csv(x = variable_df, file = paste0(output_dir_graphs, "/Soil and crop variables.csv"),
            col.names = NA, quote = FALSE, sep = ",")



### 2. DESCRIPTIVE STATISTICS ------------------------------------------------------------------------------

# read variables data frame into R memory
variable_df <- read.csv(file = paste0(output_dir_graphs, "/Soil and crop variables.csv"), header = TRUE, sep = ",", dec = ".")[, -1]

## 2.1 Descriptive statistics dataframe --------------------------------------------------------------------

# call desc_stats function to determine descriptive statistics of
# input variables and store them as a table in the given output folder
desc_stats(variables_dataframe = variable_df,
           output_dir_graphs = output_dir_graphs)


## 2.2 Create scatterplots and correlation matrix ----------------------------------------------------------

# call scatter_plots function to generate density scatterplots for each
# pair of input variables and store them in the given output folder
scatter_plots(variables_dataframe = variable_df,
              output_dir_scatter = output_dir_scatter)

# create correlation matrix and write matrix in output folder
cor_mat <- cor(variable_df[, 3:ncol(variable_df)], use = "complete.obs")
cor_mat_round <- round(cor_mat, digits = 3)
write.table(cor_mat_round, file = paste0(output_dir_graphs, "/Correlation matrix.txt"),
            col.names = NA, quote = FALSE, sep = ",")

# create correlation matrix plot and store in output folder
col <- colorRampPalette(c("darkred", "yellow", "darkgreen"))(20)
jpeg(filename = paste0(output_dir_graphs, "/Correlation matrix.jpg"), width = 750, height = 750)
corrplot(corr = cor_mat_round, method = "color", type = "upper", col = col, tl.col = "black",
         tl.srt = 45, addgrid.col = "gray", addCoef.col = "black", number.cex = 0.8,
         diag = FALSE, mar = c(0, 0, 0, 0))
dev.off()



### 3. PRINCIPAL COMPONENT ANALYSIS AND CLUSTER ANALYSIS ---------------------------------------------------

## 3.1 Spatial principal component analysis (MULTISPATI-PCA) -----------------------------------------------

# call bartlett_kmo_vars function to perform Bartlett's test of Sphericity and
# KMO MSA tests, and select variables based on KMO values larger than 0.5
pca_variables <- bartlett_kmo_vars(cor_mat = cor_mat[1:13, 1:13], variable_df = variable_df,
                                   output_dir_graphs = output_dir_graphs)

# call multi_pca function to perform PCA and create a screeplot to determine the number of PCs
multi_pca(pca_variables = pca_variables, output_dir_graphs = output_dir_graphs, screeplot = TRUE)

# call multi_pca function to perform PCA, MPCA and store output in dataframes and output folder
pca_variables_with_pcs <- multi_pca(pca_variables = pca_variables, output_dir_graphs = output_dir_graphs,
                                    num_PC = 4, eucl_dist = 10)


## 3.2 Cluster Analysis ------------------------------------------------------------------------------------

# call cluster_analysis function to perform cluster analysis and create a screeplot
cluster_analysis(cluster_variables = pca_variables_with_pcs,
                 output_dir_graphs = output_dir_graphs, screeplot = TRUE)

# call cluster_analysis function to perform cluster analysis and store clusters into a data frame
cluster_analysis(cluster_variables = pca_variables_with_pcs,
                 output_dir_graphs = output_dir_graphs, num_clus = 2:4)


## 3.3 Store raster files of PCs and clusters in output folder ---------------------------------------------

# read variables data frame into R memory
pca_variables_with_clusters <- read.csv(file = paste0(output_dir_graphs, "/pca variables with clusters.csv"),
                                        header = TRUE, sep = ",", dec = ".")[, -1]

# create loop to call function for transforming PCs and clusters into rasters
# and storing them in the output folder
for (i in 11:length(pca_variables_with_clusters)) {
  raster_pca_cluster(input_points = pca_variables_with_clusters[c(1:2, i)],
                     ref_raster = bare_soil_red,
                     output_dir_ras = output_dir_ras,
                     output_crs = prj_string_rdnew)
}



# 4. VALIDATION OF CLUSTERS (MANAGEMENT ZONES) -------------------------------------------------------------

## 4.1 Create an overlay between crop yield points and clusters -------------------------------------------

# read crop yield points and set correct project string for RDNew (ignore warning)
crop_yield_points <- readOGR(dsn = output_dir_shp, layer = "crop_yield points")
proj4string(crop_yield_points) <- CRS(prj_string_rdnew)

# read raster with four clusters and transform it into a spatial grid data frame
cluster_raster <- raster(paste0(output_dir_ras, "/kmeans_4_clus_raster.tif"))
cluster_grid <- as(cluster_raster, "SpatialGridDataFrame")

# create an overlay to attach the clusters to the crop yield points
cluster_points <- sp::over(crop_yield_points, cluster_grid)
crop_yield_points@data$clusters <- as.numeric(cluster_points[, 1])

# store the new point shapefile in the output folder
writeOGR(obj = crop_yield_points, dsn = output_dir_shp,
         layer = "crop_yield points with clusters",
         driver = "ESRI Shapefile", overwrite_layer = TRUE)


## 4.2 Fit linear models for validation of clusters  -------------------------------------------------------

# read crop yield points + clusters and set correct project string for RDNew (ignore warning)
crop_yield_clusters <- readOGR(dsn = output_dir_shp, layer = "crop_yield points with clusters")
proj4string(crop_yield_clusters) <- CRS(prj_string_rdnew)

# read fertilizer boundaries and set correct project string for RDNew (ignore warning)
fertilizer_polygons <- readOGR(dsn = output_dir_shp, layer = "fert_bound_grid_rdnew")
proj4string(fertilizer_polygons) <- CRS(prj_string_rdnew)

# initiate polygon list for extracting sub-polygons from the fertilizer boundaries
poly_list <- list(c(as.character(1:4), as.character(9:16)), "1", "3")

# initiate loop to call linear_model function for validating the selected management zones
for (i in 1:length(poly_list)) {
  linear_model(crop_yield = crop_yield_clusters,
               fertilizer_polygons = fertilizer_polygons,
               polygon_id = poly_list[[i]], iteration = i,
               output_dir_graphs = output_dir_graphs)
}


#### end ###################################################################################################