linear_model <- function(crop_yield, fertilizer_polygons, polygon_id, iteration = 1, output_dir_graphs) {

  ## This function is built for validating the derived management zones for validation areas that are based on fertilizer zones.
  ## First of all, one validation area is extracted from a fertilizer zone input dataset, and input points of crop yield within this boundary
  ## are selected. Next, random samples of 100 observations per cluster within this area are drawn, and five linear models are fitted
  ## on these samples. The model with the lowest Akaike Information Criterion (AIC) value is selected for further analysis. The results of
  ## this model are extracted and stored in the given output folder.
  
  ## The arguments of this function are:
  ## crop_yield: a point shapefile containing crop yield measurements and the cluster in which the points are located.
  ## fertilizer_polygons: a polygons shapefile containing fertilizer boundaries subdivided in 16 grids.
  ## polygon_id: the IDs of polygons which to select from the fertilizer boundaries.
  ## iteration: if the function is called in a loop, the current iteration of the loop. If no loop is present, iteration = 1.
  ## output_dir_graphs: the output directory where the linear model outputs are stored.
  
  #############################################################################################################################################
  
  # extract polygon and select crop yield points within this polygon
  select_poly <- fertilizer_polygons[fertilizer_polygons$id %in% polygon_id, ]
  crop_yield_select <- crop_yield[select_poly, ]
  
  # transform crop yield point SPDF into regular DF and turn clusters into a factor
  crop_yield_df <- as.data.frame(crop_yield_select)
  
  # initiate loop to sample crop yield points (100 per available cluster)
  set.seed(5)
  for (j in 1:length(unique(na.omit(crop_yield_df$clusters)))) {
    
    # select cluster
    select_cluster <- which(crop_yield_df$clusters == unique(na.omit(crop_yield_df$clusters))[j])
    
    # draw a random sample of 100 crop yield observations per cluster
    if (length(select_cluster) >= 100) {
      crop_yield_cluster <- crop_yield_df[select_cluster, ]
      sample_size <- sort(sample(x = nrow(crop_yield_cluster), size = 100, replace = FALSE))
      sample_per_cluster <- crop_yield_cluster[sample_size, ]
      
      # bind the samples for all clusters together into one new data frame
      if (!(exists(x = "crop_yield_sample"))) {
        crop_yield_sample <- sample_per_cluster
      } else {
        crop_yield_sample <- rbind(crop_yield_sample, sample_per_cluster)
      }
    }
  }
  
  # transform cluster variable into a factor
  crop_yield_sample$clusters <- as.factor(crop_yield_sample$clusters)
  
  # fit linear models based on generalized least squares
  model_fit_1 <- gls(model = crop_yield~clusters, data = crop_yield_sample, method = "REML", na.action = na.omit,
                     correlation = corExp(form = ~coords.x1+coords.x2, metric = "euclidean", nugget = TRUE))
  
  model_fit_2 <- gls(model = crop_yield~clusters, data = crop_yield_sample, method = "REML", na.action = na.omit,
                     correlation = corExp(form = ~coords.x1+coords.x2, metric = "euclidean", nugget = FALSE))
  
  model_fit_3 <- gls(model = crop_yield~clusters, data = crop_yield_sample, method = "REML", na.action = na.omit,
                     correlation = corSpher(form = ~coords.x1+coords.x2, metric = "euclidean", nugget = TRUE))
  
  model_fit_4 <- gls(model = crop_yield~clusters, data = crop_yield_sample, method = "REML", na.action = na.omit,
                     correlation = corSpher(form = ~coords.x1+coords.x2, metric = "euclidean", nugget = FALSE))
  
  model_fit_5 <- gls(model = crop_yield~clusters, data = crop_yield_sample, method = "REML", na.action = na.omit)
  
  # store the five models in a list
  model_list <- list(model_fit_1, model_fit_2, model_fit_3, model_fit_4, model_fit_5)
  
  # calculate Akaike Information Criterion (AIC) for each of the models
  AIC <- vector()
  for (val in 1:length(model_list)) {
    AIC[val] <- AIC(model_list[[val]])
  }
  
  # store AIC values in a data frame
  AIC_df <- round(data.frame(AIC), digits = 1)
  rownames(AIC_df) <- c("exponential correlation with nugget",
                        "exponential correlation without nugget",
                        "spherical correlation with nugget",
                        "spherical correlation without nugget",
                        "independent errors")
  
  # select model with lowest AIC
  select_model <- model_list[[which(AIC == min(AIC))]]
  
  # extract ANOVA, LSMEANS and pairwise comparisons from model
  model_anova <- anova(select_model)
  model_lsmeans <- summary(lsmeans(object = select_model, specs = pairwise ~ clusters,
                                   data = crop_yield_sample[, 1:2], adjust = "tukey"))
  model_means <- model_lsmeans$lsmeans
  model_pairs <- model_lsmeans$contrasts
  
  # store the results in output folder
  write.table(AIC_df, file = paste0(output_dir_graphs, "/Val area ", iteration, " AIC.txt"),
              col.names = NA, quote = FALSE, sep = ",")
  
  write.table(model_anova, file = paste0(output_dir_graphs, "/Val area ", iteration, " ANOVA.txt"),
              col.names = NA, quote = FALSE, sep = ",")
  
  write.table(model_means, file = paste0(output_dir_graphs, "/Val area ", iteration, " LS means.txt"),
              col.names = NA, quote = FALSE, sep = ",")
  
  write.table(model_pairs, file = paste0(output_dir_graphs, "/Val area ", iteration, " Pairs.txt"),
              col.names = NA, quote = FALSE, sep = ",")
  
  # plot side-by-side boxplot of crop yield per MZ and store in output folder
  jpeg(filename = paste0(output_dir_graphs, "/Val area ", iteration, " Boxplots.jpg"), width = 290, height = 290)
  boxplot(crop_yield_df$crop_yield ~ crop_yield_df$clusters, main = paste("Crop yield per MZ, val. area", iteration),
          xlab = "management zone", ylab = "potato crop yield (ton/ha)", col = "gray", outline = FALSE)
  dev.off()
  
  # plot barplot of crop yield means per MZ with error bars (based on confidence intervals) and store in output folder
  jpeg(filename = paste0(output_dir_graphs,"/Val area ", iteration, " LS means.jpg"), width = 290, height = 290)
  barplot(height = model_means[, 2], width = 1, names.arg = model_means$clusters, col = "gray", main = paste("Crop yield mean per MZ, val. area", iteration),
          xlab = "management zone", ylab = "least square means crop yield", xlim = c(0, nrow(model_means) + 0.5), ylim = c(0, 100))
  arrows(x0 = 1:nrow(model_means) * 1.2 - 0.5, y0 = model_means[, 2] - model_means[, 3],
         x1 = 1:nrow(model_means) * 1.2 - 0.5, y1 = model_means[, 2] + model_means[, 3],
         length = 0.05, angle = 90, code = 3, lwd = 1.5)
  dev.off()
}