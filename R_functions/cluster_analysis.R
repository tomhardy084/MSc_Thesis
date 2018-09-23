cluster_analysis <- function(cluster_variables, num_clus = 0, output_dir_graphs, screeplot = FALSE) {
  
  ## This function consists of two parts. The first part is intended to perform cluster analysis and create a screeplot to determine the number
  ## of clusters to create. The second part of the function is used to perform the actual cluster analysis based on the predefined number of
  ## clusters, and to store the actual clusters in a data frame and in the given output folder.
  
  ## The arguments of this function are:
  ## cluster_variables: the variables (including coordinates) for which cluster analysis should be performed.
  ## num_clus (only valid for actual cluster analysis): the number of output clusters.
  ## output_dir_graphs: the output directory where the sreeplot and data frame with clusters are stored.
  ## screeplot (logical): if TRUE, cluster analysis will only be performed to create a screeplot, else the complete cluster analysis will be performed.
  
  #####################################################################################################################################################
  
  # select variables excluding the coordinates
  cluster_var <- cluster_variables[, 11:14]
  
  if (screeplot == TRUE) {
    # perform cluster analysis to create a cluster screeplot
    set.seed(5)
    centers <- 10
    wss <- vector()
    for (i in 1:centers) {
      wss[i] <- sum(kmeans(x = cluster_var, centers = i, iter.max = 25, nstart = i)$withinss)
    }
    
    # plot and store cluster screeplot in output folder
    jpeg(filename = paste0(output_dir_graphs,"/Cluster scree plot.jpg"), width = 320, height = 320)
    barplot(wss, names.arg = 1:centers, main = "Clustering scree plot", xlab = "number of clusters",
            ylab = "within group sums of squares", col ="gray")
    lines(x = 1:centers * 1.2 - 0.5, y = wss, type = "b", pch = 16, col = "blue")
    dev.off()
    
  } else {
    # initiate loop to perform cluster analysis based on selected number of clusters
    for (i in num_clus) {
      set.seed(5)
      clust_kmeans <- kmeans(cluster_var, centers = i, iter.max = 25, nstart = i)
      
      # add clusters to PCA variables dataframe
      if (!(exists(x = "pca_variables_with_clusters"))) {
        pca_variables_with_clusters <- cbind(cluster_variables, clust_kmeans$cluster)
      } else {
        pca_variables_with_clusters <- cbind(pca_variables_with_clusters, clust_kmeans$cluster)
      }
      
      # rename the column with clusters
      colnames(pca_variables_with_clusters)[ncol(pca_variables_with_clusters)] <- paste0("kmeans_", i, "_clus")
    }
    
    # write data frame as csv file in output folder (to save memory in R)
    write.csv(x = pca_variables_with_clusters, file = paste0(output_dir_graphs, "/PCA variables with clusters.csv"),
              quote = FALSE, sep = ",")
  }
}