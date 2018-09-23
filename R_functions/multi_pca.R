multi_pca <- function(pca_variables, output_dir_graphs, num_PC = 0, eucl_dist = 0, screeplot = FALSE) {
  
  ## This function consists of two parts. The first part is intended to perform principal component analysis (PCA) and create a screeplot
  ## to determine the number of principal components to retain. The second part of the function is used to perform the actual PCA and
  ## Multispati PCA based on the number of components to retain, and to store the actual principal components, loadings, and loading plots in
  ## data frames and in the given output folder.
  
  ## The arguments of this function are:
  ## pca_variables: the variables (including coordinates) for which PCA and MPCA should be performed.
  ## output_dir_graphs: the output directory where the PCA outputs are stored.
  ## num_PC (only valid for MPCA): the number of principal components to retain.
  ## eucl_dist (only valid for MPCA): the Euclidean distance [m] used to identify the number of neighbouring points.
  ## screeplot (logical): if TRUE, PCA will only be performed to create a screeplot, else the complete PCA and MPCA analyses will be performed.
  
  #############################################################################################################################################
  
  # select variables excluding the coordinates
  pca_var <- pca_variables[, 3:ncol(pca_variables)]
  
  if (screeplot == TRUE) {
    # perform PCA on complete PCA variables dataframe for plotting a screeplot
    dudi_pca <- dudi.pca(pca_var, center = TRUE, scale = TRUE, scannf = FALSE,  nf = length(pca_var))
    
    # calculate Variance Accounted For (VAF) for each principal component
    VAF <- vector()
    for (i in 1:length(dudi_pca$eig)) {
      VAF[i] <- round(dudi_pca$eig[i] * 100 / length(pca_var), digits = 2)
    }
    
    # plotting PCA screeplot
    jpeg(filename = paste0(output_dir_graphs, "/PCA scree plot.jpg"), width = 320, height = 320)
    barplot(VAF, names.arg = 1:length(dudi_pca$eig), main = "PCA scree plot", xlab = "principal components",
            ylab = "variance accounted for (%)", col ="gray", ylim = c(0, 70))
    lines(x = 1:length(VAF) * 1.2 - 0.5, y = VAF, type = "b", pch = 16, col = "blue")
    dev.off()
      
  } else {
    # perform PCA and MPCA on complete PCA variables dataframe, and create output
    dudi_pca <- dudi.pca(pca_var, center = TRUE, scale = TRUE, scannf = FALSE,  nf = num_PC)
    
    # identify neighbouring points by Eucledian distance, calculate spatial weights and store them in a list
    coords <- coordinates(pca_variables[, 1:2])
    neighbours <- dnearneigh(coords, 0, eucl_dist)
    neighbours_list <- nb2listw(neighbours, style = "W", zero.policy = TRUE)
    
    # perform multivariate spatial analysis based on dudi_pca object and neighbours list
    ms_pca <- multispati(dudi_pca, neighbours_list, scannf = FALSE, nfposi = num_PC)
    
    # extract and replace eigenvalues (eig), loadings (c1) and principal components (li)
    dudi_pca$eig <- ms_pca$eig
    dudi_pca$c1 <- ms_pca$c1
    dudi_pca$li <- ms_pca$li
    
    # calculate Variance Accounted For (VAF) for each principal component
    VAF <- vector()
    for (i in 1:length(dudi_pca$eig)) {
      VAF[i] <- round(dudi_pca$eig[i] * 100 / length(pca_var), digits = 2)
    }
    
    # initiate nested loop to plot PCA biplots for each pair of PCs
    for (i in 1:length(dudi_pca$c1)) {
      for (j in 1:length(dudi_pca$c1)) {
        if (j > i) {
          
          # create pca biplot
          pca_biplot <- fviz_pca_var(X = dudi_pca, axes = c(i, j), title = paste0("PCA biplot sPC", i, " vs sPC", j),
                                     xlab = paste0("sPC", i, " (", VAF[i], "%)"), ylab = paste0("sPC", j, " (", VAF[j], "%)"),
                                     repel = TRUE, labelsize = 4)
          
          # store biplots in output folder
          jpeg(filename = paste0(output_dir_graphs, "/PCA biplot sPC", i, " vs sPC", j, ".jpg"), width = 320, height = 320)
          plot(pca_biplot)
          dev.off()
        }
      }
    }
    
    # calculate variable communalities and store them together with PCs into one data frame
    communalities <- round(apply(dudi_pca$c1^2, MARGIN = 1, FUN = sum), digits = 3)
    loadings_table <- cbind(communalities, round(dudi_pca$c1, digits = 3))
    
    # name columns of PCA loadings data frame
    for (i in 1:ncol(loadings_table)) {
      if (i == 1) {
        names(loadings_table)[i] <- (paste("communalities of first", num_PC, "sPCs"))
      } else {
        names(loadings_table)[i] <- paste0("sPC", i-1, " (", VAF[i-1], "%)")
      }
    }
    
    # write PCA loadings table in output folder
    write.table(loadings_table, file = paste0(output_dir_graphs, "/PCA loadings.txt"),
                col.names = NA, quote = FALSE, sep = ",")
    
    # add spatial principal components to PCA variables dataframe and rename columns
    pca_variables_with_pcs <- cbind(pca_variables, dudi_pca$li[, 1:num_PC])
    colnames(pca_variables_with_pcs)[(ncol(pca_variables) + 1):ncol(pca_variables_with_pcs)] <- paste0("sPC", 1:num_PC)
    return(pca_variables_with_pcs)
  }
}