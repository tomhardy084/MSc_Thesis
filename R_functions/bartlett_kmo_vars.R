bartlett_kmo_vars <- function(cor_mat, variable_df, output_dir_graphs) {
  
  ## This function performs Bartlett's test of Sphericity and calculates KMO indices
  ## for Measure of Sampling Adequacy based on an input correlation matrix. Moreover,
  ## given KMO values larger than 0.5, it selects dataframe variables that will be
  ## used for PCA and cluster analysis.
  
  ## The arguments of this function are:
  ## cor_mat: an input correlation matrix as basis for Bartlett's and KMO tests.
  ## output_crs: the project string representing the output coordinate reference system.
  ## output_dir_graphs: the output directory where the bartlett_kmo table is stored.
  
  ##############################################################################################################
  
  # determine bartlett and kmo values
  bartlett <- cortest.bartlett(R = cor_mat, n = nrow(variable_df))
  kmo_full <- KMO(cor_mat)
  
  # select variables based on KMO values larger than 0.5 and determine new KMO indices
  var_select <- which(kmo_full$MSAi >= 0.5)
  pca_variables <- variable_df[, c(1:2, (var_select + 2))]
  kmo_reduced <- KMO(cor_mat[var_select, var_select])
  
  # store Bartlett and KMO values in dataframe
  bartlett_kmo_df <- as.data.frame(matrix(round(c(bartlett$chisq, bartlett$p.value, kmo_full$MSA,
                                                  kmo_full$MSAi, kmo_reduced$MSA, kmo_reduced$MSAi),
                                                digits = 3), ncol = 1))
  
  # assign row names to dataframe
  rownames(bartlett_kmo_df) <- c("bartlett_chi2", "bartlett_pval", "kmo_full_overall",
                                 paste0("kmo_full_", names(kmo_full$MSAi)), "kmo_reduced_overall",
                                 paste0("kmo_reduced_", names(kmo_reduced$MSAi)))
  
  # write dataframe as table in output folder
  write.table(bartlett_kmo_df, file = paste0(output_dir_graphs, "/Bartlett and KMO.txt"),
              col.names = NA, quote = FALSE, sep = ",")
  
  return(pca_variables)
}