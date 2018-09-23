desc_stats <- function(variables_dataframe, output_dir_graphs) {
  
  ## This function takes a dataframe with columns representing variables as input (first two being spatial
  ## coordinates), and calculates for each variable descriptive statistics: N, min, mean, max, sd and cv.
  ## Next, these descriptive statistics are stored as columns in a second dataframe with variable names
  ## as row names. In the end, the dataframe is written in the given output folder.
  
  ## The arguments of this function are:
  ## variables_dataframe: a dataframe with variables from which descriptive statistics should be determined.
  ## output_dir_graphs: the output directory where the descriptive statistics table is stored.
  
  ###########################################################################################################
  
  # initiate loop to store descriptive statistics of all variables in one dataframe
  for (i in 3:length(names(variables_dataframe))) {
    
    # create variable object
    var <- variables_dataframe[, i]
    
    # determine descriptive statistics
    var_n <- length(var)
    var_min <- min(var)
    var_mean <- mean(var)
    var_max <- max(var)
    var_sd <- sd(var)
    var_cv <- var_sd / var_mean
    
    # create dataframe to store descriptive statistics of all variables
    desc_stats <- round(data.frame(var_n, var_min, var_mean, var_max, var_sd, var_cv), digits = 3)
    rownames(desc_stats) <- names(variable_df)[i]
    colnames(desc_stats) <- c("N", "min", "mean", "max", "sd", "cv")
    
    # extent dataframe to store descriptive statistics of all variables
    if (i == 3) {
      desc_stats_df <- desc_stats
    } else {
      desc_stats_df <- rbind(desc_stats_df, desc_stats)
    }
  }
  
  # write descriptive statistics table in output folder
  write.table(desc_stats_df, file = paste0(output_dir_graphs, "/Descriptive statistics.txt"),
              col.names = NA, quote = FALSE, sep = ",")
}