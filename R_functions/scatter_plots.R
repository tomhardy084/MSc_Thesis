scatter_plots <- function(variables_dataframe, output_dir_scatter) {
  
  ## This function takes a dataframe with columns representing variables as input (first two being spatial
  ## coordinates) and generates density scatterplots for each unique pair of variables, including a regression
  ## line and correlation coefficient. In the end, the scatterplots are stored in the given output folder.
  
  ## The arguments of this function are:
  ## variables_dataframe: a dataframe with variables from which scatterplots should be generated.
  ## output_dir_scatter: the output directory where the scatterplots are stored.
  
  ###########################################################################################################
  
  # initiate nested loop to generate and store density scatterplots for each pair of variables
  for (i in 3:length(names(variables_dataframe))) {
    for (j in 3:length(names(variables_dataframe))) {
      if (j > i) {
        # create variable objects and variable names
        var_1 <- variables_dataframe[, i]
        var_2 <- variables_dataframe[, j]
        varname_1 <- gsub("_", " ", names(variables_dataframe)[i])
        varname_2 <- gsub("_", " ", names(variables_dataframe)[j])
        
        # create linear model and determine correlation coefficient
        lm_fit <- lm(var_2 ~ var_1)
        cor_coef <- round(sqrt(summary(lm_fit)$r.squared), digits = 3)
        if (summary(lm_fit)$coefficients[2, 1] < 0) {
          cor_coef <- (-1) * cor_coef
        }
        
        # plot and store scatterplots in output folder
        col_ramp <- colorRampPalette(c("blue", "cyan", "green", "yellow", "orange", "red"))
        dens_cols <- densCols(var_1, var_2, colramp = col_ramp)
        jpeg(filename = paste0(output_dir_scatter, "/Scatterplot ", varname_1, " vs ", varname_2, ".jpg"),
             width = 480, height = 480)
        plot(var_1, var_2, col = "white", main = paste("Density scatterplot, r =", cor_coef),
             xlab = varname_1, ylab = varname_2, cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
        grid()
        points(var_1, var_2, col = dens_cols, pch = 1)
        abline(reg = lm_fit, col = "darkgray", lwd = 2, lty = "solid")
        dev.off()
      }
    }
  }
}