#' @rdname launchShiny
#' @title Launch Shiny
#'
#' is used to  launch the CytofRUV Shiny app. It examines any batch effects present when comparing CyTOF data
#' from samples replicated across batches using four different diagnostics plots:
#' Median Protein Expression, Protein Expression Distributions, Clustering Results
#' and Cluster Proportions.
#'
#' @return Opens a browser window with an interactive Shiny application.
#'
#' @export
#'

launch_Shiny<- function(){
  if (!exists("data") || (!exists("daf") || !exists("sampleID_sorted"))) {
    stop("Prior to launching the shiny application, the data needs to be loaded. Please follow the instructions in the vignette Introduction_to_CytofRUV.Rmd that explains step by step how to load the data.\n\nNote: This error assumes the output of the function load_data() was saved into a variable called 'data' which is accessible in the global environment. It also makes sure that appropriate variables like 'daf' and 'sampleID_sorted' has been defined.")
  }

  # Launch GUI
  shiny::runApp(
      appDir=system.file("shinyGUI", package="CytofRUV"),
      launch.browser=TRUE)
}
