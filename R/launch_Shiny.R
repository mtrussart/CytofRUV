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
  if (!exists("md") || (!exists("daf") || !exists("sampleID_sorted"))) {
    stop("Prior to launching the shiny application, users need to load variables as shown in the vignette Introduction_to_CytofRUV.Rmd.
         This error is thrown when variables such as 'md', 'daf' and 'sampleID_sorted' have not been defined.")
  }

  # Launch GUI
  shiny::runApp(
    appDir=system.file("shinyGUI", package="CytofRUV"),
    launch.browser=TRUE)
}
