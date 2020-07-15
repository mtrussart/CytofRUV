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
  if (!exists("data") || (length(data) != 6)) {
    stop("Prior to launching the shiny application, the data needs to be loaded. Please follow the instructions in the vignette Introduction_to_CytofRUV.Rmd that explains step by step how to load the data.\n\nNote: This error assumes the output of the function load_data() was saved into a variable called 'data' which is accessible in the global environment.")  
  }

  # Launch GUI
  shiny::runApp(
      appDir=system.file("shinyGUI", package="CytofRUV"),
      launch.browser=TRUE)

  # source("/inst/ShinyGUI/server.R")
  # source("/inst/ShinyGUI/ui.R")
  #
  # app <- shinyApp(ui = ui, server = server)
  # runApp(app,
  #        launch.browser = TRUE)
}
