#' @rdname launchShiny
#' @title Launch Shiny
#'
#' launch_Shiny is used to  launch the CytofRUV Shiny app. It examines any batch effects present when comparing CyTOF data
#' from samples replicated across batches using four different diagnostics plots:
#' Median Protein Expression, Protein Expression Distributions, Clustering Results
#' and Cluster Proportions.
#'
#' @return Opens a browser window with an interactive Shiny application.
#'
#' @export
#'

launch_Shiny<- function(){

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
