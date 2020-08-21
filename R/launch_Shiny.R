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
  if (!exists("md") || (!exists("daf"))) {
    stop("Prior to launching the shiny application, users need to load variables as shown",
         " in the vignette Introduction_to_CytofRUV.Rmd. This error is thrown when variables",
         " such as 'md' and 'daf' have not been defined.")
  }

  if (exists("daf")) {
    if (!("cluster_codes" %in% names(metadata(daf)))) {
      stop("The app was launched without running clustering. Clustering data has to be loaded into the variable 'daf'. Please refer to the vignette",
           " Introduction_to_CytofRUV.Rmd for instructions about how to run clustering.")
    }
    else if (!("TSNE" %in% reducedDimNames(daf))) {
      stop("The app was launched without running TSNE dimension reduction. TSNE data has to be loaded into the variable 'daf'. Please refer to the vignette",
           " Introduction_to_CytofRUV.Rmd for instructions about how to run them")
    }
    else if (!("UMAP" %in% reducedDimNames(daf))) {
      stop("The app was launched without running UMAP dimension reduction. UMAP data has to be loaded into the variable 'daf'. Please refer to the vignette",
           " Introduction_to_CytofRUV.Rmd for instructions about how to run them")
    }
    else if (!exists("sub_daf")) {
      stop("Prior to launching the shiny application, users need to load variables as shown",
           " in the vignette Introduction_to_CytofRUV.Rmd. This error is thrown when 'sub_daf'",
           " has not been defined.")
    }
  }

  # Launch GUI
  shiny::runApp(
    appDir=system.file("shinyGUI", package="CytofRUV"),
    launch.browser=TRUE)
}
