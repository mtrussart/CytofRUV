shinyServer(function(input, output, session) {
  # If you want to break up the file, call source functions here, local = TRUE

  session$onSessionEnded(function() {stopApp()})

  #source("server-diagnostic_plots.R", local = TRUE)
  # =================================================================================================
  # Display Normalization Pages:
  # -------------------------------------------------------------------------------------------------
  output$tab_plots <- renderUI({tab_plots})

  # =================================================================================================
  # Define Helper Functions:
  # -------------------------------------------------------------------------------------------------
  same_elements <- function(a, b) return(identical(sort(a), sort(b)))

  # =================================================================================================
  # Save Default Plot Values:
  # -------------------------------------------------------------------------------------------------

  def <- reactiveValues(
    choiceMDS        = "condition",
    MDS_update_text = "",

    choiceExprsParam = "condition",
    choiceExprsClass = sub_daf_type,
    Exprs_update_text = "",
    Exprs_patient = levels(md$patient_id)[[1]],
    Exprs_ant = levels(panel$marker_class)[[1]],

    choice_TSNE_Colour_By1 = cluster_var,
    TSNE_update_colour_by  = cluster_var,
    TSNE_update_text = "",
    TSNE_ant = panel$antigen[[1]],

    choice_UMAP_Colour_By1 = cluster_var,
    UMAP_update_colour_by  = cluster_var,
    UMAP_update_text = "",
    UMAP_ant = panel$antigen[[1]],

    # For faceted plots: to check if sampleIDs have changed, we need to compare the contents of the input and saved lists of
    # Sample_ids
    choice_TSNE_facet_colourBy   = sampleID_sorted[1:10],
    choice_TSNE_Facet_Ant_Choice = cluster_var,
    TSNE_facet_update_colour_by = cluster_var,
    TSNE_facet_update_ant = panel$antigen[[1]],
    TSNE_facet_update_text = "",

    choice_UMAP_facet_colour_by   = sampleID_sorted[1:10],
    choice_UMAP_Facet_Ant_Choice  = cluster_var,
    UMAP_facet_update_colourby = cluster_var,
    UMAP_facet_update_ant = panel$antigen[[1]],
    UMAP_facet_update_text = ""
  )
  # =================================================================================================
  # Page 1: Median Intensities
  # -------------------------------------------------------------------------------------------------
  source("server-page1-Median Intensities.R", local=TRUE)

  # =================================================================================================
  # Page 2: Markers Distribution
  # -------------------------------------------------------------------------------------------------
  source("server-page2-Markers Distribution.R", local=TRUE)

  # =================================================================================================
  # Page 3: Clustering Results
  # -------------------------------------------------------------------------------------------------
  source("server-page3-Clustering Results.R", local=TRUE)

  # =================================================================================================
  # Page 4: Cluster Proportions
  # -------------------------------------------------------------------------------------------------
  source("server-page4-Cluster Proportions.R", local=TRUE)
})
