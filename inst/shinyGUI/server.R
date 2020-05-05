shinyServer(function(input, output, session) {
  #source("server-diagnostic_plots.R", local = TRUE)
  # =================================================================================================
  # Display Normalisation Pages:
  # -------------------------------------------------------------------------------------------------
  output$tab_plots <- renderUI({tab_plots})

  # =================================================================================================
  # Save Default Plot Values:
  # -------------------------------------------------------------------------------------------------
  def <- reactiveValues(
    choiceMDS        = "condition",

    choiceExprsParam = "condition",
    choiceExprsClass = sub_daf_type,

    choice_TSNE_Colour_By1   = "meta20",

    choice_UMAP_Colour_By1   = "meta20",

    choice_TSNE_facet_colourBy  = sampleID_sorted[1:10],
    choice_TSNE_Facet_Ant_Choice  ="meta20",

    choice_UMAP_facet_colour_by   = sampleID_sorted[1:10],
    choice_UMAP_Facet_Ant_Choice   ="meta20"
  )

  # =================================================================================================
  # Page 1: Median Intensities
  # -------------------------------------------------------------------------------------------------
  observeEvent(input$mds, {
    def$choiceMDS=input$choiceMDS
  })
  # Page 1 plot 1:
  mds <- reactive({
    plotMDS(sub_daf, color_by = def$choiceMDS) +
      theme(axis.text=element_text(size=12),
            axis.title = element_text(size = 14),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12)
      ) + scale_color_manual(values = c("#0072B2","#D55E00"))
  })

  output$plotMDS <- renderPlot({
    req(mds())
    mds()
  })

  output$download_mds <- downloadHandler(
    filename = function() {
      paste("mds_plot", input$mds_tag, sep=".")
    },
    content = function(file) {
      req(mds())
      ggsave(file, plot = mds(), device = input$mds_tag)
    }
  )

  # Page 1 plot 2:
  dendogram <- reactive({
    plotExprHeatmap(sub_daf, bin_anno = FALSE, row_anno = TRUE)
  })

  output$plotDendogram <- renderPlot({
    req(dendogram())
    dendogram()
  })

  output$download_dendogram <- downloadHandler(
    filename = function() {
      paste("dendogram_plot", input$dendogram_tag, sep=".")
    },
    content = function(file) {
      req(dendogram())
      if (input$dendogram_tag == "pdf") {
        pdf(file)
      } else {
        png(file)
      }
      ComplexHeatmap::draw(dendogram())
      dev.off()
    }
  )

  # =================================================================================================
  # Page: Markers Distribution
  # -------------------------------------------------------------------------------------------------
  # First selectInput box choices: PatientIDS
  output$exprs2 <- renderUI({
    if (!input$exprs1 == "sample_id") return(NULL)
    selectInput("exprs2", "Select the patient:",
                # Unroll patient ID from a list of patient_IDs
                choices = md$patient_id)
  })

  # Second selectInput box appear when sample_id is selected in the first box, choices: antigens.
  output$exprs3 <- renderUI({
    if (!input$exprs1 == "sample_id") return(NULL)
    selectInput("exprs3", "Class of Antigen:",
                choices = levels(panel$marker_class))
  })

  # Provides appropriate data following change of parameters.
  daf_temp <- reactive({
    if (input$exprs1 == "condition" | is.null(input$exprs2) | is.null(input$exprs3)) return(sub_daf_type)
    # Marker Class: State
    if (input$exprs3 == "state") {
      patient_type = sub_daf_type[, sample_ids(sub_daf_type)%in%patient_ids(input$exprs2, sub_daf_type)]
      return(patient_type)
    }
    # Marker Class: Type
    patient_state = sub_daf_state[, sample_ids(sub_daf_state)%in%patient_ids(input$exprs2, sub_daf_state)]
    return(patient_state)
  })


  # Determines Colour-By
  observeEvent(input$exprsPlot, {
    if (input$exprs1=="condition" | is.null(input$exprs2) | is.null(input$exprs3)) {
      def$choiceExprsParam="condition"
    } else {
      def$choiceExprsParam="sample_id"
    }
  })

  # Update Plot Button
  observeEvent(input$exprsPlot, {def$choiceExprsClass=daf_temp()})

  # Plot
  exprsPlot <- reactive({
    plotExprs(def$choiceExprsClass, color_by = def$choiceExprsParam) +
      theme(axis.text=element_text(size=12),
            axis.title = element_text(size = 14),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12)
      )
  })

  output$exprsPlot  <- renderPlot({
    req(exprsPlot())
    exprsPlot()
  })

  output$download_exprsPlot <- downloadHandler(
    filename = function() {
      paste("Distribution_of_protein_expression", input$exprsPlot_tag, sep=".")
    },
    content = function(file) {
      req(exprsPlot())
      ggsave(file, plot = exprsPlot(), device = input$exprsPlot_tag)
    }
  )

  # =================================================================================================
  # Page: Clustering Results
  # -------------------------------------------------------------------------------------------------

  # Plot 1 Page 3 -
  cluster_heatmap <- reactive({
    plotClusterHeatmap(sub_daf, hm2 = NULL, k = "meta20", m = NULL, cluster_anno = TRUE, draw_freqs = TRUE)
  })

  output$cluster_heatmap <- renderPlot({
    req(cluster_heatmap())
    cluster_heatmap()
  })

  output$download_cluster_Heatmap <- downloadHandler(
    filename = function() {
      paste("ClusterHeatmap_plot", input$cluster_Heatmap_tag, sep=".")
    },
    content = function(file) {
      req(cluster_heatmap())
      if (input$cluster_Heatmap_tag == "pdf") {
        pdf(file)
      } else {
        png(file)
      }
      ComplexHeatmap::draw(cluster_heatmap()[[1]])
      dev.off()
    }
  )

  # Reactive title for plot 2 page 3 -----------------------
  TSNE_TEXT1 <- reactive({
    if (def$choice_TSNE_Colour_By1 == "meta20") {
      return("Clusters")
    }
    return(def$choice_TSNE_Colour_By1)
  })

  # Plot 2 Page 3 Antigen Selection
  output$TSNE_Ant_Choice1 <- renderUI({
    if (!input$TSNE_Colour_By1 == "Antigen") return(NULL)
    selectInput("TSNE_Ant_Choice1", "Select Antigen:", panel$antigen)
  })

  # Returns String to determine Colour-By in Plot 2 Page 3
  TSNE_grouping1 <- reactive({
    if (input$TSNE_Colour_By1 == "Antigen" & !is.null(input$TSNE_Ant_Choice1)) {
      return(input$TSNE_Ant_Choice1)
    }
    return(input$TSNE_Colour_By1)
  })

  # Update Button:
  observeEvent(input$update_TSNE1, {def$choice_TSNE_Colour_By1 = TSNE_grouping1()})

  # Define Plot 2 Page 3 -
  output$TSNE_TEXT1 <- renderText(paste0("TSNE: Coloured By ", TSNE_TEXT1()))
  plot_TSNE1 <- reactive({
    plotDR(daf, "TSNE", color_by = def$choice_TSNE_Colour_By1) +
      theme(axis.text=element_text(size=12),
            axis.title = element_text(size = 14),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12))
  })
  output$plot_TSNE1 <- renderPlot({
    req(plot_TSNE1())
    plot_TSNE1()
  })

  # Download Button
  output$download_TSNE_1 <- downloadHandler(
    filename = function() {
      paste(paste0("TSNE: Coloured By ", TSNE_TEXT1()), input$TSNE_1_tag, sep=".")
    },
    content = function(file) {
      req(plot_TSNE1())
      ggsave(file, plot = plot_TSNE1(), device = input$TSNE_1_tag)
    }
  )

  # Plot 3 page 3 -----------------------
  text_UMAP_1 <- reactive({
    if (def$choice_UMAP_Colour_By1 == "meta20") {
      return("Clusters")
    }
    return(def$choice_UMAP_Colour_By1)
  })

  # Antigen Selection
  output$UMAP_Ant_Choice1 <- renderUI({
    if (!input$Umap_Colour_By1 == "Antigen") return(NULL)
    selectInput("UMAP_Ant_Choice1", "Select Antigen:", panel$antigen)
  })

  # String to determine Colour-By
  UMAP_grouping1 <- reactive({
    if (input$Umap_Colour_By1 == "Antigen" & !is.null(input$UMAP_Ant_Choice1)) {
      return(input$UMAP_Ant_Choice1)
    }
    return(input$Umap_Colour_By1)
  })

  # Update Button:
  observeEvent(input$update_UMAP_1, {def$choice_UMAP_Colour_By1 = UMAP_grouping1()})

  # Define Plot 3 Page 3 -
  output$Umap_text_1 <- renderText(paste0("UMAP: Coloured By ", text_UMAP_1()))
  plot_UMAP1 <- reactive({
    plotDR(daf, "UMAP", color_by = def$choice_UMAP_Colour_By1) +
      theme(axis.text=element_text(size=12),
            axis.title = element_text(size = 14),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12))
  })

  output$plot_UMAP1 <- renderPlot({
    req(plot_UMAP1())
    plot_UMAP1()
  })

  # Download Button
  output$download_Umap_1 <- downloadHandler(
    filename = function() {
      paste(paste0("UMAP: Coloured By ", text_UMAP_1()), input$Umap_1_tag, sep=".")
    },
    content = function(file) {
      req(plot_UMAP1())
      ggsave(file, plot = plot_UMAP1(), device = input$Umap_1_tag)
    }
  )

  # Plot 4 Page 3 -----------------------------
  # Reactive title
  TSNE_facet_Text <- reactive({
    if (def$choice_TSNE_Facet_Ant_Choice == "meta20") {
      return("Clusters")
    }
    return(paste0("Antigen - ", def$choice_TSNE_Facet_Ant_Choice))
  })

  # Antigen Selection
  output$TSNE_Facet_Ant_Choice <- renderUI({
    if (!input$TSNE_facet_colourBy == "Antigen") return(NULL)
    selectInput("TSNE_Facet_Ant_Choice", "Select Antigen:", panel$antigen)
  })

  # String to determine Colour-by
  DR_grouping2 <- reactive({
    if (input$TSNE_facet_colourBy == "Antigen" & !is.null(input$TSNE_Facet_Ant_Choice)) {
      return(input$TSNE_Facet_Ant_Choice)
    }
    return(input$TSNE_facet_colourBy)
  })

  # Update Button
  observeEvent(input$update_TSNE_facet, {def$choice_TSNE_facet_colourBy=input$checkBox_TSNE})
  observeEvent(input$update_TSNE_facet, {def$choice_TSNE_Facet_Ant_Choice=DR_grouping2()})

  # Checbox Deselect All Button:
  observeEvent(input$deselectAll_TSNE, {
    updateCheckboxGroupInput(session,
                             "checkBox_TSNE",
                             selected = list())
  })

  # CheckBoxGroup inputs: Limit to 10 Sample_IDs to select.
  observe({
    my_min <- 0
    my_max <- 10
    if(length(input$checkBox_TSNE) > my_max){
      updateCheckboxGroupInput(session,
                               "checkBox_TSNE",
                               selected = tail(input$checkBox_TSNE,my_max))
    }
  })

  # Plot Title
  output$TSNE_facet <- renderText(paste0("TSNE: Coloured By ", UMAP_facet_Text(), ", faceted by sample_id"))

  # Define Plot
  plotTSNE_facet <- reactive({
    plotDR(daf[, sample_ids(daf)%in%def$choice_TSNE_facet_colourBy],
           "TSNE",
           color_by = def$choice_TSNE_Facet_Ant_Choice) + facet_wrap("sample_id") +
      theme(axis.text=element_text(size=12),
            axis.title = element_text(size = 14),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12))
  })

  output$plotTSNE_facet <- renderPlot({
    req(plotTSNE_facet())
    plotTSNE_facet()
  })

  # Download Button
  output$download_TSNE_facet <- downloadHandler(
    filename = function() {
      paste(paste0("TSNE: Coloured By ", textDR_2(), ", faceted by sample_id"), input$TSNE_facet_tag, sep=".")
    },
    content = function(file) {
      req(plotTSNE_facet())
      ggsave(file, plot = plotTSNE_facet(), device = input$TSNE_facet_tag)
    }
  )

  # Plot 5 page 3 -----------------------------
  UMAP_facet_Text <- reactive({
    if (def$choice_UMAP_Facet_Ant_Choice == "meta20") {
      return("Clusters")
    }
    return(paste0("Antigen - ", def$choice_UMAP_Facet_Ant_Choice))
  })

  # Antigen Selection
  output$umap_antigen_choice2 <- renderUI({
    if (!input$UMAP_facet_colour_by == "Antigen") return(NULL)
    selectInput("UMAP_Facet_Ant_Choice", "Select Antigen:", panel$antigen)
  })

  # String to determine Colour-by
  UMAP_grouping2 <- reactive({
    if (input$UMAP_facet_colour_by == "Antigen" & !is.null(input$UMAP_Facet_Ant_Choice)) {
      return(input$UMAP_Facet_Ant_Choice)
    }
    return(input$UMAP_facet_colour_by)
  })

  # Update Button
  observeEvent(input$update_UMAP_facet, {def$choice_UMAP_facet_colour_by=input$checkBox_UMAP})
  observeEvent(input$update_UMAP_facet, {def$choice_UMAP_Facet_Ant_Choice=UMAP_grouping2()})

  # Checbox Deselect All Button:
  observeEvent(input$deselectAll_UMAP, {
    updateCheckboxGroupInput(session,
                             "checkBox_UMAP",
                             selected = list())
  })

  # CheckBoxGroup inputs: Limit to 10 Sample_IDs to select.
  observe({
    my_min <- 0
    my_max <- 10
    if(length(input$checkBox_UMAP) > my_max){
      updateCheckboxGroupInput(session,
                               "checkBox_UMAP",
                               selected = tail(input$checkBox_UMAP,my_max))
    }
  })
  # Plot Title
  output$UMAP_facet_Text <- renderText(paste0("UMAP: Coloured By ", UMAP_facet_Text(), ", faceted by sample_id"))

  # Define Plot
  plot_UMAP_facet <- reactive({
    plotDR(daf[, sample_ids(daf)%in%def$choice_UMAP_facet_colour_by], "UMAP", color_by = def$choice_UMAP_Facet_Ant_Choice) +
      facet_wrap("sample_id") +
      theme(axis.text=element_text(size=12),
            axis.title = element_text(size = 14),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12))
  })

  output$plot_UMAP_facet <- renderPlot({
    req(plot_UMAP_facet())
    plot_UMAP_facet()
  })

  # Download Button
  output$download_UMAP_facet <- downloadHandler(
    filename = function() {
      paste(paste0("UMAP: Coloured By ", UMAP_facet_Text(), ", faceted by sample_id"), input$UMAP_2_tag, sep=".")
    },
    content = function(file) {
      req(plot_UMAP_facet())
      ggsave(file, plot = plot_UMAP_facet(), device = input$UMAP_2_tag)
    }
  )
  # =================================================================================================
  # Page: Cluster Proportions
  # -------------------------------------------------------------------------------------------------
  Abundance_cluster <- reactive({
    daf$sample_id<-factor(daf$sample_id,levels = sampleID_sorted)
    plotAbundances(daf, k = "meta20", by = "sample_id") +
      theme(axis.text=element_text(size=12),
            axis.title = element_text(size = 14),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            strip.text = element_blank()) +
      facet_wrap(facets = NULL, scales="fixed")
  })

  output$Abundance_cluster <- renderPlot({
    req(Abundance_cluster())
    Abundance_cluster()
  })

  output$download_Abundance_cluster <- downloadHandler(
    filename = function() {
      paste("Cluster proportions across samples", input$Abundance_cluster_tag, sep=".")
    },
    content = function(file) {
      req(Abundance_cluster())
      ggsave(file, plot = Abundance_cluster(), device = input$Abundance_cluster_tag)
    }
  )
})
