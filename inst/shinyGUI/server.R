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
    choiceMDS="condition",
    choiceExprsParam="condition",
    choiceExprsClass=sub_daf_type,
    choice_plotDR1="meta20",
    choice_plotDR2a=md$sample_id[1:10],
    choice_plotDR2b="meta20"
  )

  # =================================================================================================
  # Page 1: Median Intensities
  # -------------------------------------------------------------------------------------------------
  observeEvent(input$mds, {def$choiceMDS=input$choiceMDS})
  # Page 1 plot 1:
  mds <- reactive({
    plotMDS(sub_daf, color_by = def$choiceMDS) +
      theme(axis.text=element_text(size=12),
            axis.title = element_text(size = 14),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12)
      )
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
      paste("exprs_plot", input$exprsPlot_tag, sep=".")
    },
    content = function(file) {
      req(exprsPlot())
      ggsave(file, plot = exprsPlot(), device = input$exprsPlot_tag)
    }
  )

  # =================================================================================================
  # Page: Clustering Results
  # -------------------------------------------------------------------------------------------------
  # Plot 1 Page 3 Antigen Selection
  output$antigen_choice1 <- renderUI({
    if (!input$DR1a == "Antigen") return(NULL)
    selectInput("DR1b", "Select Antigen:", panel$antigen)
  })

  # Returns String to determine Colour-By in Plot 1 Page 3
  DR_grouping1 <- reactive({
    if (input$DR1a == "Antigen" & !is.null(input$DR1b)) {
      return(input$DR1b)
    }
    return(input$DR1a)
  })

  # Action Button Control:
  observeEvent(input$plotDR_1, {def$choice_plotDR1 = DR_grouping1()})

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

  # Reactive title for plot 2 page 3
  textDR_1 <- reactive({
    if (def$choice_plotDR1 == "meta20") {
      return("Clusters")
    }
    return(def$choice_plotDR1)
  })

  # Plot 2 Page 3 -
  output$textDR_1 <- renderText(paste0("TSNE: Coloured By ", textDR_1()))
  plotDR_1 <- reactive({
    plotDR(daf, "TSNE", color_by = def$choice_plotDR1) +
      theme(axis.text=element_text(size=12),
            axis.title = element_text(size = 14),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12))
  })
  output$plotDR_1 <- renderPlot({
    req(plotDR_1())
    plotDR_1()
  })

  output$download_plotDR_1 <- downloadHandler(
    filename = function() {
      paste(paste0("TSNE: Coloured By ", textDR_1()), input$plotDR_1_tag, sep=".")
    },
    content = function(file) {
      req(plotDR_1())
      ggsave(file, plot = plotDR_1(), device = input$plotDR_1_tag)
    }
  )

  # Plot 3 Page 3 Antigen Selection
  output$antigen_choice2 <- renderUI({
    if (!input$DR2a == "Antigen") return(NULL)
    selectInput("DR2b", "Select Antigen:", panel$antigen)
  })

  # Returns String to determine Colour-by in Plot 3 Page 3
  DR_grouping2 <- reactive({
    if (input$DR2a == "Antigen" & !is.null(input$DR2b)) {
      return(input$DR2b)
    }
    return(input$DR2a)
  })

  # Action Button Control:
  observeEvent(input$plotDR_2, {def$choice_plotDR2a=input$checkBox_sampleIDs})
  observeEvent(input$plotDR_2, {def$choice_plotDR2b=DR_grouping2()})

  # Reactive title for plot 3 page 3
  textDR_2 <- reactive({
    if (def$choice_plotDR2b == "meta20") {
      return("Clusters")
    }
    return(paste0("Antigen - ", def$choice_plotDR2b))
  })

  # Action Button Control:
  observeEvent(input$deselectAll, {
    updateCheckboxGroupInput(session,
                             "checkBox_sampleIDs",
                             selected = list())
  })

  # CheckBoxGroup inputs: Limit to 10 Sample_IDs to select.
  observe({
    my_min <- 0
    my_max <- 10
    if(length(input$checkBox_sampleIDs) > my_max){
      updateCheckboxGroupInput(session,
                               "checkBox_sampleIDs",
                               selected = tail(input$checkBox_sampleIDs,my_max))
    }
  })

  # Plot 3 Page 3 -
  output$textDR_facet <- renderText(paste0("TSNE: Coloured By ", textDR_2(), ", faceted by sample_id"))

  plotDR_facet <- reactive({
     plotDR(daf[, sample_ids(daf)%in%def$choice_plotDR2a],
           "TSNE",
           color_by = def$choice_plotDR2b) + facet_wrap("sample_id") +
      theme(axis.text=element_text(size=12),
            axis.title = element_text(size = 14),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12))
  })

  output$plotDR_facet <- renderPlot({
    req(plotDR_facet())
    plotDR_facet()
  })

  output$download_plotDR_facet <- downloadHandler(
    filename = function() {
      paste(paste0("TSNE: Coloured By ", textDR_2(), ", faceted by sample_id"), input$plotDR_facet_tag, sep=".")
    },
    content = function(file) {
      req(plotDR_facet())
      ggsave(file, plot = plotDR_facet(), device = input$plotDR_facet_tag)
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
      paste("Abundance_cluster", input$Abundance_cluster_tag, sep=".")
    },
    content = function(file) {
      req(Abundance_cluster())
      ggsave(file, plot = Abundance_cluster(), device = input$Abundance_cluster_tag)
    }
  )
})
