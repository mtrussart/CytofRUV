# Define type of markers
daf_type  <- daf[SingleCellExperiment::rowData(daf)$marker_class=="type", ]
daf_state <- daf[SingleCellExperiment::rowData(daf)$marker_class=="state", ]
sub_daf_state <- daf_state[, sample(ncol(daf_state), n_subset_marker_specific)]
sub_daf_type  <- daf_type[, sample(ncol(daf_type), n_subset_marker_specific)]
# Define batch
batch_ids <- is.factor(rep(md$batch, nrow(daf)))
sampleID_sorted <- md$sample_id[order(md$patient_id)]

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
  no_sampleIds = length(sampleID_sorted)
  nb_facets = 10
  initial_sampleIDs = if(no_sampleIds < nb_facets) as.character(sampleID_sorted[1:no_sampleIds]) else as.character(sampleID_sorted[1:nb_facets])

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
    choice_TSNE_facet_colourBy   = initial_sampleIDs,
    choice_TSNE_Facet_Ant_Choice = cluster_var,
    TSNE_facet_update_colour_by = cluster_var,
    TSNE_facet_update_ant = panel$antigen[[1]],
    TSNE_facet_update_text = "",

    choice_UMAP_facet_colour_by   = initial_sampleIDs,
    choice_UMAP_Facet_Ant_Choice  = cluster_var,
    UMAP_facet_update_colourby = cluster_var,
    UMAP_facet_update_ant = panel$antigen[[1]],
    UMAP_facet_update_text = ""
  )

  # =================================================================================================
  # Page 1: Median Intensities
  # -------------------------------------------------------------------------------------------------


  ### PlotMDS function

  plot_MDS <- function (x, color_by="condition") {

    # compute medians across samples
    cs_by_s <- split(seq_len(ncol(x)), x$sample_id)
    es <- as.matrix(assay(x, "exprs"))
    ms <- vapply(cs_by_s, function(cs)
      rowMedians(es[, cs, drop = FALSE]),
      numeric(nrow(x)))
    rownames(ms) <- rownames(x)

    # get MDS coordinates
    mds <- limma::plotMDS(ms, plot=FALSE)
    df <- data.frame(MDS1=mds$x, MDS2=mds$y)

    # add metadata information
    md <- metadata(x)$experiment_info
    m <- match(rownames(df), md$sample_id)
    df <- data.frame(df, md[m, ])

    ggplot(df, aes_string(x="MDS1", y="MDS2", col=color_by)) +
      geom_label_repel(aes_string(label="sample_id"),
                       show.legend=FALSE) + geom_point(alpha=.8, size=1.2) +
      guides(col=guide_legend(overide.aes=list(alpha=1, size=3))) +
      theme_void() + theme(aspect.ratio=1,
                           panel.grid.minor=element_blank(),
                           panel.grid.major=element_line(color='lightgrey', size=.25),
                           axis.title=element_text(face='bold'),
                           axis.title.y=element_text(angle=90),
                           axis.text=element_text())
  }

  # Update Button Logic:
  observeEvent(input$mds, {
    def$choiceMDS = input$choiceMDS
    def$MDS_update_text = ""
  })

  # Logic: Reminder Text to Press Update Button.
  observeEvent(input$choiceMDS, {
    if (input$choiceMDS != def$choiceMDS){
      def$MDS_update_text <- "Press the Update Button."
    } else {
      def$MDS_update_text <- ""
    }
  })

  # Renders Reminder Text
  output$MDS_updateReminder <- renderText(def$MDS_update_text)

  # Page 1 plot 1:
  mds <- reactive({
    plot_MDS(sub_daf, color_by = def$choiceMDS) +
      theme(axis.text=element_text(size=12),
            axis.title = element_text(size = 14),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12))
    #+  (if (length(levels(sub_daf[[def$choiceMDS]])) <= 8) scale_color_manual(values = brewer.pal(n = 8, name = "Dark2")))
    # If more than 8 options, default to default ggplot colour scheme. If < 8, colour blind friendly palette is used.
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
        pdf(file, width = 8)
      } else {
        png(file, width=720, units = "px")
      }
      ComplexHeatmap::draw(dendogram())
      dev.off()
    }
  )

  # =================================================================================================
  # Page 2: Markers Distribution
  # -------------------------------------------------------------------------------------------------
  # Returns sample_IDs related to a given patient, found through patient_ID.
  patient_ids <- function(patient_id, dframe){
    return(levels(factor(sample_ids(dframe)[grepl(patient_id,sample_ids(dframe))])))
  }

  # First selectInput box choices: PatientIDS
  output$exprs2 <- renderUI({
    if (!(input$exprs1 == "sample_id")) return(NULL)
    selectInput("exprs2", "Select the patient:",
                choices = levels(md$patient_id),
                selected = levels(md$patient_id)[[1]])
  })

  # Second selectInput box appear when sample_id is selected in the first box, choices: antigens.
  output$exprs3 <- renderUI({
    if (!(input$exprs1 == "sample_id")) return(NULL)
    selectInput("exprs3", "Select the class of markers:",
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

  # Upon Update Button Press:
  observeEvent(input$exprsPlot, {
    # Determines Colour-By parameter
    if (input$exprs1=="condition" | is.null(input$exprs2) | is.null(input$exprs3)) {
      def$choiceExprsParam="condition"
    } else {
      def$choiceExprsParam="sample_id"
    }
    def$choiceExprsClass = daf_temp()
    def$Exprs_update_text = ''
    def$Exprs_patient = input$exprs2
    def$Exprs_ant = input$exprs3
  })

  # Logic for Update Reminder Text:
  observeEvent({
    input$exprs1
    input$exprs2
    input$exprs3
  },
  {
    if (input$exprs1 != def$choiceExprsParam) {
      def$Exprs_update_text <- "Press the update button."
    }
    else if (!is.null(input$exprs2) && (input$exprs2 != def$Exprs_patient)) {
      def$Exprs_update_text <- "Press the update button."
    }
    else if (!is.null(input$exprs3) && (input$exprs3 != def$Exprs_ant)) {
      def$Exprs_update_text <- "Press the update button."
    }
    else {
      def$Exprs_update_text <- ""
    }
  })

  output$Exprs_update_text <- renderText(def$Exprs_update_text)

  # Define the Plot
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
      ggsave(file, plot = exprsPlot(), device = input$exprsPlot_tag, width = 34, height = 18, units = "cm")
    }
  )

  # =================================================================================================
  # Page 3: Clustering Results
  # -------------------------------------------------------------------------------------------------
  # =================================================================================================
  #   Plot 1: Cluster Heatmap
  # -------------------------------------------------------------------------------------------------
  cluster_heatmap <- reactive({
    plotClusterHeatmap(sub_daf, hm2 = NULL, k = cluster_var, m = NULL, cluster_anno = TRUE, draw_freqs = TRUE)
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
        pdf(file, width = 10)
      } else {
        png(file, width=720, units = "px")
      }
      ComplexHeatmap::draw(cluster_heatmap()[[1]])
      dev.off()
    }
  )

  # =================================================================================================
  #   Plot 2: TSNE Non Facetted
  # -------------------------------------------------------------------------------------------------
  TSNE_TEXT1 <- reactive({
    if (def$choice_TSNE_Colour_By1 == cluster_var) {
      return("Clusters")
    } else if (def$choice_TSNE_Colour_By1 == "batch") {
      return("batch")
    }
    return(paste0("Antigen - ", def$choice_TSNE_Colour_By1))
  })

  # Antigen Selection
  output$TSNE_Ant_Choice1 <- renderUI({
    if (!(input$TSNE_Colour_By1 == "Antigen")) return(NULL)
    selectInput("TSNE_Ant_Choice1", "Select Antigen:", panel$antigen)
  })

  # String to determine Colour-By
  TSNE_grouping1 <- reactive({
    if (input$TSNE_Colour_By1 == "Antigen" & !is.null(input$TSNE_Ant_Choice1)) {
      return(input$TSNE_Ant_Choice1)
    }
    return(input$TSNE_Colour_By1)
  })

  # Update Button:
  observeEvent(input$update_TSNE1, {
    def$choice_TSNE_Colour_By1 = TSNE_grouping1()
    def$TSNE_update_colour_by = input$TSNE_Colour_By1
    def$TSNE_update_text = ""
    def$TSNE_ant = input$TSNE_Ant_Choice1
  })

  # Logic for Update Reminder Text:
  observeEvent({
    input$TSNE_Ant_Choice1
    input$TSNE_Colour_By1
  },
  {
    if (input$TSNE_Colour_By1 != def$TSNE_update_colour_by) {
      def$TSNE_update_text <- "Press the update button."
    }
    else if (!is.null(input$TSNE_Ant_Choice1) && (input$TSNE_Ant_Choice1 != def$TSNE_ant)) {
      def$TSNE_update_text <- "Press the update button."
    }
    else {
      def$TSNE_update_text <- ""
    }
  })

  # Renders Reminder Text
  output$TSNE_update_text <- renderText({ def$TSNE_update_text })

  output$TSNE_TEXT1 <- renderText(paste0("TSNE: Coloured By ", TSNE_TEXT1()))

  # Define Plot 2 Page 3 -
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

  # =================================================================================================
  #   Plot 3: UMAP Non-Facetted
  # -------------------------------------------------------------------------------------------------
  text_UMAP_1 <- reactive({
    if (def$choice_UMAP_Colour_By1 == cluster_var) {
      return("Clusters")
    } else if (def$choice_UMAP_Colour_By1 == "batch") {
      return("Batch")
    }
    return(paste0("Antigen - ", def$choice_UMAP_Colour_By1))
  })

  # Antigen Selection
  output$UMAP_Ant_Choice1 <- renderUI({
    if (!(input$Umap_Colour_By1 == "Antigen")) return(NULL)
    selectInput("UMAP_Ant_Choice1", "Select Antigen:", panel$antigen)
  })

  # String to determine Colour-By
  UMAP_grouping1 <- reactive({
    if (input$Umap_Colour_By1 == "Antigen" && (!is.null(input$UMAP_Ant_Choice1))) {
      return(input$UMAP_Ant_Choice1)
    }
    return(input$Umap_Colour_By1)
  })

  # Update Button:
  observeEvent(input$update_UMAP_1, {
    def$choice_UMAP_Colour_By1 = UMAP_grouping1()
    def$UMAP_update_colour_by = input$Umap_Colour_By1
    def$UMAP_update_text = ""
    def$UMAP_ant = input$UMAP_Ant_Choice1
  })

  # Logic for Update Reminder Text:
  observeEvent({
    input$UMAP_Ant_Choice1
    input$Umap_Colour_By1
  },
  {
    if (input$Umap_Colour_By1 != def$UMAP_update_colour_by) {
      def$UMAP_update_text <- "Press the update button."
    }
    else if (!is.null(input$UMAP_Ant_Choice1) && (input$UMAP_Ant_Choice1 != def$UMAP_ant)) {
      def$UMAP_update_text <- "Press the update button."
    }
    else {
      def$UMAP_update_text <- ""
    }
  })

  # Renders Reminder Text
  output$UMAP_update_text <- renderText({ def$UMAP_update_text })

  # Reactive Title:
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

  # =================================================================================================
  #   Plot 4: TSNE Facetted
  # -------------------------------------------------------------------------------------------------
  # Reactive title
  TSNE_facet_Text <- reactive({
    if (def$choice_TSNE_Facet_Ant_Choice == cluster_var) {
      return("Clusters")
    } else if (def$choice_TSNE_Facet_Ant_Choice == "batch") {
      return("Batch")
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
    if (input$TSNE_facet_colourBy == "Antigen" && (!is.null(input$TSNE_Facet_Ant_Choice))) {
      return(input$TSNE_Facet_Ant_Choice)
    }
    return(input$TSNE_facet_colourBy)
  })

  # Update Button
  observeEvent(input$update_TSNE_facet, {
    def$choice_TSNE_facet_colourBy = input$checkBox_TSNE
    def$choice_TSNE_Facet_Ant_Choice = DR_grouping2()

    def$TSNE_facet_update_text = ""
    def$TSNE_facet_update_ant = input$TSNE_Facet_Ant_Choice
    def$TSNE_facet_update_colour_by = input$TSNE_facet_colourBy
  })

  # Logic for Update Reminder Text:
  observeEvent({
    input$TSNE_facet_colourBy
    input$TSNE_Facet_Ant_Choice
    input$checkBox_TSNE
  },
  {
    if (input$TSNE_facet_colourBy != def$TSNE_facet_update_colour_by) {
      def$TSNE_facet_update_text <- "Press the update button."
    }
    else if ((!is.null(input$TSNE_Facet_Ant_Choice)) && (input$TSNE_Facet_Ant_Choice != def$TSNE_facet_update_ant)) {
      def$TSNE_facet_update_text <- "Press the update button."
    }
    else if (!same_elements(input$checkBox_TSNE, def$choice_TSNE_facet_colourBy)) {
      def$TSNE_facet_update_text <- "Press the update button."
    }
    else {
      def$TSNE_facet_update_text <- ""
    }
  })

  # Renders Reminder Text
  output$TSNE_facet_update_text <- renderText({ def$TSNE_facet_update_text })

  # Checbox Deselect All Button:
  observeEvent(input$deselectAll_TSNE, {
    updateCheckboxGroupInput(session, "checkBox_TSNE", selected = list())
    if (length(def$choice_TSNE_facet_colourBy) > 0) {def$TSNE_facet_update_text <- "Select at least one sample id."}
  })

  # CheckBoxGroup inputs: Limit to 10 Sample_IDs to select.
  observe({
    my_min <- 0
    my_max <- 10
    if(length(input$checkBox_TSNE) > my_max){
      updateCheckboxGroupInput(session, "checkBox_TSNE", selected = tail(input$checkBox_TSNE,my_max))
    }
  })

  # Plot Title
  output$TSNE_facet_Text <- renderText(paste0("TSNE: Coloured By ", TSNE_facet_Text(), ", Separated by Sample_id"))


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
      paste(paste0("TSNE: Coloured By ", TSNE_facet_Text(), ", Separated by Sample Ids"), input$TSNE_facet_tag, sep=".")
    },
    content = function(file) {
      req(plotTSNE_facet())
      ggsave(file, plot = plotTSNE_facet(), device = input$TSNE_facet_tag)
    }
  )

  # =================================================================================================
  #   Plot 5: UMAP Facetted
  # -------------------------------------------------------------------------------------------------
  UMAP_facet_Text <- reactive({
    if (def$choice_UMAP_Facet_Ant_Choice == cluster_var) {
      return("Clusters")
    } else if (def$choice_UMAP_Facet_Ant_Choice == "batch") {
      return("Batch")
    }
    return(paste0("Antigen - ", def$choice_UMAP_Facet_Ant_Choice))
  })

  # Antigen Selection
  output$UMAP_Facet_Ant_Choice <- renderUI({
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
  observeEvent(input$update_UMAP_facet, {
    def$choice_UMAP_facet_colour_by = input$checkBox_UMAP
    def$choice_UMAP_Facet_Ant_Choice = UMAP_grouping2()
    def$UMAP_facet_update_colourby = input$UMAP_facet_colour_by
    def$UMAP_facet_update_ant = input$UMAP_Facet_Ant_Choice
    def$UMAP_facet_update_text <- ""
  })

  # Logic for Update Reminder Text:
  observeEvent({
    input$UMAP_facet_colour_by
    input$UMAP_Facet_Ant_Choice
    input$checkBox_UMAP
  },
  {
    if (input$UMAP_facet_colour_by != def$UMAP_facet_update_colourby) {
      def$UMAP_facet_update_text <- "Press the update button."
    }
    else if ((!is.null(input$UMAP_Facet_Ant_Choice)) && (input$UMAP_Facet_Ant_Choice != def$UMAP_facet_update_ant)) {
      def$UMAP_facet_update_text <- "Press the update button."
    }
    else if (!same_elements(input$checkBox_UMAP, def$choice_UMAP_facet_colour_by)) {
      def$UMAP_facet_update_text <- "Press the update button."
    }
    else {
      def$UMAP_facet_update_text <- ""
    }
  })

  # Renders Reminder Text
  output$UMAP_facet_update_text <- renderText({ def$UMAP_facet_update_text })

  # Checbox Deselect All Button:
  observeEvent(input$deselectAll_UMAP, {
    updateCheckboxGroupInput(session, "checkBox_UMAP", selected = list())
    if (length(def$choice_UMAP_facet_colour_by) > 0) {def$UMAP_facet_update_text <- "Select at least one sample id."}
  })

  # CheckBoxGroup inputs: Limit to 10 Sample_IDs to select.
  observe({
    my_min <- 0
    my_max <- 10
    if(length(input$checkBox_UMAP) > my_max){
      updateCheckboxGroupInput(session, "checkBox_UMAP", selected = tail(input$checkBox_UMAP,my_max))
    }
  })
  # Plot Title
  output$UMAP_facet_Text <- renderText(paste0("UMAP: Coloured By ", UMAP_facet_Text(), ", Separated by Sample_id"))

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
      paste(paste0("UMAP: Coloured By ", UMAP_facet_Text(), ", Separated by Sample Ids"), input$UMAP_2_tag, sep=".")
    },
    content = function(file) {
      req(plot_UMAP_facet())
      ggsave(file, plot = plot_UMAP_facet(), device = input$UMAP_2_tag)
    }
  )

  # =================================================================================================
  # Page 4: Cluster Proportions
  # -------------------------------------------------------------------------------------------------
  Abundance_cluster <- reactive({
    daf$sample_id<-factor(daf$sample_id,levels = sampleID_sorted)
    plotAbundances(daf, k = cluster_var, by = "sample_id") +
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
      ggsave(file, plot = Abundance_cluster(), device = input$Abundance_cluster_tag, width = 2*nlevels(md$sample_id), height = 21, units = "cm")
    }
)})
