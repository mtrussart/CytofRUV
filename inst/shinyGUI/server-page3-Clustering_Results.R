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
    paste(paste0("TSNE: Coloured By ", textDR_2(), ", Separated by Sample_id"), input$TSNE_facet_tag, sep=".")
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
    paste(paste0("UMAP: Coloured By ", UMAP_facet_Text(), ", Separated by Sample_id"), input$UMAP_2_tag, sep=".")
  },
  content = function(file) {
    req(plot_UMAP_facet())
    ggsave(file, plot = plot_UMAP_facet(), device = input$UMAP_2_tag)
  }
)
