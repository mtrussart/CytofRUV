# =================================================================================================
# Page 1: Median Intensities
# -------------------------------------------------------------------------------------------------
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
  plotMDS(sub_daf, color_by = def$choiceMDS) +
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
