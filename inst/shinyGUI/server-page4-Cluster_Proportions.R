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
)
