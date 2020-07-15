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
    multiplier = 0
    if (nlevels(md$sample_id)%%5 == 0) {
      multiplier = 5
    } else {
      multiplier = nlevels(md$sample_id)%%5
    }
    ggsave(file, plot = exprsPlot(), device = input$exprsPlot_tag, width = (6.5 * multiplier) + 6, height = 18, units = "cm")
  }
)
