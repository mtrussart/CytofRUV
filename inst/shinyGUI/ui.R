#source("ui-diagnostic_plots.R", local = TRUE)

# =================================================================================================
# Renders Tab Panel
# -------------------------------------------------------------------------------------------------
ui_diagnostic_plots <- fluidPage(
  h1("Diagnostic Plots"),
  fluidRow(
    column(
      12,
      uiOutput("tab_plots")
    )
  )
)

# =================================================================================================
# Page 1: Median Intensities
# -------------------------------------------------------------------------------------------------
medIntensities <-fluidPage(
  fluidRow(
    column(
      5,
      ## Plot 1
      h2("Multi-dimensional scaling plot computed using median protein expression"),
      plotOutput(outputId="plotMDS", width="100%", height = "600px"),
      selectInput("choiceMDS", "Select to colour by:",
                  choices = list("condition", "batch")),  ##Colour By is limited to condition and batch!
      actionButton("mds", "update")
    ),
    column(
      7,
      ## Plot 2
      h2("Heatmap of the median protein expression"),
      plotOutput(outputId="plotDendogram", width="80%", height = "600px"),
      align = "center"
    )
  )
)

# =================================================================================================
# Page 2: Markers Distribution
# -------------------------------------------------------------------------------------------------
markerDistribution <-fluidPage(
  fluidRow(
    column(
      12,
      ## Plot 3
      h2("Distribution of protein expression"),
      plotOutput("exprsPlot", height="700px"),
    )
  ),
  fluidRow(
    column(4,
           h4("Select Parameters"),
           selectInput("exprs1", "Select your parameters:",
                       list("Condition" = "condition",
                            "Sample ID" = "sample_id")),
           uiOutput("exprs2"),
           uiOutput("exprs3"),
           actionButton("exprsPlot", "update")
    )
  )
)

# =================================================================================================
# Page 3: Clustering Results
# -------------------------------------------------------------------------------------------------
controls <-
  fluidPage(
    fluidRow(
      column(
        12,
        list(tags$div(align = 'left',
                      class = 'multicol',
                      checkboxGroupInput("checkBox_sampleIDs", "Select Up To 10 SampleIDs",
                                         choices = sampleID_sorted,
                                         selected = sampleID_sorted[1:10],
                                         inline = TRUE,
                                         width = "100%"))),
        actionButton("deselectAll", "Deselect All Options")
      )
    )
  )

checkBox_SampleIDs <- shinydashboard::box(
  title="Select Sample IDs",
  solidHeader=TRUE,
  status="warning",
  id="box_2",
  width=12,
  collapsible=TRUE,
  controls,
  inline = FALSE, choiceNames = NULL,
  choiceValues = NULL
)

clusteringResults <-fluidPage(
  fluidRow(
    column(
      8,
      ## Plot 1
      h2("Heatmap of the median protein expression per cluster"),
      plotOutput("cluster_heatmap", height="650px")
    )
  ),
  fluidRow(
    column(
      12,
      ## Plot 2
      hr(),
      h2(textOutput("textDR_1")),
      plotOutput("plotDR_1", width = "750px", height = "550px"),
      fluidRow(
        column(
          6,
          selectInput("DR1a", "Colour By:", list("Cluster"="meta20", "Antigen"="Antigen", "Batch"="batch"))
        ),
        column(
          6,
          uiOutput("antigen_choice1")
        )
      )
    )
  ),
  fluidRow(
    column(
      4,
      actionButton("plotDR_1", "update"),
      hr()
    )
  ),
  fluidRow(
    column(
      12,
      ## Plot 3
      h2(textOutput("textDR_facet")),
      plotOutput("plotDR_facet", width = "700px", height = "600px"),
      fluidRow(
        column(
          6,
          selectInput("DR2a", "Colour By:", list("Cluster"="meta20", "Antigen"="Antigen", "Batch"="batch"))
        ),
        column(
          6,
          uiOutput("antigen_choice2")
        )
      )
    )
  ),
  fluidRow(
    column(
      6,
      actionButton("plotDR_2", "update")
    )
  ),
  fluidRow(
    column(
      12,
      checkBox_SampleIDs
    )
  )
)
# =================================================================================================
# Page 4: Cluster Proportions
# -------------------------------------------------------------------------------------------------
clusterProportions <- fluidPage(
  fluidRow(
    column(
      12,
      align="center",
      h2("Cluster proportions across sample"),
      ## Plot 1
      plotOutput("Abundance_cluster", height="800px")
    )
  )
)

# =================================================================================================
# Define Design of Tab Panel
# -------------------------------------------------------------------------------------------------
tab_plots <- shinydashboard::tabBox(
  side = "left",
  # selected = "Median Intensities",
  width = 12,
  tabPanel(h4("Median Protein Expression"), medIntensities),
  tabPanel(h4("Protein Expression Distributions"), markerDistribution),
  tabPanel(h4("Clustering Results"), clusteringResults),
  tabPanel(h4("Cluster Proportions"), clusterProportions)
)


##############3
header <- dashboardHeader(
  title = a(href="https://www.wehi.edu.au/",
            img(src="https://www.wehi.edu.au/sites/default/files/WEHI_logo_2016_0.png",
                style="padding-top:5px;padding-right:5px;",
                height="70px")
  ),
  titleWidth = 250,
  # title = span(img(src="Logo1.jpg", width = 190)
  tags$li(class = "dropdown",
          tags$style(".main-header {max-height: 65px}"),
          tags$style(".main-header .logo {height: 75px}"),
          tags$style(
            HTML(".box {
                    margin: 10px;
                    }
                    .checkbox-inline {
                        margin-left: 10px;
                        margin-right: 10px;
                    }
                    .checkbox-inline+.checkbox-inline {
                        margin-left: 10px;
                        margin-right: 10px;
                    }
                    .multicol {
                        height: 150px;
                        -webkit-column-count: 4; /* Chrome, Safari, Opera */
                        -moz-column-count: 4;    /* Firefox */
                        column-count: 4;
                        -moz-column-fill: auto;
                        -column-fill: auto;
                     }")
          )
  )
)

sidebar <- dashboardSidebar(disable = TRUE)

body <- dashboardBody(
  fluidRow(
    tabsetPanel(id = "tabset",
                tabPanel(h4("Diagnostic Plots"), ui_diagnostic_plots)
    )
  ),
)

ui <- dashboardPage(header, sidebar, body, skin="black")


