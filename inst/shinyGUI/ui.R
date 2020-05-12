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
      actionButton("mds", "update"),
      downloadButton(outputId = "download_mds", label = "Download Plot"),
      radioButtons("mds_tag", "Select the file type:", choices = list("pdf", "png"))
    ),
    column(
      7,
      ## Plot 2
      h2("Heatmap of the median protein expression"),
      plotOutput(outputId="plotDendogram", width="80%", height = "600px"),
      downloadButton(outputId = "download_dendogram", label = "Download Plot"),
      radioButtons("dendogram_tag", "Select the file type:", choices = list("pdf", "png"))
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
      ## Page 2 plot 1
      h2("Distribution of protein expression"),
      plotOutput("exprsPlot", height="700px")
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
           actionButton("exprsPlot", "update"),
           downloadButton(outputId = "download_exprsPlot", label = "Download Plot"),
           radioButtons("exprsPlot_tag", "Select the file type:", choices = list("pdf", "png"))
    )
  )
)

# =================================================================================================
# Page 3: Clustering Results
# -------------------------------------------------------------------------------------------------
TNSE_check_box <-fluidPage(
  fluidRow(
    column(
      12,
      list(tags$div(align = 'left',
                    class = 'multicol',
                    checkboxGroupInput("checkBox_TSNE", "Select Up To 10 SampleIDs",
                                       choices = sampleID_sorted,
                                       selected = sampleID_sorted[1:10],
                                       inline = TRUE,
                                       width = "100%"))),
      actionButton("deselectAll_TSNE", "Deselect All Options")
    )
  )
)

checkBox_TSNE <- shinydashboard::box(
  title="Select Sample IDs",
  solidHeader=TRUE,
  status="warning",
  id="box_2",
  width=12,
  collapsible=TRUE,
  TNSE_check_box,
  inline = FALSE, choiceNames = NULL,
  choiceValues = NULL
)

UMAP_check_box <- fluidPage(
  fluidRow(
    column(
      12,
      list(tags$div(align = 'left',
                    class = 'multicol',
                    checkboxGroupInput("checkBox_UMAP", "Select Up To 10 SampleIDs",
                                       choices = sampleID_sorted,
                                       selected = sampleID_sorted[1:10],
                                       inline = TRUE,
                                       width = "100%"))),
      actionButton("deselectAll_UMAP", "Deselect All Options")
    )
  )
)

checkBox_UMAP <- shinydashboard::box(
  title="Select Sample IDs",
  solidHeader=TRUE,
  status="warning",
  id="box_3",
  width=12,
  collapsible=TRUE,
  UMAP_check_box,
  inline = FALSE, choiceNames = NULL,
  choiceValues = NULL
)

clusteringResults <-fluidPage(
  fluidRow(
    column(
      8,
      ## Plot 1
      h2("Heatmap of the median protein expression per cluster"),
      plotOutput("cluster_heatmap", height="650px"),
      downloadButton(outputId = "download_cluster_Heatmap", label = "Download Plot"),
      radioButtons("cluster_Heatmap_tag", "Select the file type:", choices = list("pdf", "png")),
      hr()
    )
  ),
  fluidRow(
    column(
      6,
      ## Plot 2
      hr(),
      h2(textOutput("TSNE_TEXT1")),
      plotOutput("plot_TSNE1", width = "750px", height = "550px"),
      fluidRow(
        column(
          6,
          selectInput("TSNE_Colour_By1", "Colour By:", list("Cluster"="meta20", "Antigen"="Antigen", "Batch"="batch"))
        ),
        column(
          6,
          uiOutput("TSNE_Ant_Choice1")
        )
      )
    ),
    column(
      6,
      ## Plot 3
      hr(),
      h2(textOutput("Umap_text_1")),
      plotOutput("plot_UMAP1", width = "750px", height = "550px"),
      fluidRow(
        column(
          6,
          selectInput("Umap_Colour_By1", "Colour By:", list("Cluster"="meta20", "Antigen"="Antigen", "Batch"="batch"))
        ),
        column(
          6,
          uiOutput("UMAP_Ant_Choice1")
        )
      )
    )
  ),
  fluidRow(
    column(
      6,
      actionButton("update_TSNE1", "update"),
      downloadButton(outputId = "download_TSNE_1", label = "Download Plot"),
      radioButtons("TSNE_1_tag", "Select the file type:", choices = list("pdf", "png")),
      hr()
    ),
    column(
      6,
      actionButton("update_UMAP_1", "update"),
      downloadButton(outputId = "download_Umap_1", label = "Download Plot"),
      radioButtons("Umap_1_tag", "Select the file type:", choices = list("pdf", "png")),
      hr()
    )
  ),
  fluidRow(
    column(
      6,
      ## Plot 4
      h2(textOutput("TSNE_facet_Text")),
      plotOutput("plotTSNE_facet", width = "700px", height = "600px"),
      fluidRow(
        column(
          6,
          selectInput("TSNE_facet_colourBy", "Colour By:", list("Cluster"="meta20", "Antigen"="Antigen", "Batch"="batch"))
        ),
        column(
          6,
          uiOutput("TSNE_Facet_Ant_Choice")
        )
      )
    ),
    column(
      6,
      ## Plot 5
      h2(textOutput("UMAP_facet_Text")),
      plotOutput("plot_UMAP_facet", width = "700px", height = "600px"),
      fluidRow(
        column(
          6,
          selectInput("UMAP_facet_colour_by", "Colour By:", list("Cluster"="meta20", "Antigen"="Antigen", "Batch"="batch"))
        ),
        column(
          6,
          uiOutput("UMAP_Facet_Ant_Choice")
        )
      )
    )
  ),
  fluidRow(
    column(
      6,
      actionButton("update_TSNE_facet", "update"),
      downloadButton(outputId = "download_TSNE_facet", label = "Download Plot"),
      radioButtons("TSNE_facet_tag", "Select the file type:", choices = list("pdf", "png")),
    ),
    column(
      6,
      actionButton("update_UMAP_facet", "update"),
      downloadButton(outputId = "download_UMAP_facet", label = "Download Plot"),
      radioButtons("UMAP_2_tag", "Select the file type:", choices = list("pdf", "png")),
    )
  ),
  fluidRow(
    column(
      6,
      checkBox_TSNE
    ),
    column(
      6,
      checkBox_UMAP
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
      h2("Cluster proportions across samples"),
      ## Plot 1
      plotOutput("Abundance_cluster", height="800px"),
      downloadButton(outputId = "download_Abundance_cluster", label = "Download Plot"),
      radioButtons("Abundance_cluster_tag", "Select the file type:", choices = list("pdf", "png")),
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


