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
medIntensities <- fluidPage(
  fluidRow(
    column(
      5,
      ## Plot 1
      h2("Multi-dimensional scaling plot computed using median protein expression"),
      plotOutput(outputId="plotMDS", width="100%", height = "600px"),
      h5(strong("Select the parameter to color by and press update")),
      fluidRow(
        column(3, selectInput("choiceMDS", NULL, choices = list("condition", "batch"))),
        column(3, actionButton("mds", "update"))
      ),
      h5(strong("Select the file type and Download Plot:")),
      fluidRow(
        column(3, radioButtons("mds_tag", NULL, choices = list("pdf", "png"))),
        column(3, downloadButton(outputId = "download_mds", label = "Download Plot"))
      )
    ),
    column(
      7,
      ## Plot 2
      h2("Heatmap of the median protein expression"),
      plotOutput(outputId="plotDendogram", width="80%", height = "600px"),
      h5(strong("Select the file type and Download Plot:")),
      fluidRow(
        column(2, radioButtons("dendogram_tag", NULL, choices = list("pdf", "png"))),
        column(2, downloadButton(outputId = "download_dendogram", label = "Download Plot"))

      )
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
    column(
      3,
      h5(strong("Select the parameter to color by and press update")),
      fluidRow(
         column(
           6,
           selectInput("exprs1", NULL, list("Condition" = "condition", "Sample ID" = "sample_id")),
           uiOutput("exprs2"),
           uiOutput("exprs3")
         ),
         column(
           6,
           actionButton("exprsPlot", "update")
         )
      )
    ),
    column(
      3,
      h5(strong("Select the file type and Download Plot:")),
      fluidRow(
        column(6, radioButtons("exprsPlot_tag", NULL, choices = list("pdf", "png"))),
        column(6, downloadButton(outputId = "download_exprsPlot", label = "Download Plot"))
      )
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
      h5(strong("Select the file type and Download Plot:")),
      fluidRow(
        column(3, radioButtons("cluster_Heatmap_tag", NULL, choices = list("pdf", "png"))),
        column(3, downloadButton(outputId = "download_cluster_Heatmap", label = "Download Plot"))
      )
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
        column(3, selectInput("TSNE_Colour_By1", "Colour By:", list("Cluster"="meta20", "Antigen"="Antigen", "Batch"="batch"))),
        column(3, uiOutput("TSNE_Ant_Choice1")),
        column(4,
               h5(strong("Press after updating parameter:")),
               div(style="margin-top:0px; margin-bottom:0px;", actionButton("update_TSNE1", "update"))
        )
      ),
      fluidRow(h5(strong("Select the file type and Download Plot:"))),
      fluidRow(
        column(3, radioButtons("TSNE_1_tag", NULL, choices = list("pdf", "png"))),
        column(3, downloadButton(outputId = "download_TSNE_1", label = "Download Plot"))
      )
    ),
    column(
      6,
      ## Plot 3
      hr(),
      h2(textOutput("Umap_text_1")),
      plotOutput("plot_UMAP1", width = "750px", height = "550px"),
      fluidRow(
        column(3, selectInput("Umap_Colour_By1", "Colour By:", list("Cluster"="meta20", "Antigen"="Antigen", "Batch"="batch"))),
        column(3, uiOutput("UMAP_Ant_Choice1")),
        column(4,
               h5(strong("Press after updating parameter:")),
               actionButton("update_UMAP_1", "update"))
      ),
      fluidRow(h5(strong("Select the file type and Download Plot:"))),
      fluidRow(
        column(3, radioButtons("Umap_1_tag", NULL, choices = list("pdf", "png"))),
        column(3, downloadButton(outputId = "download_Umap_1", label = "Download Plot"))
      )
    )
  ),
  fluidRow(
    column(
      6,
      ## Plot 4
      h2(textOutput("TSNE_facet_Text")),
      plotOutput("plotTSNE_facet", width = "700px", height = "600px"),
      fluidRow(
        column(3, selectInput("TSNE_facet_colourBy", "Colour By:", list("Cluster"="meta20", "Antigen"="Antigen", "Batch"="batch"))),
        column(3, uiOutput("TSNE_Facet_Ant_Choice")),
        column(4,
               h5(strong("Press after updating parameter:")),
               actionButton("update_TSNE_facet", "update"))
      ),
      fluidRow(h5(strong("Select the file type and Download Plot:"))),
      fluidRow(
        column(3, radioButtons("TSNE_facet_tag", NULL, choices = list("pdf", "png"))),
        column(3, downloadButton(outputId = "download_TSNE_facet", label = "Download Plot"))

      )
    ),
    column(
      6,
      ## Plot 5
      h2(textOutput("UMAP_facet_Text")),
      plotOutput("plot_UMAP_facet", width = "700px", height = "600px"),
      fluidRow(
        column(3, selectInput("UMAP_facet_colour_by", "Colour By:", list("Cluster"="meta20", "Antigen"="Antigen", "Batch"="batch"))),
        column(3, uiOutput("UMAP_Facet_Ant_Choice")),
        column(4,
               h5(strong("Press after updating parameter:")),
               actionButton("update_UMAP_facet", "update"))
      ),
      fluidRow(h5(strong("Select the file type and Download Plot:"))),
      fluidRow(
        column(3, radioButtons("UMAP_2_tag", NULL, choices = list("pdf", "png"))),
        column(3, downloadButton(outputId = "download_UMAP_facet", label = "Download Plot"))

      )
    )
  ),
  fluidRow(
    column(6, checkBox_TSNE),
    column(6, checkBox_UMAP)
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
      h5(strong("Select the file type and Download Plot:")),
      fluidRow(
        column(2, radioButtons("Abundance_cluster_tag", NULL, choices = list("pdf", "png")),
        column(2, downloadButton(outputId = "download_Abundance_cluster", label = "Download Plot"))
        )
      )
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
    tabsetPanel(id = "tabset", tabPanel(h4("Diagnostic Plots"), ui_diagnostic_plots))
  )
)

ui <- dashboardPage(header, sidebar, body, skin="black")


