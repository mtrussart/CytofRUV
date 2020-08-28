# =================================================================================================
# Options for loaders -> This requires package - shinycssloaders to be added to build package.
# -------------------------------------------------------------------------------------------------
options(spinner.color="#0275D8", spinner.color.background="#ffffff", spinner.size=2)

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
      ## Need to incl library("shinycssloaders") into the package
      withSpinner(plotOutput(outputId="plotMDS", height = "600px"), type = 2),
      fluidRow(
        column(6,
               h5(strong("Select the parameter to color by:")),
               selectInput("choiceMDS", NULL, choices = list("condition", "batch")),
               h5(strong("Press after updating parameter:")),
               column(6, uiOutput("MDS_updateReminder")),
               column(6, actionButton("mds", "update"))
        ),
        column(6,
               h5(strong("Select the file type and Download Plot:")),
               radioButtons("mds_tag", NULL, choices = list("pdf", "png")),
               downloadButton(outputId = "download_mds", label = "Download Plot")
        )
      ),
    ),
    column(
      7,
      ## Plot 2
      h2("Heatmap of the median protein expression"),
      withSpinner(plotOutput(outputId="plotDendogram", width="80%", height = "600px"), type = 2),
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
      uiOutput("exprsPlot"),
    )
  ),
  fluidRow(
    column(5,
           selectInput("exprs1", "Select the parameter to color by:", list("Condition" = "condition", "Sample ID" = "sample_id")),
           uiOutput("exprs2"),
           uiOutput("exprs3"),
           h5(strong("Press after updating parameter:")),
           column(3,uiOutput("Exprs_update_text")),
           column(3,actionButton("exprsPlot", "update"))
    ),
    column(4,
           h5(strong("Select the file type and Download Plot:")),
           radioButtons("exprsPlot_tag", NULL, choices = list("pdf", "png")),
           downloadButton(outputId = "download_exprsPlot", label = "Download Plot")
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
                    checkboxGroupInput("checkBox_TSNE", NULL,
                                       choices = sampleID_sorted,
                                       selected = sampleID_sorted[1:10],
                                       inline = FALSE,
                                       width = "100%"))),
      actionButton("deselectAll_TSNE", "Deselect All Options")
    )
  )
)

checkBox_TSNE <- shinydashboard::box(
  title="Select Up To 10 Sample IDs",
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
                    checkboxGroupInput("checkBox_UMAP", NULL,
                                       choices = sampleID_sorted,
                                       selected = sampleID_sorted[1:10],
                                       inline = FALSE,
                                       width = "100%"))),
      actionButton("deselectAll_UMAP", "Deselect All Options")
    )
  )
)

checkBox_UMAP <- shinydashboard::box(
  title="Select Up To 10 Sample IDs",
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
      withSpinner(plotOutput("cluster_heatmap", height="650px"), type = 2),
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
      withSpinner(plotOutput("plot_TSNE1", height = "550px"), type = 2),
      fluidRow(
        column(6,
               selectInput("TSNE_Colour_By1", "Select the parameter to color by:", list("Cluster"=cluster_var, "Antigen"="Antigen", "Batch"="batch")),
               uiOutput("TSNE_Ant_Choice1"),
               h5(strong("Press after updating parameter:")),
               column(6,uiOutput("TSNE_update_text")),
               column(6, actionButton("update_TSNE1", "update"))
        ),
        column(6,
               h5(strong("Select the file type and Download Plot:")),
               radioButtons("TSNE_1_tag", NULL, choices = list("pdf", "png")),
               downloadButton(outputId = "download_TSNE_1", label = "Download Plot")
        )
      )
    ),
    column(
      6,
      ## Plot 3
      hr(),
      h2(textOutput("Umap_text_1")),
      withSpinner(plotOutput("plot_UMAP1", height = "550px"), type = 2),
      fluidRow(
        column(6,
               selectInput("Umap_Colour_By1", "Select the parameter to color by:", list("Cluster"=cluster_var, "Antigen"="Antigen", "Batch"="batch")),
               uiOutput("UMAP_Ant_Choice1"),
               h5(strong("Press after updating parameter:")),
               column(6,uiOutput("UMAP_update_text")),
               column(6,div(style="margin-top:0px; margin-bottom:0px;", actionButton("update_UMAP_1", "update")))
        ),
        column(6,
               h5(strong("Select the file type and Download Plot:")),
               radioButtons("Umap_1_tag", NULL, choices = list("pdf", "png")),
               downloadButton(outputId = "download_Umap_1", label = "Download Plot")
        )
      )
    )
  ),
  fluidRow(
    column(
      6,
      hr(),
      ## Plot 4
      hr(),
      h2(textOutput("TSNE_facet_Text")),
      withSpinner(plotOutput("plotTSNE_facet", height = "600px"), type = 2),
      fluidRow(
        column(6,
               selectInput("TSNE_facet_colourBy", "Select the parameter to color by:", list("Cluster"=cluster_var, "Antigen"="Antigen", "Batch"="batch")),
               uiOutput("TSNE_Facet_Ant_Choice"),
               h5(strong("Press after updating parameter:")),
               column(6,uiOutput("TSNE_facet_update_text")),
               column(6,div(style="margin-top:0px; margin-bottom:0px;", actionButton("update_TSNE_facet", "update")))
        ),
        column(6,
               h5(strong("Select the file type and Download Plot:")),
               radioButtons("TSNE_facet_tag", NULL, choices = list("pdf", "png")),
               downloadButton(outputId = "download_TSNE_facet", label = "Download Plot")
        )
      )
    ),
    column(
      6,
      hr(),
      ## Plot 5
      hr(),
      h2(textOutput("UMAP_facet_Text")),
      withSpinner(plotOutput("plot_UMAP_facet", height = "600px"), type = 2),
      fluidRow(
        column(6,
               selectInput("UMAP_facet_colour_by", "Select the parameter to color by:", list("Cluster"=cluster_var, "Antigen"="Antigen", "Batch"="batch")),
               uiOutput("UMAP_Facet_Ant_Choice"),
               h5(strong("Press after updating parameter:")),
               column(6,uiOutput("UMAP_facet_update_text")),
               column(6,div(style="margin-top:0px; margin-bottom:0px;", actionButton("update_UMAP_facet", "update")))
        ),
        column(6,
               h5(strong("Select the file type and Download Plot:")),
               radioButtons("UMAP_2_tag", NULL, choices = list("pdf", "png")),
               downloadButton(outputId = "download_UMAP_facet", label = "Download Plot")
        )
      )
    )
  ),
  fluidRow(
    column(4, checkBox_TSNE),
    column(4, offset = 2, checkBox_UMAP)
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
      withSpinner(plotOutput("Abundance_cluster", height="800px"), type = 2),
      h5(strong("Select the file type and Download Plot:")),
      fluidRow(
        column(2, radioButtons("Abundance_cluster_tag", NULL, choices = list("pdf", "png"))),
        column(2, downloadButton(outputId = "download_Abundance_cluster", label = "Download Plot"))
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

header <- dashboardHeader(
  title = a(href="https://www.wehi.edu.au/",
            img(src="https://www.wehi.edu.au/sites/default/files/WEHI_logo_2016_0.png",
                style="padding-top:5px;padding-right:5px;",
                height="70px")
  ),
  tags$li(class = "dropdown",
          tags$style(".main-header {max-height: 65px}"),
          tags$style(".main-header .logo {height: 75px}"),
          # CheckBox inputs
          tags$style(
            HTML(".box {
                    margin-top: 10px;
                  }
                  .checkbox-inline {
                    margin-left: 10px;
                    margin-right: 10px;
                  }
                  .multicol {
                    height: 150px;
                    overflow-y: auto;
                   }")),
          # Blue Update Text for Plots.
          tags$style(
            "#MDS_updateReminder{color: blue; font-size: 13px;}
             #Exprs_update_text{color: blue; font-size: 13px;}
             #TSNE_update_text{color: blue; font-size: 13px;}
             #UMAP_update_text{color: blue; font-size: 13px;}
             #TSNE_facet_update_text{color: blue; font-size: 13px;}
             #UMAP_facet_update_text{color: blue; font-size: 13px;}"
          ))
)

sidebar <- dashboardSidebar(disable = TRUE)

body <- dashboardBody(
  tags$head(tags$script(HTML(freezeButtonsScript))),
  fluidRow(
    tabsetPanel(id = "tabset", tabPanel(h4("Diagnostic Plots"), ui_diagnostic_plots))
  )
)

ui <- dashboardPage(header, sidebar, body, skin="black")
