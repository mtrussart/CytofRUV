library(shiny)
library(shinydashboard)
library(shinyjs)
library(ggplot2)
library(flowCore)
library(purrr)
library(CATALYST)
library(CytofRUV)



source("server.R")
#source("server-diagnostic_plots.R")
#source("ui-diagnostic_plots.R")
source("ui.R")


# Allows for Collapsible Boxes
collapseBox <- "shinyjs.collapse=function(id){
$('#'+id).closest('.box').not('.collapsed-box')
.find('[data-widget=collapse]').click();}"

