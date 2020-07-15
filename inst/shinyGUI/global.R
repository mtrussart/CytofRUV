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

## All Code Below Gets The Scope of the Global Environment.

# This gets the largest "meta#" term from the list of cluster_codes
cluster_var = names(metadata(data$daf)$cluster_codes)[length(names(metadata(data$daf)$cluster_codes))]

# Allows for Collapsible Boxes
collapseBox <- "shinyjs.collapse=function(id){
$('#'+id).closest('.box').not('.collapsed-box')
.find('[data-widget=collapse]').click();}"

