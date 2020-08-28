## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(CytofRUV)
library(CATALYST)
library(flowCore)
library(ggplot2)
library(readxl)
library(ruv)
library(purrr)
library(FlowSOM)
library(SummarizedExperiment)
library(ConsensusClusterPlus)
library(SingleCellExperiment)
library(shiny)
library(shinyjs)
library(shinydashboard)
library(writexl)
library(ComplexHeatmap)
library(shinycssloaders)

## ----loading example dataset--------------------------------------------------
  output_dir="CytofRUV_output"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  wd_data=file.path(getwd(),output_dir)
  write.FCS(x=CytofRUV::A1,filename =file.path(wd_data,"A1.fcs"))
  write.FCS(x=CytofRUV::A2,filename = file.path(wd_data,"A2.fcs"))
  write.FCS(x=CytofRUV::Run3_A1,filename = file.path(wd_data,"Run3_A1.fcs"))
  write.FCS(x=CytofRUV::Run3_A2,filename = file.path(wd_data,"Run3_A2.fcs"))
  write_xlsx(x=CytofRUV::md,path = file.path(wd_data,"Metadata.xlsx"))
  write_xlsx(x=CytofRUV::panel,path = file.path(wd_data,"Panel.xlsx"))

## ----load and cluster data before CytofRUV normalisation----------------------
  ## Define parameters to load and cluster the data
  output_dir="CytofRUV_output"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  wd_data=file.path(getwd(),output_dir)
  metadata_filename="Metadata.xlsx"
  panel_filename="Panel.xlsx"
  seed=1234
  clusters_nb=20

  ## Loading the data
  data=load_data(wd_data,metadata_filename,panel_filename)

  ## Cluster the data
  data$daf=cluster_data(data$daf,seed,markers_to_use=data$lineage_markers,clusters_nb)


## ----ShinyApp on the raw data,eval=FALSE--------------------------------------
#    ## Define parameters to use for R-Shiny
#    daf=data$daf
#    md=data$md
#    seed=1234
#    # Number of cells displayed for the Distribution of protein expression plot
#    n_subset_marker_specific <- 10000
#  
#    # Running Dimension Reduction -> TSNE & UMAP
#    set.seed(seed)
#  
#    # Number of cells for tSNE plots marker specific
#    TSNE_subset <- 2000
#    print("Running TSNE")
#    daf <- CATALYST::runDR(daf, dr = "TSNE", cells = TSNE_subset)
#  
#    # Number of cells for UMAP plots marker specific
#    UMAP_subset <- 2000
#    print("Running UMAP")
#    daf <- runDR(daf, "UMAP", cells = UMAP_subset)
#  
#    ## Launch Shiny
#    # For a subset of the data, define the number of cells for diagnostic plots
#    n_subset <- 5000
#    sub_daf <- daf[, sample(ncol(daf), n_subset)]
#  
#    ## For the full dataset:
#    # sub_daf <- daf
#  
#    panel=data$panel
#  
#    CytofRUV::launch_Shiny()

## ----CytofRUV normalisation---------------------------------------------------
#dir_name_norm_data = get_directory()
dir_name_norm_data="CytofRUV_Norm_data_HC2_all_cl_20"
raw_data <- data.frame(sample = data$daf$sample_id, cluster=cluster_ids(data$daf,"meta20"), t(SummarizedExperiment::assay(data$daf)))
colnames(raw_data) <- gsub("^X", "",  colnames(raw_data))
rep_samples=list(c("HC2_B1","HC2_B2"))
cluster_list_rep_samples <- list(seq(1,20))
k_value <- 5
seed=1234

normalise_data(data=data,raw_data=raw_data,rep_samples=rep_samples, norm_clusters=cluster_list_rep_samples, k=k_value, num_clusters=clusters_nb,wd_data=wd_data,dir_norm_data=dir_name_norm_data)


## ----load and cluster data after CytofRUV normalisation-----------------------
  ## Define parameters to load and cluster the data
  wd_norm=file.path(wd_data,dir_name_norm_data)
  metadata_norm_filename="Norm_Metadata.xlsx"
  panel_norm_filename="Norm_Panel.xlsx"
  seed=1234
  clusters_nb=20

  ## Loading the data
  norm_data=load_data(wd_norm,metadata_norm_filename,panel_norm_filename,cofactor = NULL)

  ## Cluster the data
  norm_data$daf=cluster_data(norm_data$daf,seed,markers_to_use=norm_data$lineage_markers,clusters_nb)


## ----Shiny App on the normalised data,eval=FALSE------------------------------
#    ## Define parameters to use for R-Shiny
#    daf=norm_data$daf
#    md=norm_data$md
#    seed=1234
#    # Number of cells for diagnostic plots marker specific
#    n_subset_marker_specific <- 10000
#  
#    # Define type of markers
#    daf_type <- daf[SingleCellExperiment::rowData(daf)$marker_class=="type", ]
#    daf_state <- daf[SingleCellExperiment::rowData(daf)$marker_class=="state", ]
#    sub_daf_state <- daf_state[, sample(ncol(daf_state), n_subset_marker_specific)]
#    sub_daf_type <- daf_type[, sample(ncol(daf_type), n_subset_marker_specific)]
#    # Define batch
#    batch_ids <- is.factor(rep(md$batch, nrow(daf)))
#    sampleID_sorted <- md$sample_id[order(md$patient_id)]
#  
#    ## Running Dimension Reduction -> TSNE & UMAP
#    set.seed(seed)
#    # Number of cells for tSNE plots marker specific
#    TSNE_subset <- 2000
#    print("Running TSNE")
#    daf <- runDR(daf, "TSNE", cells = TSNE_subset)
#  
#    # Number of cells for UMAP plots marker specific
#    UMAP_subset <- 2500
#    print("Running UMAP")
#    daf <- runDR(daf, "UMAP", cells = UMAP_subset)
#  
#    # Launch Shiny
#    # For a subset of the data, define the number of cells for diagnostic plots
#    n_subset <- 5000
#    sub_daf <- daf[, sample(ncol(daf), n_subset)]
#  
#    # # For the full dataset:
#    # sub_daf <- daf
#  
#    panel=data$panel
#  
#    CytofRUV::launch_Shiny()

