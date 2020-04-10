#' @rdname load_data
#' @title load_data
#'
#' is used to load the fcs files from all the samples in the study,
#' the metadata file containing the details of each sample
#' and a panel file containing the details of all proteins used in the study.
#'
#'
#' @param wd_data Path to the directory containing all raw fcs files, the metadata file
#' and the panel file
#' @param metadata_name Metadata filename containing the details of each sample
#' @param panel_filename Panel filename containing the details of each marker
#'
#'
#'
#' @return Datasets before normalisation
#' @export
#'

get_directory <- function() {
  print("Please refer to the console.")
  directory<-""
  while (nchar(directory) == 0) {
    directory <- readline("Name the directory you want the normalized files to be saved in: ")
    if (nchar(directory) == 0) {
      print("Error: Nothing was provided, please enter a working directory.")
    }
  }
  if (!dir.exists(directory)) {
    dir.create(directory)
    print(sprintf("The folder: %s was created as the directory.", directory))
    return(directory)
  }
  print(sprintf("%s was chosen as the directory.", directory))
  return(directory)
}


#load_data("/Users/trussart.m/WEHI/CytofRUV/CytofRUV/data/",metadata_filename="Metadata.xlsx",panel_filename="panel.xlsx")
load_data<- function(wd_data,metadata_filename,panel_filename){
  data=read_data(wd_data,metadata_filename,panel_filename)
  return(data)
}


read_data<- function(wd_data,metadata_filename,panel_filename){

  ## Load the metadata file
  print("Reading MetaData")
  md <- readxl::read_excel(file.path(wd_data, metadata_filename))

  ## Check the metadata file
  #print("Checking md ColNames")
  checkColnames_md(md)
  #print("Checking md FileNames")
  checkFileNames_md(md)
  #print("Check md Values")
  checkAllValues_md(md)

  ## Load the fcs files
  setwd(wd_data)
  print("Reading fcs files")
  fcs_raw <- flowCore::read.flowSet(file.path(md$file_name), transformation = FALSE,
                          truncate_max_range = FALSE)
  #print(flowCore::colnames(fcs_raw[[1]]))
  #setwd(wd_data)
  #print("Read Fcs Files.")

  ## Load the panel file
  panel <- readxl::read_excel(file.path(wd_data, panel_filename))
  panel_fcs <- pData(parameters(fcs_raw[[1]]))
  #print(panel_fcs)
  print("Read Panel File")

  ## Check the panel file:
  #print("Checking panel colNames")
  checkColnames_panel(panel)
  #print("Checking panel values")
  checkAllValues_panel(panel, fcs_raw)

  # ===========
  # Format Data
  # -----------
  print("Formatting Data")
  ## Make sure condition variables are factors with the right levels
  md$condition  <- factor(md$condition)
  md$sample_id  <- factor(md$sample_id, levels = md$sample_id[order(md$condition)])
  md$patient_id <- factor(md$patient_id)
  md$batch      <- factor(md$batch)

  panel$marker_class <- factor(panel$marker_class)

  # Replace problematic characters
  panel_fcs$desc <- gsub("-", "_", panel_fcs$desc)
  panel$antigen <- gsub("-", "_", panel$antigen)

  # Can be used in features of cluster
  # Lineage markers
  (lineage_markers <- panel$antigen[panel$marker_class == "type"])
  (lineage_markers_fullname <- panel$fcs_colname[panel$marker_class == "type"])
  # Functional markers
  (functional_markers <- panel$antigen[panel$marker_class == "state"])
  (functional_markers_fullname <- panel$fcs_colname[panel$marker_class == "state"])

  daf <- prepData(fcs_raw, panel, md, md_cols = list(file="file_name", id="sample_id", factors = c("batch", "condition", "patient_id")))
  print("Completed creating DaFrame Objects")

  ## All data
  data=list(fcs_raw=fcs_raw,md=md,panel=panel,lineage_markers=lineage_markers,functional_markers=functional_markers,daf=daf)
  return(data)
}


transform_data = function(fcs_raw,lineage_markers_fullname,functional_markers_fullname) {
  ## arcsinh transformation and column subsetting
  print("Arcsinh Transformation")
  fcs <- fsApply(fcs_raw, function(x, cofactor = 5){
    expr <- exprs(x)
    expr <- expr[, c(lineage_markers_fullname, functional_markers_fullname)]
    expr <- asinh(expr/ cofactor)
    exprs(x) <- expr
    x
  })
  return (fcs)
  #print("Arcsinh Transformation Complete!")
}

# ========================================================================
# Returns sample_IDs related to a given patient, found through patient_ID.
# ------------------------------------------------------------------------
patient_ids <- function(patient_id, dframe){
  return(levels(factor(sample_ids(dframe)[grepl(patient_id,sample_ids(dframe))])))
}


# ===========================================
# Makes sure file_names in the metadata files have the correct extension "*.fcs"
# -------------------------------------------
checkFileNames_md <- function(md) {
  cond <- unlist(lapply(md$file_name, grepl, pattern = ".fcs"))
  if(!all(cond)) {
    stop({
      print("Error: The following files are not .fcs files: ")
      print(md$file_name[which(cond %in% FALSE)])
    })
  }
}

# ==============================================================================
# Helper Function used in checkColnames_md and checkColnames_panel to findf if a
# string is part of another list of strings.
# ------------------------------------------------------------------------------
colname_IN <- function(colName, existing_colNames) {
  return(colName %in% existing_colNames)
}

# =====================================================================
# Check if the mandatory column names are present in the metadata file.
# ---------------------------------------------------------------------
checkColnames_md <- function(md) {
  colnames_required <- c("file_name", "sample_id", "condition", "patient_id", "batch")
  cond <- unlist(lapply(colnames_required, colname_IN, existing_colNames=colnames(md)))
  if(!all(cond)) {
    stop({
      print("Error: The following column names were missing from the metadata file: ")
      print(colnames_required[which(cond %in% FALSE)])
    })
  }
}

# ===================================================================
# Check if there are any empty cells in the metadata file.
# -------------------------------------------------------------------
checkAllValues_md <- function(md) {
  if(any(unlist(purrr::flatten(lapply(md, is.na))))) {
    stop(sprintf("Error: There are empty values in the metadata file."))
  }
}

# ==================================================================
# Check if the mandatory column names are present in the panel file.
# ------------------------------------------------------------------
checkColnames_panel <- function(panel) {
  colnames_required <- c("fcs_colname", "antigen", "marker_class")
  cond <- unlist(lapply(colnames_required, colname_IN, existing_colNames=colnames(panel)))
  if(!all(cond)) {
    stop({
      print("Error: The following column names were missing from the panel file: ")
      print(colnames_required[which(cond %in% FALSE)])
    })
  }
}

# ==========================================================================================================
# Checks if there are colnames (elements and antigens) in the panel file that do not exist in the
# FCS files.
# ----------------------------------------------------------------------------------------------------------
checkAllValues_panel <- function(panel, fcs_raw) {
  print(panel$fcs_colname)
  print(colnames(fcs_raw))
  cond1 <- unlist(lapply(panel$fcs_colname,colname_IN, existing_colNames=flowCore::colnames(fcs_raw)))
  if(!all(cond1)) {
    stop({
      print("Error: The following colnames from the panel file were not found in the .fcs files column data.")
      print(panel$fcs_colname[which(cond1 %in% FALSE)])
    })
  }
  if(length(panel$fcs_colname)!=length(unique(panel$fcs_colname))) {
    stop(print("Error: There are repeated fcs_colnames in the panel file."))
  }
  if(length(panel$antigen)!=length(unique(panel$antigen))) {
    stop(print("Error: There are repeated antigens in the panel file."))
  }
}


