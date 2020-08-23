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
#' @param cofactor Cofactor for asinh transformation, default is 5 and set to NULL for untransformed data
#'
#'
#' @return Datasets before normalisation
#' @export
#'


#load_data("/Users/trussart.m/WEHI/CytofRUV/CytofRUV/data/",metadata_filename="Metadata.xlsx",panel_filename="panel.xlsx")
load_data<- function(wd_data,metadata_filename,panel_filename,cofactor=5){
  if (!is.null(cofactor)){
    transform=TRUE
  }else{
    transform=FALSE
    cofactor=5
  }
  data=read_data(wd_data,metadata_filename,panel_filename,transform,cofact=cofactor)
  return(data)
}


read_data<- function(wd_data,metadata_filename,panel_filename,transform,cofact=5){

  ## Load the metadata file
  print("Reading MetaData")
  md <- readxl::read_excel(file.path(wd_data, metadata_filename))

  ## Check the metadata file
  # print("Checking md ColNames")
  checkColnames_md(md)
  # print("Checking md FileNames")
  checkFileNames_md(md)
  # print("Check md Values")
  checkAllValues_md(md)

  ## Load the fcs files
  setwd(wd_data)
  print("Reading fcs files")
  fcs_raw <- flowCore::read.flowSet(file.path(wd_data, md$file_name), transformation = FALSE,
                          truncate_max_range = FALSE)
  #print(flowCore::colnames(fcs_raw[[1]]))
  #setwd(wd_data)
  #print("Read Fcs Files.")

  ## Load the panel file
  panel <- readxl::read_excel(file.path(wd_data, panel_filename))
  panel_fcs <- pData(parameters(fcs_raw[[1]]))
  print(panel_fcs)
  print("Read Panel File")

  ## Check the panel file:
  print("Checking panel colNames")
  checkColnames_panel(panel)
  print("Checking panel values")
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

  daf <- prepData(fcs_raw, panel, md, md_cols = list(file="file_name", id="sample_id", factors = c("batch", "condition", "patient_id")),transform=transform,cofactor = cofact)
  if (isFALSE(transform)){
    assayNames(daf) = "exprs"
  }
  print("Completed creating DaFrame Objects")

  ## All data
  data=list(fcs_raw=fcs_raw,md=md,panel=panel,lineage_markers=lineage_markers,functional_markers=functional_markers,daf=daf)
  return(data)
}
#
# create_data=function (x, panel, md, features = NULL, cofactor = 5, panel_cols = list(channel = "fcs_colname",
#   antigen = "antigen", class = "marker_class"), md_cols = list(file = "file_name", id = "sample_id", factors = c("condition", "patient_id"))){
#   if (!is(panel, "data.frame"))
#     panel <- data.frame(panel, check.names = FALSE, stringsAsFactors = FALSE)
#   if (!is(md, "data.frame"))
#     md <- data.frame(md, check.names = FALSE, stringsAsFactors = FALSE)
#   stopifnot(is.list(panel_cols), is.list(md_cols), c("channel", "antigen") %in% names(panel_cols), c("file", "id", "factors") %in% names(md_cols))
#   if (!is.null(cofactor))
#     stopifnot(is.numeric(cofactor), length(cofactor) == 1,
#               cofactor > 0)
#   if (is(x, "flowSet")) {
#     fs <- x
#   }
#   else if (is.character(x)) {
#     stopifnot(dir.exists(x))
#     fcs <- list.files(x, ".fcs$", full.names = TRUE, ignore.case = TRUE)
#     if (length(fcs) < 2)
#       stop("The specified directory contains", " none or only a single FCS file.")
#     stopifnot(all(vapply(fcs, isFCSfile, logical(1))))
#     fs <- read.flowSet(fcs, transformation = FALSE, truncate_max_range = FALSE)
#   }
#   else {
#     stop("Invalid argument 'x'; should be either a flowSet",
#          " or a character string specifying the path to",
#          " a directory containing a set of FCS files.")
#   }
#   #stopifnot(panel[[panel_cols$channel]] %in% colnames(fs))
#   if (is.null(features)) {
#     features <- as.character(panel[[panel_cols$channel]])
#   }
#   else {
#     chs <- colnames(fs)
#     check1 <- is.logical(features) && length(features) ==
#       length(chs)
#     check2 <- is.integer(features) && all(features %in% seq_along(chs))
#     check3 <- all(features %in% chs)
#     if (!any(check1, check2, check3))
#       stop("Invalid argument 'features'. Should be either",
#            " a logial vector,\n  a numeric vector of indices, or",
#            " a character vector of column names.")
#   }
#   ids <- c(keyword(fs, "FILENAME"))
#   if (is.null(unlist(ids)))
#     ids <- c(flowCore::fsApply(fs, identifier))
#   stopifnot(all(ids %in% md[[md_cols$file]]))
#   fs <- fs[match(ids, md[[md_cols$file]])]
#   if (!is.null(cofactor))
#     fs <- flowCore::fsApply(fs, function(ff) {
#       exprs(ff) <- asinh(exprs(ff)/cofactor)
#       return(ff)
#     })
#   k <- c(md_cols$id, md_cols$factors)
#   md <- data.frame(md)[, k] %>% dplyr::mutate_all(factor) %>% dplyr::rename(sample_id = md_cols$id)
#   o <- order(md[[md_cols$factors[1]]])
#   md$sample_id <- factor(md$sample_id, levels = md$sample_id[o])
#   antigens <- panel[[panel_cols$antigen]]
#   antigens <- gsub("-", "_", antigens)
#   antigens <- gsub(":", ".", antigens)
#   fs <- fs[, features]
#   chs0 <- colnames(fs)
#   m1 <- match(panel[[panel_cols$channel]], chs0, nomatch = 0)
#   m2 <- match(chs0, panel[[panel_cols$channel]], nomatch = 0)
#   flowCore::colnames(fs)[m1] <- antigens[m2]
#   chs <- colnames(fs)
#   es <- matrix(flowCore::fsApply(fs, flowSet::exprs), byrow = TRUE, nrow = length(chs),
#                dimnames = list(chs, NULL))
#   md$n_cells <- as.numeric(flowCore::fsApply(fs, nrow))
#   valid_mcs <- c("type", "state", "none")
#   if (is.null(panel_cols$class)) {
#     mcs <- factor("none", levels = valid_mcs)
#   }
#   else {
#     mcs <- factor(panel[[panel_cols$class]], levels = valid_mcs)
#     mcs <- mcs[match(chs0, panel[[panel_cols$channel]])]
#     if (any(is.na(mcs)))
#       stop("Invalid marker classes detected.", " Valid classes are 'type', 'state', and 'none'.")
#   }
#   rd <- DataFrame(row.names = chs, channel_name = chs0, marker_name = chs,
#                   marker_class = mcs)
#   k <- setdiff(names(md), "n_cells")
#   cd <- DataFrame(lapply(md[k], function(u) {
#     v <- as.character(rep(u, md$n_cells))
#     factor(v, levels = levels(u))
#   }), row.names = NULL)
#   SingleCellExperiment(assays = list(exprs = es), rowData = rd,
#                        colData = cd, metadata = list(experiment_info = md, cofactor = cofactor))
# }

# ==============================================================================
# split cell indices by cell metadata factor(s)
#   - x:   a SCE with rows = cells, columns = features
#   - by:  colData columns specifying factor(s) to aggregate by
# ------------------------------------------------------------------------------
#' @importFrom data.table data.table
#' @importFrom SummarizedExperiment colData
#' @importFrom purrr map_depth
.split_cells <- function(x, by) {
  stopifnot(is.character(by), by %in% colnames(colData(x)))
  cd <- data.frame(colData(x))
  dt <- data.table(cd, i = seq_len(ncol(x)))
  dt_split <- split(dt, by = by, sorted = TRUE, flatten = FALSE)
  map_depth(dt_split, length(by), "i")
}

# ==============================================================================
# aggregation of single-cell to pseudobulk data;
# e.g., median expression by cluster- or cluster-sample
#   - x:   a SCE with rows = cells, columns = features
#   - by:  colData columns specifying factor(s) to aggregate by
#   - fun: aggregation function specifying the
#          summary statistic, e.g., sum, mean, median
# ------------------------------------------------------------------------------
#' @importFrom dplyr bind_rows
#' @importFrom Matrix rowMeans rowSums
#' @importFrom matrixStats rowMedians
#' @importFrom purrr map_depth
.agg <- function(x, by, fun = c("median", "mean", "sum")) {
  fun <- switch(match.arg(fun),
                median = rowMedians, mean = rowMeans, sum = rowSums)
  cs <- .split_cells(x, by)
  pb <- map_depth(cs, -1, function(i) {
    if (length(i) == 0) return(numeric(nrow(x)))
    fun(assay(x, "exprs")[, i, drop = FALSE])
  })
  map_depth(pb, -2, function(u) as.matrix(data.frame(
    u, row.names = rownames(x), check.names = FALSE)))
}

#' @importFrom ComplexHeatmap columnAnnotation rowAnnotation
#' @importFrom grid gpar
#' @importFrom methods is
#' @importFrom scales hue_pal
#' Function From CATALYST
.anno_factors <- function(df, type = c("row", "column")) {
  # check that all data.frame columns are factors
  stopifnot(is(df, "data.frame"))
  stopifnot(all(vapply(as.list(df), is.factor, logical(1))))
  # for ea. factor, extract levels & nb. of levels
  lvls <- lapply(as.list(df), levels)
  nlvls <- vapply(lvls, length, numeric(1))
  # cols <- pal.safe(parula, n = sum(nlvls), main = NULL)
  names(cols) <- unlist(lvls)
  cols <- split(cols, rep.int(seq_len(ncol(df)), nlvls))
  names(cols) <- names(df)
  HeatmapAnnotation(which = match.arg(type),
                    df = df, col = cols, gp = gpar(col = "white"))
}

get_directory <- function() {
  print("Please refer to the console")
  directory<-""
  while (nchar(directory) == 0) {
    directory <- readline("Please provide the directory name where the the normalized data will be saved: ")
    if (nchar(directory) == 0) {
      print("Error: Nothing was provided. Please enter the directory name")
    }
  }
  if (!dir.exists(directory)) {
    dir.create(directory)
    print(sprintf("The folder: %s was created.", directory))
    return(directory)
  }
  print(sprintf("%s was chosen as the directory.", directory))
  return(directory)
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
  print(colnames(flowCore::colnames(fcs_raw)))
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


