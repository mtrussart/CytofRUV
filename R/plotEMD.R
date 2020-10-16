#' @rdname plotEMD
#' @title plotEMD
#'
#' Insert description of plotEMD here
#'
#' @param daf a \code{\link[SingleCellExperiment]{SingleCellExperiment}}, this was generated after running examples present
#'   in the vignette.
#' @param raw_dir String: Directory to file containing raw data.
#' @param method_names Vector: Strings that describe method names.
#'   e.g. c(cytonorm, batchadj, cytofruv)
#' @param data_dir Vector: Strings of directories that contain appropriate data.
#'   IMPORTANT: Order of directory strings should correlate with order of method names in method_names.
#'   Note the directory for cytofRUV processed data is a vector of strings.
#'   e.g. c("../Cytofnorm_Cytof_Package/EMD_metric_comp_norm_all_samples_all_cells.rds",
#'          "../BatchAdjust/BatchAdjust_95p_and_DG23_new_trial/EMD_metric_comp_norm_all_samples_all_cells.rds"
#'          c("CytofRUV_Norm_data_CLL2_HC1_all_cl_20_k5/EMD_metric_comp_norm_all_samples_all_cells.rds",
#'            "CytofRUV_Norm_data_CLL2_HC1_all_cl_20_k15/EMD_metric_comp_norm_all_samples_all_cells.rds",
#'            "CytofRUV_Norm_data_CLL2_HC1_all_cl_20_k20/EMD_metric_comp_norm_all_samples_all_cells.rds"
#'          )
#'          )
#' @param k Vector: cluster numbers to be inspected, default val: c(5,10,15)
#' @param filter_patients Vector: relevant patient IDs  e.g. val = c("LL1_B1","LL2_B1","LL3_B1") that is,
#'   display patient IDs that contain: "LL1_B1","LL2_B1","LL3_B1".
#' @param x_axis_name String: Will be shared by all x-axis labels will share in common. E.g. "CLL", resulting
#'   x-axis items = "CLL1, CLL2, CLL3"
#' @param cols Vector: Strings of colours for plot annotations. Default val:
#'   c("raw"="black","batchadj"="indianred4","cytonorm"="darkmagenta","cytofruv_k1"="#084081","cytofruv_k3"="#0868AC","cytofruv_k5"="#2B8CBE","cytofruv_k7"="#4EB3D3","cytofruv_k10"="#7BCCC4","cytofruv_k12"= "#A8DDB5","cytofruv_k15"="darkgreen","cytofruv_k20"="seagreen","cytofruv_k25"="yellowgreen")
#'
#' @export

# cytofruv_k5, _k10 ...
# K is an optional parameter
plotEMD <- function(daf, raw_dir, method_names, data_dir, k=NULL, filter_patients, x_axis_name,
                    cols=c("raw"="black","batchadj"="indianred4","cytonorm"="darkmagenta","cytofruv_k1"="#084081","cytofruv_k3"="#0868AC","cytofruv_k5"="#2B8CBE","cytofruv_k7"="#4EB3D3","cytofruv_k10"="#7BCCC4","cytofruv_k12"= "#A8DDB5","cytofruv_k15"="darkgreen","cytofruv_k20"="seagreen","cytofruv_k25"="yellowgreen")
) {

  ## ======= Function Parameter Assertions! =======
  if(!file.exists(raw_dir)) {stop("The directory that was passed and meant to contain raw data does not exist!")}

  # Ensure num method names & k match num directories for data given.
  if(length(method_names)!=length(data_dir)) {stop(paste("Number of method names does not match number of directory strings provided", length(method_names),length(data_dir)))}

  # check that all directories provided actually exist.
  for (i in (1:length(method_names))) {
    if (method_names[i]=="cytofruv") {
      # Ensure num cytofruv files & k match.
      if(length(data_dir[[i]])!=length(k)) {stop(paste("Number of cytofRUV files does not match number of K values provided", length(data_dir[[i]]),length(k)))}

      for (j in (1:length(k))) {
        if(!file.exists(data_dir[[i]][[j]])) {stop(paste("The directory", data_dir[[i]][[j]], "does not exist!"))}
      }
    } else {
      if(!file.exists(data_dir[[i]])) {stop(paste("The directory", data_dir[[i]], "does not exist!"))}
    }
  }

  # Ensure patient id is in daf
  # for (i in (1:length(filter_patients))) {
  #   if (!filter_patients[[i]]%in%daf$patient_id) {stop(paste("The patient to filter by does not exist in the current dataset."))}
  # }


  # TODO: ASSERT FOR INPUTS!
    # i.e. CYTOFRUV is required as a method - people can
    # method names in method_names and cols have to match!!!
    # cytofruv_k# is the format!
    # raw is req
    # TODO: name given (patients) should be part of daf (add daf as a new param)

  ## ======= Main Function ======= ##

  EMD=list()

  ## Reading the raw file EMD
  EMD_raw_file=readRDS(raw_dir)
  EMD_raw=data.frame(patient=row.names(EMD_raw_file),EMD_raw_file)
  colnames(EMD_raw) <- gsub("^X", "",  colnames(EMD_raw))
  EMD[["raw"]]=EMD_raw

  cytof_methods <- vector()
  for (i in (1:length(method_names)) ){
    if (method_names[i]=="cytofruv") {
      ## Reading the CytofRUV files EMD
      for (j in (1:length(k)) ){
        # TODO: Marie will decide instead of taking k values, read all dirs inside the cytofRUV directory.
        # TODO: Users have to pass in files individually with full names
        # i.e. "CytofRUV_Norm_data_CLL2_HC1_all_cl_20_k1/...", "CytofRUV_Norm_data_CLL2_HC1_all_cl_20_k5/..." rather than
        # "CytofRUV_Norm_data_CLL2_HC1_all_cl_20_k"
        # wd=paste(data_dir[i],k[j],"/",sep="")
        # EMD_metric=readRDS(paste(wd,"EMD_metric_comp_norm_all_samples_all_cells.rds",sep=""))
        EMD_metric=readRDS(data_dir[[i]][[j]])
        EMD_tmp=data.frame(patient=row.names(EMD_metric),EMD_metric)
        colnames(EMD_tmp) <- gsub("^X", "",  colnames(EMD_tmp))
        cytof_methods <- append(cytof_methods, paste0("cytofruv_", "k", k[j]))
        EMD[[paste0("cytofruv_", "k", k[j])]]=EMD_tmp
      }

    } else {
      ## Reading additional data files for EMD, e.g. cytonorm, batchadjust
      EMD_data_file=readRDS(data_dir[[i]])
      EMD_data=data.frame(patient=row.names(EMD_data_file),EMD_data_file)
      colnames(EMD_data) <- gsub("^X", "",  colnames(EMD_data))
      EMD[[method_names[i]]]=EMD_data
    }
  }

  ## Merging all the data
  EMD_all=melt(EMD,id.var = c("patient"))
  EMD_all$method=EMD_all$L1
  ggdf <- melt(EMD, id.var = c("patient"),
               value.name = "EMD", variable.name = "antigen")

  method_no_cytofruv=method_names[method_names%!in%c("cytofruv")]

  ggdf$method=factor(ggdf$L1,levels=c("raw", method_no_cytofruv, cytof_methods))

  ## Selecting CLL patients only
  ggdf2=ggdf[ggdf$patient%in%filter_patients,]

  ggdf2s=ggdf2[ggdf2$method%in%c("raw", method_no_cytofruv, cytof_methods),]
  ggdf2s$sample=""

  for (i in (1:length(filter_patients))){
    ggdf2s$sample[ggdf2s$patient==filter_patients[i]]=paste0(x_axis_name,i)
  }

  ## Draw Plot
  ggplot(ggdf2s,aes(y=EMD,x=sample,color=method))+geom_boxplot(aes(colour = method)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    theme(text = element_text(size=40),
          axis.text.x = element_text(angle=90, hjust=1)) + scale_color_manual(values = cols) +
    theme(aspect.ratio=1)
}

# Get a boolean vector - true if the element of x IS NOT y.
'%!in%' <- function(x,y) return(!('%in%'(x,y)))
