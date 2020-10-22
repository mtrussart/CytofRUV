#' @rdname plotHellingerDistance
#' @title plotHellingerDistance
#'
#' Insert description of plotHellingerDistance here
#'
#' @param daf a \code{\link[SingleCellExperiment]{SingleCellExperiment}}, this was generated after running examples present
#'   in the vignette.
#' @param raw_dist_dir String: Directory to file containing raw distance data.
#' @param method_names Vector: Strings that describe method names.
#' @param data_dir List: Strings of directories that contain appropriate data.
#'   IMPORTANT: Order of directory strings should correlate with order of method names in method_names.
#'   Note the directory for cytofRUV processed data is a vector of strings.
#'   e.g. c("../Cytofnorm_Cytof_Package/Dist_props_table_norm_all_samples.rds",
#'          "../BatchAdjust/BatchAdjust_95p_and_DG23_new_trial/Dist_props_table_norm_all_samples.rds"
#'          c("CytofRUV_Norm_data_CLL2_HC1_all_cl_20_k5/Dist_props_table_norm_all_samples.rds",
#'            "CytofRUV_Norm_data_CLL2_HC1_all_cl_20_k15/Dist_props_table_norm_all_samples.rds",
#'            "CytofRUV_Norm_data_CLL2_HC1_all_cl_20_k20/Dist_props_table_norm_all_samples.rds"
#'          )
#'          )
#' @param k Vector: cluster numbers to be inspected, default val: c(5,10,15)
#' @param filter_patients Vector: relevant patient IDs  e.g. val = c("LL1_B1","LL2_B1","LL3_B1") that is,
#'   display patient IDs that contain: "LL1_B1","LL2_B1","LL3_B1".
#' @param x_axis_name String: Will be shared by all x-axis labels will share in common. E.g. "CLL", results in
#'   x-axis items = "CLL1, CLL2, CLL3"
#' @param cols Vector: Strings of colours for plot annotations. Default val:
#'   c("raw"="black","batchadj"="indianred4","cytonorm"="darkmagenta","cytofruv_k1"="#084081","cytofruv_k3"="#0868AC","cytofruv_k5"="#2B8CBE","cytofruv_k7"="#4EB3D3","cytofruv_k10"="#7BCCC4","cytofruv_k12"= "#A8DDB5","cytofruv_k15"="darkgreen","cytofruv_k20"="seagreen","cytofruv_k25"="yellowgreen")
#'
#' @export

plotHellingerDistance <- function(daf, raw_dist_dir, method_names, data_dir, k=NULL, filter_patients, x_axis_name,
                                  cols=c("raw"="black","batchadj"="indianred4","cytonorm"="darkmagenta","cytofruv_k1"="#084081","cytofruv_k3"="#0868AC","cytofruv_k5"="#2B8CBE","cytofruv_k7"="#4EB3D3","cytofruv_k10"="#7BCCC4","cytofruv_k12"= "#A8DDB5","cytofruv_k15"="darkgreen","cytofruv_k20"="seagreen","cytofruv_k25"="yellowgreen")) {

  ## ======= Function Parameter Assertions! =======
  if(!file.exists(raw_dist_dir)) {stop("The directory that was passed and meant to contain raw distance data does not exist!")}

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

  # Ensure that k is set when cytofruv is stated to be used.
  if ("cytofruv"%in%method_names && is.null(k)) {stop(paste("The method cytofruv was attempted to be used but k values were not provided!"))}

  # Ensure patient id is in daf
  # for (i in (1:length(filter_patients))) {
  #   if (!filter_patients[[i]]%in%daf$patient_id) {stop(paste("The patient to filter by does not exist in the current dataset."))}
  # }

  ## ======= Main Function ======= ##

  Dist=list()

  ## Reading the raw file Hellinger
  Dist_raw_file=readRDS(raw_dist_dir)
  Dist_raw=data.frame(patient=colnames(Dist_raw_file),Hell=Dist_raw_file[1,])
  Dist[["raw"]]=Dist_raw

  cytof_methods <- vector()
  for (i in (1:length(method_names)) ){
    if (method_names[i]=="cytofruv") {
      ## Reading the cytofruv files Hellinger
      for (j in (1:length(k)) ){
        # wd=paste(dir,"CytofRUV_Norm_data_CLL2_HC1_all_cl_20_k",k[i],"/",sep="")
        # Dist_metric=readRDS(paste(wd,"Dist_props_table_norm_all_samples.rds",sep=""))
        Dist_metric=readRDS(data_dir[[i]][[j]])
        Dist_tmp=data.frame(patient=colnames(Dist_metric),Hell=Dist_metric[1,])
        cytof_methods <- append(cytof_methods, paste0("cytofruv_", "k", k[j]))
        Dist[[paste0("cytofruv_", "k", k[j])]]=Dist_tmp
        ## pass in prefix for directory online
      }

    } else {
      ## Reading additional Hellinger containing files e.g. cytonorm, batchadjust
      Dist_data_file=readRDS(data_dir[[i]])
      Dist_data=data.frame(patient=colnames(Dist_data_file),Hell=Dist_data_file[1,])
      Dist[[method_names[i]]]=Dist_data
    }
  }

  method_no_cytofruv=method_names[method_names%!in%c("cytofruv")]

  ## Merging all the data
  Dist_all=melt(Dist,id.var = c("patient"))
  Dist_all$method=factor(Dist_all$L1,levels=c("raw", method_no_cytofruv, cytof_methods))

  ## Plotting CLL samples
  ggdf2=Dist_all[Dist_all$patient%in%filter_patients,]
  ggdf3=ggdf2[ggdf2$method %in% c("raw", method_no_cytofruv, cytof_methods),]
  ggdf3$sample=""

  for (i in (1:length(filter_patients))){
    ggdf3$sample[ggdf3$patient==filter_patients[i]]=paste0(x_axis_name,i)
  }

  ggdf3$distance=ggdf3$value

  # method refers to the method column in ggdf3!
  gg=ggplot(ggdf3,aes(y=distance, x=method, color=method))+ geom_point(aes(colour = method),size=3) + scale_color_manual(values = cols)
  gg+  theme_bw() + theme(aspect.ratio=1) +
    theme(text = element_text(size=40),
          axis.text.x = element_text(angle=90, hjust=1)) + theme(panel.background = element_blank()) +
    facet_wrap(~sample) +
    theme(axis.text.x = element_blank())
}
