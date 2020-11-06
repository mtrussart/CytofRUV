#' @rdname plotSilhouette
#' @title plotSilhouette
#'
#' Insert description of plotSilhouette here
#'
#' @param daf a \code{\link[SingleCellExperiment]{SingleCellExperiment}}, this was generated after running examples present
#'   in the vignette.
#' @param raw_dir String: Directory to file containing raw data.
#' @param silh_fink_dir String: Directory to file containing silhouette fink data.
#' @param data_dir List: Strings of directories that contain appropriate data.
#'   IMPORTANT: Order of directory strings should correlate with order of method names in method_names.
#'   Note the directory for cytofRUV processed data is a vector of strings.
#'   e.g. c("../Cytofnorm_Cytof_Package/EMD_metric_comp_norm_all_samples_all_cells.rds",
#'          "../BatchAdjust/BatchAdjust_95p_and_DG23_new_trial/EMD_metric_comp_norm_all_samples_all_cells.rds"
#'          c("CytofRUV_Norm_data_CLL2_HC1_all_cl_20_k5/EMD_metric_comp_norm_all_samples_all_cells.rds",
#'            "CytofRUV_Norm_data_CLL2_HC1_all_cl_20_k15/EMD_metric_comp_norm_all_samples_all_cells.rds",
#'            "CytofRUV_Norm_data_CLL2_HC1_all_cl_20_k20/EMD_metric_comp_norm_all_samples_all_cells.rds"
#'          )
#'          )
#' @param method_names Vector: Strings that describe method names.
#'   e.g. c(cytonorm, batchadj, cytofruv)
#' @param k Vector: cluster numbers to be inspected, default val: c(5,10,15)
#' @param filter_patients Vector: relevant patient IDs  e.g. val = c("LL1_B1","LL2_B1","LL3_B1") that is,
#'   display patient IDs that contain: "LL1_B1","LL2_B1","LL3_B1".
#' @param x_axis_name String: Will be shared by all x-axis labels will share in common. E.g. "CLL", results in
#'   x-axis items = "CLL1, CLL2, CLL3"
#' @param cols Vector: Strings of colours for plot annotations. Default val:
#'   c("raw"="black","batchadj"="indianred4","cytonorm"="darkmagenta","cytofruv_k1"="#084081","cytofruv_k3"="#0868AC","cytofruv_k5"="#2B8CBE","cytofruv_k7"="#4EB3D3","cytofruv_k10"="#7BCCC4","cytofruv_k12"= "#A8DDB5","cytofruv_k15"="darkgreen","cytofruv_k20"="seagreen","cytofruv_k25"="yellowgreen")
#'
#' @export

# filter_patients should replace x_axis_name.

# Test on: All_analysis.Rdata
plotSilhouette <- function(raw_dir, silh_fink_dir, data_dir, method_names, k, filter_patients, x_axis_name,
                           cols=c("raw"="black","batchadj"="indianred4","cytonorm"="darkmagenta","cytofruv_k1"="#084081","cytofruv_k3"="#0868AC","cytofruv_k5"="#2B8CBE","cytofruv_k7"="#4EB3D3","cytofruv_k10"="#7BCCC4","cytofruv_k12"= "#A8DDB5","cytofruv_k15"="darkgreen","cytofruv_k20"="seagreen","cytofruv_k25"="yellowgreen")) {
  ## Reading the raw file Silhouette for all samples
  silh_finck<-readRDS(silh_fink_dir)
  CLL<-list()

  for (p in 1:length(filter_patients)){
    CLL[[filter_patients[p]]]<-data.frame()
    CLL[[filter_patients[p]]][1,'method']<-"raw"
    CLL[[filter_patients[p]]][1,'bio']<-as.numeric(silh_finck[p,2])
    CLL[[filter_patients[p]]][1,'batch']<-as.numeric(silh_finck[p,3])
  }

  cytof_methods <- vector()
  ## Reading the Silhouette for all samples
  for (p in (1:length(filter_patients))){
    for (i in (1:length(method_names))) {

      row_num=dim(CLL[[filter_patients[p]]])[1]+1

      if (method_names[i]=="cytofruv") {
        ## Reading the CytofRUV files EMD
        for (j in (1:length(k))){
          if (length(cytof_methods)!=length(k)) {cytof_methods <- append(cytof_methods, paste0("cytofruv_", "k", k[j]))}

          silh<-readRDS(data_dir[[i]][[j]])

          CLL[[filter_patients[p]]][row_num,'method']<-paste("cytofruv_k",k[j],sep="")
          CLL[[filter_patients[p]]][row_num,'bio']<-as.numeric(silh[p,2])
          CLL[[filter_patients[p]]][row_num,'batch']<-as.numeric(silh[p,3])
          row_num=row_num+1
        }
      } else {
        ## Reading the files for additional methods
        temp_data<-readRDS(data_dir[[i]])
        CLL[[filter_patients[p]]][(row_num),'method']<-method_names[i]
        CLL[[filter_patients[p]]][(row_num),'bio']<-as.numeric(temp_data[p,2])
        CLL[[filter_patients[p]]][(row_num),'batch']<-as.numeric(temp_data[p,3])
      }
    }

    # Formating
    method_no_cytofruv<-method_names[method_names%!in%c("cytofruv")]

    CLL[[filter_patients[p]]]$method<-factor(CLL[[filter_patients[p]]]$method,levels=c("raw", method_no_cytofruv, cytof_methods))
    CLL[[filter_patients[p]]]$bio<-as.numeric(CLL[[filter_patients[p]]]$bio)
    CLL[[filter_patients[p]]]$batch<-as.numeric(CLL[[filter_patients[p]]]$batch)
  }

  method_no_cytofruv<-method_names[method_names%!in%c("cytofruv")]

  ## Merging all the data
  All_patients_CLL<-melt(CLL,id.var=c("bio","batch","method"))

  ## Plotting CLL samples
  All_patients_CLL_select<-All_patients_CLL[All_patients_CLL$method%in%c("raw", method_no_cytofruv, cytof_methods),]

  for (i in (1:length(filter_patients))){
    All_patients_CLL_select$L1[All_patients_CLL_select$L1==filter_patients[i]]=paste0(x_axis_name,i)
  }

  gg=ggplot(data = All_patients_CLL_select, aes(x = bio, y = batch))+ geom_point(aes(colour = method)) + facet_wrap(~L1) + scale_color_manual(values = cols)#+ scale_color_manual(values=c("black","indianred4","darkmagenta", "#084081","#0868AC", "#2B8CBE", "#4EB3D3", "#7BCCC4", "#A8DDB5" ,"darkgreen","seagreen"))
  gg+ theme(aspect.ratio=1)+theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(text = element_text(size=20))
}
