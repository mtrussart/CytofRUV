#' @rdname plotEMD
#' @title plotEMD
#'
#' Insert description of plotEMD here
#'
#' @param raw_dir String - Directory to file containing raw data.
#' @param cytonorm_dir String - Directory to file containing cytonorm data
#' @param batchadj_dir String - Directory to file containing cytonorm data
#' @param cytofRUV_dir String - Directory to FOLDER that appropriate files containing cytofRUV processed data.
#' @param k Vector of cluster numbers to be inspected, default val = c(5,10,15)
#' @param filter_patients Vector of relevant patient IDs  e.g. val = c("LL1_B1","LL2_B1","LL3_B1") that is, display patient IDs that contain: "LL1_B1","LL2_B1","LL3_B1"
#' @param x_axis_name String which all x labels will share in common. E.g. "CLL", resulting x-axis items = "CLL1, CLL2, CLL3"
#' @param cols Vector of colours for plot annotations. Default val= c("raw"="black","batchadj"="indianred4","cytonorm"="darkmagenta","cytofruv_k1"="#084081","cytofruv_k3"="#0868AC","cytofruv_k5"="#2B8CBE","cytofruv_k7"="#4EB3D3","cytofruv_k10"="#7BCCC4","cytofruv_k12"= "#A8DDB5","cytofruv_k15"="darkgreen","cytofruv_k20"="seagreen","cytofruv_k25"="yellowgreen")
#'
#' @export
#'
############################### Plotting the metrics ###########################
library(ggplot2)
library(reshape2)

'%!in%' <- function(x,y)!('%in%'(x,y))

##################### EMD plot ##################

### Dataset from CLL+HC samples with 20 clusters


# Example Inputs
# TODO: Change outdir -> How Necessary is outdir? What does it do?
outdir="/stornext/Home/data/allstaff/t/trussart.m/R_scripts/RUV_design_experiment/RUV1b_RUV3b/Metrics_of_performances_all_methods/All_samples_DG23_VDBR996_anchor_CytofPackage_new_trial/"
my_dir="/stornext/Bioinf/data/lab_speed/Marie/RUV_experiment_design/Norm_Raw_RUV1b_RUV3b_simult/CytofRUV/"
raw_dir=paste0(my_dir,"CytofRUV_Norm_data_CLL2_HC1_all_cl_20_k5/EMD_metric_comp_raw_all_samples_all_cells.rds")
cytonorm_dir="/stornext/Home/data/allstaff/t/trussart.m/R_scripts/RUV_design_experiment/RUV1b_RUV3b/Cytofnorm_Cytof_Package/EMD_metric_comp_norm_all_samples_all_cells.rds"
batchadj_dir="/stornext/Bioinf/data/lab_speed/Marie/RUV_experiment_design/Norm_Raw_RUV1b_RUV3b_simult/BatchAdjust/BatchAdjust_95p_and_DG23_new_trial/EMD_metric_comp_norm_all_samples_all_cells.rds"
cytofRUV_dir=paste0(my_dir,"CytofRUV_Norm_data_CLL2_HC1_all_cl_20_k");
  # cytofRUV_dir is the folder that contains all files with appropriate K values.

filter_patients = c("LL1_B1","LL2_B1","LL3_B1")
x_axis_name="CLL"
cols <- c("raw"="black","batchadj"="indianred4","cytonorm"="darkmagenta","cytofruv_k1"="#084081","cytofruv_k3"="#0868AC","cytofruv_k5"="#2B8CBE","cytofruv_k7"="#4EB3D3","cytofruv_k10"="#7BCCC4","cytofruv_k12"= "#A8DDB5","cytofruv_k15"="darkgreen","cytofruv_k20"="seagreen","cytofruv_k25"="yellowgreen")

# plotEMD(raw_dir,cytonorm_dir,batchadj_dir,cytofRUV_dir,filter_patients,x_axis_name)

method_names=c("cytonorm", "batchadjust", "cytofruv")

## Note that cytofRUV_dir is a directory to a folder and not a singular file!
data_dir=c(cytonorm_dir, batchadj_dir, cytofRUV_dir)

plotEMD <- function(raw_dir, method_names=c("cytofRUV"), data_dir, k=c(5,10,15), filter_patients, x_axis_name,
                    cols=c("raw"="black","batchadj"="indianred4","cytonorm"="darkmagenta","cytofruv_k1"="#084081","cytofruv_k3"="#0868AC","cytofruv_k5"="#2B8CBE","cytofruv_k7"="#4EB3D3","cytofruv_k10"="#7BCCC4","cytofruv_k12"= "#A8DDB5","cytofruv_k15"="darkgreen","cytofruv_k20"="seagreen","cytofruv_k25"="yellowgreen")
) {
  EMD=list()

  # TODO: ASSERT FOR INPUTS!
    # i.e. CYTOFRUV is required as a method

  ### Loading EMD for each method
  # TODO: method is hardcoded
  # method=c("raw","cytofruv_k5","cytofruv_k10","cytofruv_k15","cytofruv_k20","cytofruv_k25","cytonorm","batchadjust")

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
        wd=paste(data_dir[i],k[i],"/",sep="")
        EMD_metric=readRDS(paste(wd,"EMD_metric_comp_norm_all_samples_all_cells.rds",sep=""))
        EMD_tmp=data.frame(patient=row.names(EMD_metric),EMD_metric)
        colnames(EMD_tmp) <- gsub("^X", "",  colnames(EMD_tmp))
        append(cytof_methods, paste0("cytofruv_", i+1))
        EMD[[paste0("cytofruv_", j+1)]]=EMD_tmp
      }

    } else {
      ## Reading additional data files for EMD, e.g. cytonorm, batchadjust
      EMD_data_file=readRDS(data_dir[i])
      EMD_data=data.frame(patient=row.names(EMD_data_file),EMD_data_file)
      colnames(EMD_data) <- gsub("^X", "",  colnames(EMD_data))
      EMD[[method_names[i]]]=EMD_data
    }
  }

  # ## Reading the cytonorm file EMD
  # EMD_cytofnorm_file=readRDS(cytonorm_dir)
  # EMD_cytofnorm=data.frame(patient=row.names(EMD_cytofnorm_file),EMD_cytofnorm_file)
  # colnames(EMD_cytofnorm) <- gsub("^X", "",  colnames(EMD_cytofnorm))
  # EMD[['cytonorm']]=EMD_cytofnorm
  #
  # ## Reading the batchadjust file EMD
  # EMD_batchadj_file=readRDS(batchadj_dir)
  # EMD_batchadj=data.frame(patient=row.names(EMD_batchadj_file),EMD_batchadj_file)
  # colnames(EMD_batchadj) <- gsub("^X", "",  colnames(EMD_batchadj))
  # EMD[["batchadj"]]=EMD_batchadj

  ## Merging all the data
  EMD_all=melt(EMD,id.var = c("patient"))
  EMD_all$method=EMD_all$L1
  ggdf <- melt(EMD, id.var = c("patient"),
               value.name = "EMD", variable.name = "antigen")

  method_no_cytofruv=method_names[method_names%!in%c("cytofruv")]

  # TODO: Do we need k1-k20 if user selects specific Ks only?
  #c(c("raw","batchadj","cytonorm"), cytof_methods, "cytofruv_k1","cytofruv_k3","cytofruv_k5","cytofruv_k7","cytofruv_k10","cytofruv_k12","cytofruv_k15","cytofruv_k20","cytofruv_k25")
  ggdf$method=factor(ggdf$L1,levels=c("raw", method_no_cytofruv, cytof_methods))

  ## Selecting CLL patients only
  ggdf2=ggdf[ggdf$patient%in%filter_patients,]

  ## Plotting EMD comparing raw, batchadj,cytonorm and 3 cytofruv k=5,10,15
  #  ggdf2s=ggdf2[ggdf2$method%in%c("raw","batchadj","cytonorm","cytofruv_k5","cytofruv_k10","cytofruv_k15"),]
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
