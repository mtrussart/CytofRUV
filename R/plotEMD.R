############################### Plotting the metrics ###########################
library(ggplot2)
library(reshape2)
cols <- c("raw"="black","batchadj"="indianred4","cytonorm"="darkmagenta","cytofruv_k1"="#084081","cytofruv_k3"="#0868AC","cytofruv_k5"="#2B8CBE","cytofruv_k7"="#4EB3D3","cytofruv_k10"="#7BCCC4","cytofruv_k12"= "#A8DDB5","cytofruv_k15"="darkgreen","cytofruv_k20"="seagreen","cytofruv_k25"="yellowgreen")


##################### EMD plot ##################

### Dataset from CLL+HC samples with 20 clusters

outdir="/stornext/Home/data/allstaff/t/trussart.m/R_scripts/RUV_design_experiment/RUV1b_RUV3b/Metrics_of_performances_all_methods/All_samples_DG23_VDBR996_anchor_CytofPackage_new_trial/"
my_dir="/stornext/Bioinf/data/lab_speed/Marie/RUV_experiment_design/Norm_Raw_RUV1b_RUV3b_simult/CytofRUV/"
raw_dir=paste(dir,"CytofRUV_Norm_data_CLL2_HC1_all_cl_20_k5/EMD_metric_comp_raw_all_samples_all_cells.rds",sep="")
raw_dir=paste(dir,"CytofRUV_Norm_data_CLL2_HC1_all_cl_20_k5/EMD_metric_comp_raw_all_samples_all_cells.rds",sep="")
cytonorm_dir="/stornext/Home/data/allstaff/t/trussart.m/R_scripts/RUV_design_experiment/RUV1b_RUV3b/Cytofnorm_Cytof_Package/EMD_metric_comp_norm_all_samples_all_cells.rds"
batchadj_dir

# TODO: Change outdir -> How Necessary is outdir? What does it do?
# TODO: params - dir-> raw_data_dir, k->k,

plotEMD <- function(raw_dir, cytonorm_dir, batchadj_dir, raw_data_dir, k=c(5,10,15,20,25)) {

  ### Loading EMD for each method
  # Vector of k values to use
  EMD=list()
  # TODO: method is hardcoded
  method=c("raw","cytofruv_k5","cytofruv_k10","cytofruv_k15","cytofruv_k20","cytofruv_k25","cytonorm","batchadjust")
  method=c("raw","cytofruv_k5","cytofruv_k10","cytofruv_k15","cytofruv_k20","cytofruv_k25","cytonorm","batchadjust")
  ## Reading the raw file EMD
  EMD_raw_file=readRDS()
  EMD_raw=data.frame(patient=row.names(EMD_raw_file),EMD_raw_file)
  colnames(EMD_raw) <- gsub("^X", "",  colnames(EMD_raw))
  EMD[["raw"]]=EMD_raw
  ## Reading the cytonorm file EMD
  EMD_cytofnorm_file=readRDS(cytonorm_dir)
  EMD_cytofnorm=data.frame(patient=row.names(EMD_cytofnorm_file),EMD_cytofnorm_file)
  colnames(EMD_cytofnorm) <- gsub("^X", "",  colnames(EMD_cytofnorm))
  EMD[['cytonorm']]=EMD_cytofnorm
  ## Reading the batchadjust file EMD
  EMD_batchadj_file=readRDS("/stornext/Bioinf/data/lab_speed/Marie/RUV_experiment_design/Norm_Raw_RUV1b_RUV3b_simult/BatchAdjust/BatchAdjust_95p_and_DG23_new_trial/EMD_metric_comp_norm_all_samples_all_cells.rds")
  EMD_batchadj=data.frame(patient=row.names(EMD_batchadj_file),EMD_batchadj_file)
  colnames(EMD_batchadj) <- gsub("^X", "",  colnames(EMD_batchadj))
  EMD[["batchadj"]]=EMD_batchadj
  ## Reading the CytofRUV files EMD
  for (i in (1:length(k)) ){
    wd=paste(dir,"CytofRUV_Norm_data_CLL2_HC1_all_cl_20_k",k[i],"/",sep="")
    EMD_metric=readRDS(paste(wd,"EMD_metric_comp_norm_all_samples_all_cells.rds",sep=""))
    EMD_tmp=data.frame(patient=row.names(EMD_metric),EMD_metric)
    colnames(EMD_tmp) <- gsub("^X", "",  colnames(EMD_tmp))
    EMD[[method[i+1]]]=EMD_tmp
  }
  ## Merging all the data
  EMD_all=melt(EMD,id.var = c("patient"))
  EMD_all$method=EMD_all$L1
  ggdf <- melt(EMD, id.var = c("patient"),
               value.name = "EMD", variable.name = "antigen")
  ggdf$method=factor(ggdf$L1,levels=c("raw","batchadj","cytonorm","cytofruv_k1","cytofruv_k3","cytofruv_k5","cytofruv_k7","cytofruv_k10","cytofruv_k12","cytofruv_k15","cytofruv_k20","cytofruv_k25"))
  ## Selecting CLL patients only
  ggdf2=ggdf[ggdf$patient%in%c("LL1_B1","LL2_B1","LL3_B1"),]
  ## Plotting EMD comparing raw, batchadj,cytonorm and 3 cytofruv k=5,10,15
  ggdf2s=ggdf2[ggdf2$method%in%c("raw","batchadj","cytonorm","cytofruv_k5","cytofruv_k10","cytofruv_k15"),]
  ggdf2s$sample=""
  ggdf2s$sample[ggdf2s$patient=="LL1_B1"]="CLL1"
  ggdf2s$sample[ggdf2s$patient=="LL2_B1"]="CLL2"
  ggdf2s$sample[ggdf2s$patient=="LL3_B1"]="CLL3"
  ggplot(ggdf2s,aes(y=EMD,x=sample,color=method))+geom_boxplot(aes(colour = method)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    theme(text = element_text(size=40),
          axis.text.x = element_text(angle=90, hjust=1)) + scale_color_manual(values = cols) +
    theme(aspect.ratio=1)
}
