############################ METRICS FOR EVALUATION OF METHODS ############################
run_metrics <- function(raw_files,cluster_to_select_raw_data,norm_files,outdir,cluster_to_select_norm_data,rep_samples,nb_cells=10000,props_table_raw,props_table_norm){
  ### METRICS TO ANALYSE PERFORMANCE OF THE NORMALISATION
  library(cluster)
  library(MASS)
  library(ruv)
  library(rsvd)
  library("gridExtra")
  library("ggpubr")
  library("tidyverse")
  library("RColorBrewer")
  library(purrr)
  library(matrixStats)
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)
  library(pheatmap)

  files <- c(raw_files, norm_files)
  file=c(1,2)
  data_list <- map(files, ~readRDS(.))
  raw_data=data_list[[1]]
  norm_data=data_list[[2]]
  colnames(norm_data) <- gsub("^X", "",  colnames(norm_data))
  colnames(raw_data) <- gsub("^X", "",  colnames(raw_data))
  nb_markers=dim(raw_data)[2]-2
  nclust=dim(props_table_norm)[1]
  nsample=dim(props_table_norm)[2]
  all_samples=colnames(props_table_norm)
  all_cluster=rownames(props_table_norm)

  ############################ 1.Silhouette and LDA on the Silhouette
  silhouette_specific_clusters <- function (cluster_to_select,rep_samples,nb_cells,data)
  {
    sample_ids=data$sample
    cell_clustering1=data$cluster
    expr= as.matrix(data[, 3:ncol(data)])
    colnames(expr) <- gsub("^X", "",  colnames(expr))
    rep_sample_ids=sample_ids[sample_ids%in%rep_samples]
    rep_expr=expr[sample_ids%in%rep_samples,]
    rep_cell_clustering=cell_clustering1[sample_ids%in%rep_samples]
    rep_sample_ids=rep_sample_ids[rep_cell_clustering%in%cluster_to_select]
    rep_expr=rep_expr[rep_cell_clustering%in%cluster_to_select,]
    rep_cell_clustering=rep_cell_clustering[rep_cell_clustering%in%cluster_to_select]
    ## Overall effect of bio vs batch with subsampling
    rep_data=data.frame(rep_expr,clust=as.numeric(rep_cell_clustering),batch=as.numeric(substr(rep_sample_ids,nchar(as.character(rep_sample_ids)), nchar(as.character(rep_sample_ids)))))
    colnames(rep_data) <- gsub("^X", "",  colnames(rep_data))
    base::set.seed(1234)
    subs_rep_data=rep_data[sample(nrow(rep_data), nb_cells), ]
    ## Silhouette
    over_ss_biology=silhouette(subs_rep_data$clust,dist(subs_rep_data[,1:ncol(data)]))
    over_ss_batch=silhouette(subs_rep_data$batch,dist(subs_rep_data[,1:ncol(data)]))
    over_ss<- c(bio = mean(over_ss_biology[,"sil_width"]), batch = mean(over_ss_batch[,"sil_width"]))
    #saveRDS(over_ss, paste(outdir,"Sil_Rep_",paste(c(rep_samples),collapse="_"),"_Clust_",paste(c(cluster_to_select),collapse="_"),"_",datatype,"_data_over_ss.rds",sep=""))
    return (over_ss)
  }

  silhouette_and_LDA_specific_clusters <- function (cluster_to_select,rep_samples,nb_cells,data,datatype){
    sample_ids=data$sample
    cell_clustering1=data$cluster
    expr= as.matrix(data[, 3:ncol(data)])
    colnames(expr) <- gsub("^X", "",  colnames(expr))
    rep_sample_ids=sample_ids[sample_ids%in%rep_samples]
    rep_expr=expr[sample_ids%in%rep_samples,]
    rep_cell_clustering=cell_clustering1[sample_ids%in%rep_samples]
    rep_sample_ids=rep_sample_ids[rep_cell_clustering%in%cluster_to_select]
    rep_expr=rep_expr[rep_cell_clustering%in%cluster_to_select,]
    rep_cell_clustering=rep_cell_clustering[rep_cell_clustering%in%cluster_to_select]
    ## Overall effect of bio vs batch with subsampling
    rep_data=data.frame(rep_expr,clust=as.numeric(rep_cell_clustering),batch=as.numeric(substr(rep_sample_ids,nchar(as.character(rep_sample_ids)), nchar(as.character(rep_sample_ids)))))
    colnames(rep_data) <- gsub("^X", "",  colnames(rep_data))
    base::set.seed(1234)
    subs_rep_data=rep_data[sample(nrow(rep_data), nb_cells), ]
    ## Silhouette
    over_ss_biology=silhouette(subs_rep_data$clust,dist(subs_rep_data[,1:ncol(data)]))
    over_ss_batch=silhouette(subs_rep_data$batch,dist(subs_rep_data[,1:ncol(data)]))
    over_ss<- c(bio = mean(over_ss_biology[,"sil_width"]), batch = mean(over_ss_batch[,"sil_width"]))
    saveRDS(over_ss, paste(outdir,"Sil_Rep_",paste(c(rep_samples),collapse="_"),"_Clust_",paste(c(cluster_to_select),collapse="_"),"_",datatype,"_data_over_ss.rds",sep=""))
    ### LDS plot
    batches=c(1,2)
    batch_data=subs_rep_data[subs_rep_data$batch%in%batches,]
    subs.lda2batch <- lda(clust + batch ~., data =batch_data )
    subs.lda.values <- predict(subs.lda2batch, batch_data[,1:ncol(data)])
    #convert to data frame
    color_batch=c("#0072B2","#D55E00")
    newbatch_data <- data.frame(batch = factor(batch_data$batch),cluster = factor(batch_data$clust), LDA1 = subs.lda.values$x[,1],LDA2 = subs.lda.values$x[,2],samp=paste("Rep_",batch_data$batch,"_Clust_",batch_data$clust,sep=""))
    ggplot(newbatch_data) + geom_point(aes(LDA1,LDA2,colour = batch,shape=cluster), size = 2.5) +
      scale_color_manual(values = color_batch) +
      guides(col = guide_legend(overide.aes = list(alpha = 1, size = 5))) +
      theme_bw() +
      theme(aspect.ratio = 1,
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "lightgrey",size = 0.5))+
      theme(text = element_text(size=20), legend.title = element_text(size = 20),
            legend.text = element_text(size = 20),
            axis.text.x = element_text(angle=90, hjust=1))
    ggsave(paste(outdir,"LDA_Rep_",paste(c(rep_samples),collapse="_"),"_Clust_",paste(c(cluster_to_select),collapse="_"),"_",datatype,"_data.png",sep=""))
  }

  # ## All samples
  ind=1
  sil_raw_all_samples=matrix(nrow=nsample/2,ncol=3)
  colnames(sil_raw_all_samples)=c("samples","bio","batch")
  sil_raw_all_samples[,1]=substring(colnames(props_table_norm)[seq(1,nsample,2)],first=1,last = 4)
  sil_norm_all_samples=sil_raw_all_samples
  for (s in seq(1,nsample,2)){
    rep_samples_sel=colnames(props_table_norm)[c(s,s+1)]
    sraw=silhouette_specific_clusters(all_cluster,rep_samples_sel,nb_cells,data=data_list[[1]])
    sil_raw_all_samples[ind,2:3]=sraw
    snorm=silhouette_specific_clusters(all_cluster,rep_samples_sel,nb_cells,data=data_list[[2]])
    sil_norm_all_samples[ind,2:3]=snorm
    ind=ind+1
  }
  saveRDS(sil_raw_all_samples, paste(outdir,"Sil_all_samples_all_clust_raw_data_over_ss.rds",sep=""))
  saveRDS(sil_norm_all_samples, paste(outdir,"Sil_all_samples_all_clust_norm_data_over_ss.rds",sep=""))

  ############################ 2. EMD
  EMD_metric <- function (data,binSize=0.005){
    nb_markers=dim(data[[1]])[2]
    distr <- list()
    for(file in c(1:length(files))){
      distr[[file]] <-  apply(data[[file]][,3:nb_markers],2,function(x){
        graphics::hist(x,
                       breaks = seq(-500,500,by=binSize),
                       plot = FALSE)$counts
      })
    }

    distances <- matrix(nrow=1,ncol=(nb_markers-3+1))
    colnames(distances)=colnames(data[[1]])[3:nb_markers]
    for(marker in c(1:(nb_markers-3+1))){
      distances[marker] <-
        emdist::emd2d(
          matrix(distr[[1]][,marker]),
          matrix(distr[[2]][,marker]))
    }
    return(distances)
  }

  ## EMD metric for rep samples per cluster
  EMD_per_cluster <- function(rep1,rep2){
    rep_comp_clust=matrix(nrow=nclust,ncol=nb_markers)
    colnames(rep_comp_clust)= colnames(data_list[[1]])[3:ncol(data_list[[1]])]
    row.names(rep_comp_clust)=paste("Cl_",seq(1,nclust),sep="")
    for (c in c(1:nclust)){
      rep1_clust_c=rep1[rep1$clust%in%c,]
      rep2_clust_c=rep2[rep2$clust%in%c,]
      rep_clust_c=list(rep1_clust_c,rep2_clust_c)
      rep_comp_clust[c,]=EMD_metric(rep_clust_c)
    }
    return(rep_comp_clust)
  }

  ## EMD metric for rep samples per cluster
  cluster_to_use=norm_data$cluster # !!!! IMPORTANT
  raw_data_tmp=raw_data
  raw_data_tmp$cluster=cluster_to_use
  raw_rep1=raw_data_tmp[raw_data$sample%in%rep_samples[1],]
  raw_rep2=raw_data_tmp[raw_data$sample%in%rep_samples[2],]
  norm_rep1=norm_data[norm_data$sample%in%rep_samples[1],]
  norm_rep2=norm_data[norm_data$sample%in%rep_samples[2],]
  raw_rep_comp_clust=EMD_per_cluster(raw_rep1,raw_rep2)
  norm_rep_comp_clust=EMD_per_cluster(norm_rep1,norm_rep2)
  saveRDS(raw_rep_comp_clust, paste(outdir,"EMD_metric_raw_data_per_cluster_",paste(c(rep_samples),collapse="_"),".rds",sep=""))
  saveRDS(norm_rep_comp_clust, paste(outdir,"EMD_metric_norm_data_per_cluster_",paste(c(rep_samples),collapse="_"),".rds",sep=""))
  ## EMD metric for all rep samples all cells
  ind=1
  EMD_metric_comp_raw_all_samples=matrix(nrow=nsample/2,ncol=nb_markers)
  colnames(EMD_metric_comp_raw_all_samples)=colnames(data_list[[1]])[3:ncol(data_list[[1]])]
  row.names(EMD_metric_comp_raw_all_samples)=substring(colnames(props_table_norm)[seq(1,nsample,2)],2)
  EMD_metric_comp_norm_all_samples=EMD_metric_comp_raw_all_samples
  for (s in seq(1,nsample,2)){
    rep_samples_sel=colnames(props_table_norm)[c(s,s+1)]
    raw_rep1=data_list[[1]][data_list[[1]]$sample%in%rep_samples_sel[1],]
    raw_rep2=data_list[[1]][data_list[[1]]$sample%in%rep_samples_sel[2],]
    norm_rep1=data_list[[2]][data_list[[1]]$sample%in%rep_samples_sel[1],]
    norm_rep2=data_list[[2]][data_list[[1]]$sample%in%rep_samples_sel[2],]
    raw_rep=list(raw_rep1,raw_rep2)
    norm_rep=list(norm_rep1,norm_rep2)
    EMD_metric_comp_raw_all_samples[ind,]=EMD_metric(raw_rep)
    EMD_metric_comp_norm_all_samples[ind,]=EMD_metric(norm_rep)
    ind=ind+1
  }
  saveRDS(EMD_metric_comp_raw_all_samples, paste(outdir,"EMD_metric_comp_raw_all_samples_all_cells.rds",sep=""))
  saveRDS(EMD_metric_comp_norm_all_samples, paste(outdir,"EMD_metric_comp_norm_all_samples_all_cells.rds",sep=""))
  EMD_all=rbind(EMD_metric_comp_raw_all_samples,EMD_metric_comp_norm_all_samples)
  ref=c(rep("raw",nsample/2),rep("norm",nsample/2))
  df=data.frame(EMD_all,ref,sample=rownames(EMD_all))
  ggdf <- melt(df, id.var = c("sample","ref"),
               value.name = "expression", variable.name = "antigen")
  #ggplot(ggdf,aes(x=expression,y=sample,color=ref))+geom_point(aes(colour = ref))+facet_wrap(~ref)
  ggplot(ggdf,aes(y=expression,x=sample,color=ref))+geom_boxplot(aes(colour = ref)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste(outdir,"EMD_metric_comp_raw_and_norm_all_samples_all_cells.png",sep=""))

  # ############################ 3. Clusters proportions accross rep
  # Compute different distances
  nclust=dim(props_table_norm)[1]
  nsample=dim(props_table_norm)[2]

  compute_dist <- function (props_table){
    nsample=dim(props_table)[2]
    nclust=dim(props_table)[1]
    Dist_props_table=matrix(nrow=5,ncol=nsample/2)
    colnames(Dist_props_table)=substring(colnames(props_table)[seq(1,nsample,2)],2)
    row.names(Dist_props_table)=c("Hell","TV","KL","BC","CW1")
    ind=1
    for (s in seq(1,nsample,2)){
      Dist_props_table[1,ind]=1/sqrt(2)*sqrt(sum((sqrt(props_table[,s])-sqrt(props_table[,s+1]))^2))
      Dist_props_table[2,ind]=1/2*sum(abs(props_table[,s]-props_table[,s+1]))
      Dist_props_table[3,ind]=sum(props_table[,s]*log(props_table[,s]/props_table[,s+1]))
      Dist_props_table[4,ind]=sum(sqrt(props_table[,s]*props_table[,s+1]))
      Dist_props_table[5,ind] = 1 - 4*sum(props_table[,s+1]/(props_table[,s]+props_table[,s+1])*
                                            (1-props_table[,s+1]/(props_table[,1]+props_table[,s+1])))/nclust
      ind=ind+1
    }
    return(Dist_props_table)
  }
  Dist_props_table_raw=compute_dist(props_table_raw)
  Dist_props_table_norm=compute_dist(props_table_norm)
  saveRDS(Dist_props_table_raw,paste(outdir,"Dist_props_table_raw_all_samples.rds",sep=""))
  saveRDS(Dist_props_table_norm,paste(outdir,"Dist_props_table_norm_all_samples.rds",sep=""))
}

#################### Example of computing the metrics comparing Raw data vs Norm data #################

# TODO: Create a new function starting here to plot!
# TODO: Create three functions for plotting - EMD, Hellinger, Sihlouette!

########################### Raw data  ##########################
wd_data="/stornext/Bioinf/data/lab_speed/Marie/CytofRUV_Figures/"
metadata_filename="Metadata.xlsx"
panel_filename="Panel.xlsx"
seed=1234
clusters_nb=20
#Loading the data
data=CytofRUV::load_data(wd_data,metadata_filename,panel_filename)

## Cluster the data
set.seed(seed)
data$daf=cluster_data(data$daf,seed,markers_to_use=data$lineage_markers,clusters_nb)
saveRDS(data,paste(wd_data,"Raw_data_clustered.rds",sep=""))

## Raw data
raw_data <- data.frame(sample = data$daf$sample_id, cluster=cluster_ids(data$daf,"meta20"), t(SummarizedExperiment::assay(data$daf)))
saveRDS(raw_data,paste(wd_data,"Raw_Data.rds",sep=""))

## Proportions raw data
md=data$md
samples_order=md$sample_id[order(md$patient_id)]
raw_data$sample=factor(raw_data$sample,levels=samples_order)
counts_table <- table(raw_data$cluster, raw_data$sample)
props_table_raw <- t(t(counts_table) / colSums(counts_table)) * 100
saveRDS(props_table_raw ,paste(wd_data,"props_table_raw.rds",sep=""))


########################### Norm data  ##########################

## Define parameters for normalisation
data=readRDS(paste(wd_data,"Raw_data_clustered.rds",sep=""))
daf=data$daf
md=data$md
dir_name_norm_data="CytofRUV_Norm_data_CLL2_HC1_all_cl_20_k5/"
raw_data=readRDS(paste(wd_data,"Raw_Data.rds",sep=""))
colnames(raw_data) <- gsub("^X", "",  colnames(raw_data))
rep_samples=list(c("CLL2_B1","CLL2_B2"),c("HC1_B1","HC1_B2"))
cluster_list_rep_samples <- list(seq(1,20),seq(1,20))
k_value <- 5
seed=1234

## CytofRUV normalisation
normalise_data(data=data,raw_data=raw_data,rep_samples=rep_samples, norm_clusters=cluster_list_rep_samples, k=k_value, num_clusters=clusters_nb,wd_data=wd_data,dir_norm_data=dir_name_norm_data)

## Define parameters to load and cluster the data
wd_norm=paste(wd_data,dir_name_norm_data,sep="")
metadata_norm_filename="Norm_Metadata.xlsx"
panel_norm_filename="Norm_Panel.xlsx"
seed=1234
clusters_nb=20

## Loading the norm data
norm_data=load_data(wd_norm,metadata_norm_filename,panel_norm_filename,cofactor=NULL)

## Cluster the norm data
set.seed(seed)
norm_data$daf=cluster_data(norm_data$daf,seed,markers_to_use=norm_data$lineage_markers,clusters_nb)
data <- data.frame(sample = norm_data$daf$sample_id, cluster=cluster_ids(norm_data$daf,"meta20"), t(SummarizedExperiment::assay(norm_data$daf)))
saveRDS(data,paste(wd_norm,"Norm_Data.rds",sep=""))

## Props table norm data
md=norm_data$md
samples_order=md$sample_id[order(md$patient_id)]
data$sample=factor(data$sample,levels=samples_order)
counts_table_norm <- table(data$cluster, data$sample)
props_table_norm <- t(t(counts_table_norm) / colSums(counts_table_norm)) * 100
saveRDS(props_table_norm ,paste(wd_norm,"props_table_norm.rds",sep=""))


################################# Parameters to compute metrics ###############################

## Raw files
raw_files=paste(wd_data,"Raw_Data.rds",sep="")
props_table_raw=readRDS(paste(wd_data,"props_table_raw.rds",sep=""))
## Norm data k3
wd_norm=paste(wd_data,dir_name_norm_data,sep="")
props_table_norm=readRDS(paste(wd_norm,"props_table_norm.rds",sep=""))
## Define parameters and sample to compute the silhouette and LDA
nb_cells=10000
rep_samples=c("CLL2_B1","CLL2_B2")
## On specific cluster
cluster_to_select_raw_data=c(9,2)
cluster_to_select_norm_data=c(17,14)

run_metrics(raw_files,cluster_to_select_raw_data,paste(wd_norm,"Norm_Data.rds",sep=""),wd_norm,cluster_to_select_norm_data,rep_samples,nb_cells,props_table_raw,props_table_norm)


############################### Plotting the metrics ###########################
library(ggplot2)
library(reshape2)
cols <- c("raw"="black","batchadj"="indianred4","cytonorm"="darkmagenta","cytofruv_k1"="#084081","cytofruv_k3"="#0868AC","cytofruv_k5"="#2B8CBE","cytofruv_k7"="#4EB3D3","cytofruv_k10"="#7BCCC4","cytofruv_k12"= "#A8DDB5","cytofruv_k15"="darkgreen","cytofruv_k20"="seagreen","cytofruv_k25"="yellowgreen")


##################### EMD plot ##################

### Dataset from CLL+HC samples with 20 clusters

# TODO: Change outdir
outdir="/stornext/Home/data/allstaff/t/trussart.m/R_scripts/RUV_design_experiment/RUV1b_RUV3b/Metrics_of_performances_all_methods/All_samples_DG23_VDBR996_anchor_CytofPackage_new_trial/"
dir="/stornext/Bioinf/data/lab_speed/Marie/RUV_experiment_design/Norm_Raw_RUV1b_RUV3b_simult/CytofRUV/"
### Loading EMD for each method
# Vector of k values to use
k=c(5,10,15,20,25)
EMD=list()
method=c("raw","cytofruv_k5","cytofruv_k10","cytofruv_k15","cytofruv_k20","cytofruv_k25","cytonorm","batchadjust")
## Reading the raw file EMD
# TODO: Make explicit that file names should not be changed! Only directories can be changed
EMD_raw_file=readRDS(paste(dir,"CytofRUV_Norm_data_CLL2_HC1_all_cl_20_k5/EMD_metric_comp_raw_all_samples_all_cells.rds",sep=""))
EMD_raw=data.frame(patient=row.names(EMD_raw_file),EMD_raw_file)
colnames(EMD_raw) <- gsub("^X", "",  colnames(EMD_raw))
EMD[["raw"]]=EMD_raw
## Reading the cytonorm file EMD
EMD_cytofnorm_file=readRDS("/stornext/Home/data/allstaff/t/trussart.m/R_scripts/RUV_design_experiment/RUV1b_RUV3b/Cytofnorm_Cytof_Package/EMD_metric_comp_norm_all_samples_all_cells.rds")
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


##################### Hellinger distance plot ##################

### Only CLL samples
# TODO: Change outdir
outdir="/stornext/Home/data/allstaff/t/trussart.m/R_scripts/RUV_design_experiment/RUV1b_RUV3b/Metrics_of_performances_all_methods/All_samples_DG23_VDBR996_anchor_CytofPackage_new_trial/"
dir="/stornext/Bioinf/data/lab_speed/Marie/RUV_experiment_design/Norm_Raw_RUV1b_RUV3b_simult/CytofRUV/"
### Loading Hellinger for each method
# Vector of k values to use
k=c(5,10,15,20,25)
Dist=list()
method=c("raw","cytofruv_k5","cytofruv_k10","cytofruv_k15","cytofruv_k20","cytofruv_k25","cytonorm","batchadjust")
## Reading the raw file Hellinger
Dist_raw_file=readRDS(paste(dir,"CytofRUV_Norm_data_CLL2_HC1_all_cl_20_k5/Dist_props_table_raw_all_samples.rds",sep=""))
Dist_raw=data.frame(patient=colnames(Dist_raw_file),Hell=Dist_raw_file[1,])
Dist[["raw"]]=Dist_raw
## Reading the cytonorm file Hellinger
Dist_cytofnorm_file=readRDS("/stornext/Home/data/allstaff/t/trussart.m/R_scripts/RUV_design_experiment/RUV1b_RUV3b/Cytofnorm_Cytof_Package/Dist_props_table_norm_all_samples.rds")
Dist_cytofnorm=data.frame(patient=colnames(Dist_cytofnorm_file),Hell=Dist_cytofnorm_file[1,])
Dist[['cytonorm']]=Dist_cytofnorm
## Reading the batchadjust file Hellinger
Dist_batchadj_file=readRDS("/stornext/Bioinf/data/lab_speed/Marie/RUV_experiment_design/Norm_Raw_RUV1b_RUV3b_simult/BatchAdjust/BatchAdjust_95p_and_DG23_new_trial/Dist_props_table_norm_all_samples.rds")
Dist_batchadj=data.frame(patient=colnames(Dist_batchadj_file),Hell=Dist_batchadj_file[1,])
Dist[["batchadj"]]=Dist_batchadj
## Reading the cytofruv files Hellinger
for (i in (1:length(k)) ){
  wd=paste(dir,"CytofRUV_Norm_data_CLL2_HC1_all_cl_20_k",k[i],"/",sep="")
  Dist_metric=readRDS(paste(wd,"Dist_props_table_norm_all_samples.rds",sep=""))
  Dist_tmp=data.frame(patient=colnames(Dist_metric),Hell=Dist_metric[1,])
  Dist[[method[i+1]]]=Dist_tmp
}
## Merging all the data
Dist_all=melt(Dist,id.var = c("patient"))
Dist_all$method=factor(Dist_all$L1,levels=c("raw","batchadj","cytonorm","cytofruv_k1","cytofruv_k3","cytofruv_k5","cytofruv_k7","cytofruv_k10","cytofruv_k12","cytofruv_k15","cytofruv_k20","cytofruv_k25"))
## Plotting CLL samples
ggdf2=Dist_all[Dist_all$patient%in%c("LL1_B1","LL2_B1","LL3_B1"),]
ggdf3=ggdf2[ggdf2$method %in% c("raw","batchadj","cytonorm","cytofruv_k5","cytofruv_k10","cytofruv_k15"),]
ggdf3$samp=""
ggdf3$samp[ggdf3$patient=="LL1_B1"]="CLL1"
ggdf3$samp[ggdf3$patient=="LL2_B1"]="CLL2"
ggdf3$samp[ggdf3$patient=="LL3_B1"]="CLL3"
ggdf3$distance=ggdf3$value
gg=ggplot(ggdf3,aes(y=distance,x=method,color=method))+ geom_point(aes(colour = method),size=3) + scale_color_manual(values = cols)#+ scale_color_manual(values=c("black","indianred4","darkmagenta", "#084081","#0868AC", "#2B8CBE", "#4EB3D3", "#7BCCC4", "#A8DDB5" ,"darkgreen","seagreen"))
gg+  theme_bw() + theme(aspect.ratio=1) +
  theme(text = element_text(size=40),
        axis.text.x = element_text(angle=90, hjust=1)) + theme(panel.background = element_blank()) +
  facet_wrap(~samp) +
  theme(axis.text.x = element_blank())



##################### Silhouette plot ##################

### Only CLL samples
# TODO: Change outdir
outdir="/stornext/Home/data/allstaff/t/trussart.m/R_scripts/RUV_design_experiment/RUV1b_RUV3b/Metrics_of_performances_all_methods/All_samples_DG23_VDBR996_anchor_CytofPackage_new_trial/"
dir="/stornext/Bioinf/data/lab_speed/Marie/RUV_experiment_design/Norm_Raw_RUV1b_RUV3b_simult/CytofRUV/"
### Loading Hellinger for each method
# Vector of k values to use
k=c(5,10,15,20,25)
CLL=list()
patient=c("DG23","DG33","DG27")
## Reading the raw file Silhouette for all samples
silh_finck=readRDS(paste(dir,"CytofRUV_Norm_data_CLL2_HC1_all_cl_20_k5/Sil_all_samples_all_clust_raw_data_over_ss.rds",sep=""))
CLL[[patient[1]]]=data.frame()
CLL[[patient[1]]][1,'method']="raw"
CLL[[patient[1]]][1,'bio']=as.numeric(silh_finck[1,2])
CLL[[patient[1]]][1,'batch']=as.numeric(silh_finck[1,3])
for (p in 2:length(patient)){
  CLL[[patient[p]]]=CLL[[1]]
  CLL[[patient[p]]][1,'bio']=as.numeric(silh_finck[p,2])
  CLL[[patient[p]]][1,'batch']=as.numeric(silh_finck[p,3])
}
## Reading the cytofnorm file Silhouette
cytofnorm=readRDS("/stornext/Home/data/allstaff/t/trussart.m/R_scripts/RUV_design_experiment/RUV1b_RUV3b/Cytofnorm_Cytof_Package/Sil_all_samples_all_clust_norm_data_over_ss.rds")
## Reading the batchadjust file Silhouette
batchadj=readRDS("/stornext/Bioinf/data/lab_speed/Marie/RUV_experiment_design/Norm_Raw_RUV1b_RUV3b_simult/BatchAdjust/BatchAdjust_95p_and_DG23_new_trial/Sil_all_samples_all_clust_norm_data_over_ss.rds")#,sep=""))
## Reading the Silhouette for all samples
for (p in c(1:3)){
  for (i in (1:length(k)) ){
    wd=paste(dir,"CytofRUV_Norm_data_CLL2_HC1_all_cl_20_k",k[i],"/",sep="")
    silh=readRDS(paste(wd,"Sil_all_samples_all_clust_norm_data_over_ss.rds",sep=""))
    CLL[[patient[p]]][i+1,'method']=paste("cytofruv_k",k[i],sep="")
    CLL[[patient[p]]][i+1,'bio']=as.numeric(silh[p,2])
    CLL[[patient[p]]][i+1,'batch']=as.numeric(silh[p,3])
  }
  ### Cytofnorm
  CLL[[patient[p]]][(length(k)+2),'method']="cytonorm"
  CLL[[patient[p]]][(length(k)+2),'bio']=as.numeric(cytofnorm[p,2])
  CLL[[patient[p]]][(length(k)+2),'batch']=as.numeric(cytofnorm[p,3])
  # Batchadjust
  CLL[[patient[p]]][(length(k)+3),'method']="batchadj"
  CLL[[patient[p]]][(length(k)+3),'bio']=as.numeric(batchadj[p,2])
  CLL[[patient[p]]][(length(k)+3),'batch']=as.numeric(batchadj[p,3])
  # ## Formating
  CLL[[patient[p]]]$method=factor(CLL[[patient[p]]]$method,levels=c("raw","batchadj","cytonorm","cytofruv_k1","cytofruv_k3","cytofruv_k5","cytofruv_k7","cytofruv_k10","cytofruv_k12","cytofruv_k15","cytofruv_k20","cytofruv_k25"))
  CLL[[patient[p]]]$bio=as.numeric(CLL[[patient[p]]]$bio)
  CLL[[patient[p]]]$batch=as.numeric(CLL[[patient[p]]]$batch)
}
## Merging all the data
All_patients_CLL=melt(CLL,id.var=c("bio","batch","method"))
## Plotting CLL samples
All_patients_CLL_select=All_patients_CLL[All_patients_CLL$method%in%c("raw","batchadj","cytonorm","cytofruv_k5","cytofruv_k10","cytofruv_k15"),]
All_patients_CLL_select$L1[All_patients_CLL_select$L1=="DG33"]="CLL1"
All_patients_CLL_select$L1[All_patients_CLL_select$L1=="DG23"]="CLL2"
All_patients_CLL_select$L1[All_patients_CLL_select$L1=="DG27"]="CLL3"
gg=ggplot(data = All_patients_CLL_select, aes(x = bio, y = batch))+ geom_point(aes(colour = method)) + facet_wrap(~L1) + scale_color_manual(values = cols)#+ scale_color_manual(values=c("black","indianred4","darkmagenta", "#084081","#0868AC", "#2B8CBE", "#4EB3D3", "#7BCCC4", "#A8DDB5" ,"darkgreen","seagreen"))
gg+ theme(aspect.ratio=1)+theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size=20))
