#' @rdname normalise_data
#' @title normalise_data
#'
#' contains all the functions to normalise several Cytof datasets by removing the
#' unwanted variation between datasets arising from experimental artefacts.
#'
#' @param daf Datasets before normalisation
#' @param raw_data Raw data
#' @param rep_samples Replicated samples
#' @param norm_clusters Clusters to be normalised
#' @param k Dimension of the unwanted variation
#' @param num_clusters Total number of clusters
#' @param wd_data Path to the directory containing all raw fcs files, the metadata file
#' and the panel file
#' @param dir_norm_data Directory name containing all norm fcs files, the metadata file
#' and the panel file
#'
#' @return Normalised metadata file
#' @export
#'

normalise_data <- function(data,raw_data,rep_samples, norm_clusters, k, num_clusters,wd_data,dir_norm_data){

  # Normalise the cells
  norm_cells <- run_RUVIII(raw_data, norm_clusters,k,rep_samples)
  # Output to save files
  output_dir <- file.path(wd_data, dir_norm_data)
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  # Save data into fcs files
  #print(data$md)
  new_md=save_norm_files(data,norm_cells, data$fcs_raw, data$md,data$panel,output_dir,k)
  # Save Metadata file
  writexl::write_xlsx(new_md, path=file.path(output_dir,"Norm_Metadata.xlsx"),
                      col_names=TRUE, format_headers = TRUE)

  writexl::write_xlsx(data$panel, path=file.path(output_dir,"Norm_Panel.xlsx"),
                      col_names=TRUE, format_headers = TRUE)
  ## Save normalisation details
  norm_setup=list(rep_samples=rep_samples,cluster_to_norm=norm_clusters,k=k,metadata_norm=new_md)
  saveRDS(norm_setup,file =file.path(output_dir,"Norm_setup.rds"))
  ## Save norm data
  saveRDS(norm_cells,file = file.path(output_dir,"Norm_data.rds"))

}

save_norm_files<- function(data,norm_cells, fcs_raw, md,panel,output_dir,k){

  # All the file names
  file_ids <- rep(md$file_name, flowCore::fsApply(fcs_raw, nrow))
  #Number of files
  num_files <- length(data$fcs_raw)
  # Panel markers fullname
  all_fullname_markers=c("Time",c(panel$fcs_colname[panel$antigen%in%data$lineage_markers],
                                  panel$fcs_colname[panel$antigen%in%data$functional_markers]))

  fcs_raw_asinh=fsApply(fcs_raw,function(x,cofactor = 5){
    expr <- exprs(x)
    expr <- expr[,all_fullname_markers]
    expr[,all_fullname_markers[2:length(all_fullname_markers)]] <- asinh(expr[,all_fullname_markers[2:length(all_fullname_markers)]] / cofactor)
    exprs(x) <- expr
    x
  })

  # Extract how many cells per sample
  n_cells <- fsApply(fcs_raw_asinh, function(x) nrow(exprs(x)))
  # Calculate the correction vector - how much we should substract from the inds to line them up the start of the flowframe
  correction_vec <- c(0, cumsum(n_cells))

  # # Helper function which returns a modified flowframe with only inds indexed cells
  subset_flowframe <- function(flow_frame, inds){
    expr <- exprs(flow_frame)
    expr <- expr[inds, ]
    exprs(flow_frame) <- expr
    flow_frame
  }

  ### Save Norm data
  fcs_norm=fcs_raw_asinh
  new_md=md
  dir_new_md=md
  for (i in 1:num_files) {
    # The inds of the cells in the right file
    inds <- which(file_ids == md$file_name[i])
    # Line the inds up with the flowframe
    if (length(inds)>1){
      corrected_inds <- inds - correction_vec[i]
      dir_new_md$file_name[i] <- file.path(output_dir,paste("Norm_RUVIII_k",k,"_",md$file_name[i],sep=""))
      new_md$file_name[i] <- file.path(paste("Norm_RUVIII_k",k,"_",md$file_name[i],sep=""))
      tmp_exp=exprs(fcs_norm[[i]])
      tmp_exp[,all_fullname_markers[2:length(all_fullname_markers)]]= norm_cells[inds,c(data$lineage_markers,data$functional_markers)]
      exprs(fcs_norm[[i]])=tmp_exp
      write.FCS(subset_flowframe(fcs_norm[[i]], corrected_inds), dir_new_md$file_name[i])
    }
  }
  return(new_md)
}

### Functions
make_residual_mat_several_rep<- function(data,clusters, norm_clus_list,samples,rep_samples_list){
  #make_residual_mat(raw_Y,data$cluster, norm_clusters,norm_clusters_second,data$sample,rep_samples,second_rep_samples)
  res_mat=data
  res_mat[]=0
  for (r in 1:length(rep_samples_list)){
    norm_clus=norm_clus_list[[r]]
    rep_samples=rep_samples_list[[r]]
    mean_pseud=matrix(nrow = length(norm_clus),ncol=dim(data)[2])
    for (i in 1:length(norm_clus)){
      tmp=((clusters == norm_clus[i])&(samples%in%rep_samples))
      mean_pseud[i,]=colMeans(data[tmp,])
      res_mat[tmp,]<- t(apply(data[tmp,], 1, function(x) x-mean_pseud[i,]))
    }
  }
  return(res_mat)
}


fastRUVIII = function(Y, M, ctl,res_mat,k=NULL, eta=NULL, average=FALSE, fullalpha=NULL){
  # Assumes good input
  if (!(k > 0)) stop("Bad input - read the documentation")
  Y = ruv::RUV1(Y,eta,ctl)
  m = nrow(Y)
  #Y0 = fast_residop(Y, M)
  Y0=res_mat
  fullalpha = diag(rsvd::rsvd(Y0)$d) %*% t(rsvd::rsvd(Y0)$v)
  alpha = fullalpha[1:k,,drop=FALSE]
  ac = alpha[,ctl,drop=FALSE]
  W = Y[,ctl] %*% t(ac) %*% solve(ac %*% t(ac))
  newY = Y - W %*% alpha
  return(list(newY = newY, fullalpha=fullalpha))
}


fast_residop <- function (A, B) {
  return(A - B %*% solve(t(B) %*% B) %*% (t(B) %*% A))
}


run_RUVIII <- function(data, norm_clusters, k,rep_samples){
  raw_Y <- as.matrix(data[3:ncol(data)])
  # Standardise the input and then compensate output
  col_means <- colMeans(raw_Y)
  col_sds <- apply(raw_Y, 2, function(x) sd(x))

  for(i in 1:ncol(raw_Y)){
    raw_Y[,i] <- (raw_Y[,i] - col_means[i])/col_sds[i]
  }
  # Run the actual RUVIII
  res_mat<-make_residual_mat_several_rep(raw_Y,data$cluster, norm_clusters,data$sample,rep_samples)
  norm_Y <- fastRUVIII(Y = raw_Y, M, ctl = c(1:ncol(raw_Y)), res_mat=res_mat,k = k)$newY

  for(i in 1:ncol(norm_Y)){
    norm_Y[,i] <- norm_Y[,i]*col_sds[i] + col_means[i]
  }

  return(norm_Y)
}

