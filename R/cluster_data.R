#' @rdname cluster_data
#' @title cluster_data
#'
#' is used to cluster all the samples in the study using FlowSOM clustering
#' and metaclustering into with ConsensusClusterPlus.
#'
#'
#' @param daf data single cell experiment
#' @param seed seed used to reproduce the results
#' @param markers_to_use markers to use for clustering
#' @param clusters_nb number of clusters
#'
#'
#'
#' @return daf data single cell experiment after clustering
#' @export

cluster_data <- function(daf,seed,markers_to_use,clusters_nb){
  set.seed(seed)
  daf_cluster <- CATALYST::cluster(daf, features = markers_to_use, xdim = 10, ydim = 10, maxK = clusters_nb, seed = seed)
  return(daf_cluster)
}
