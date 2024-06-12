

#' clustering function
#' methods: seurat, kmeans, kmedoids, hclust, dbscan,simlr

#' import Rtsne
#' import umap
#' import cluster
#' import Seurat
#' import SIMLR
#' import dbscan






#suppressPackageStartupMessages(library(Seurat))
#suppressPackageStartupMessages(library(Rtsne))
#suppressPackageStartupMessages(library(umap))
#suppressPackageStartupMessages(library(cluster))

#library(mclust)
#library(SC3)
#library(SingleCellExperiment)
#library("devtools")
#install_github("BatzoglouLabSU/SIMLR", ref = 'master')

#library('dplyr')
#library(reshape2)
# devtools::install_github("yycunc/SAFEclustering")
#library(cidr)

# cal_cidr = function(data.mat,cell_type, data_type){
#   cidr.data <- scDataConstructor(as.matrix(data.mat), tagType="cpm")
#   cidr.data <- determineDropoutCandidates(cidr.data)
#   cidr.data <- wThreshold(cidr.data)
#   cidr.data <- scDissim(cidr.data)
#   cidr.data <- scPCA(cidr.data)
#   cidr.data <- nPC(cidr.data)
#   nPC.cidr <- cidr.data@nPC
#   cidr.data <- scCluster(cidr.data, nPC = nPC.cidr)
#   return(adjustedRandIndex(cidr.data@clusters, cell_type))
# }

setGeneric("clustering", function(expr,...) standardGeneric("clustering"))

setMethod("clustering", signature(expr="matrix"),
          function(expr,
                   method=c("seurat",  "hclust", "simlr", "kmedoids", "kmeans", "dbscan"),
                   seurat_dims,
                   seurat_resolution,
                   k_cluster,
                   min_pts,
                   verbose=FALSE
)
{
  method <- match.arg(method)

  print(method)
  if (nrow(expr) < 2)
    stop("Less than two rows in the input data matrix\n")



  if (method == "seurat") {
      if(verbose)
            cat("clustering data matrix using PCA+KNN performed by Seurat\n")
    return(cal_seurat(expr,
                      seurat_dims,
                      seurat_resolution))
  }

  if (method == "tsne_kmeans") {
      if(verbose)
            cat("clustering data matrix using tsne+kmeans\n")
    return(cal_tsne_kmeans(expr,k_cluster))
  }

  if (method == "umap_kmeans") {
      if(verbose)
            cat("clustering data matrix using umap+kmeans\n")
    return(cal_umap_kmeans(expr,k_cluster))
  }

  if (method == "hclust") {
      if(verbose)
            cat("clustering data matrix using hclust\n")
    return(cal_hclust(expr,k_cluster))
  }

  if (method == "simlr") {
      if(verbose)
	        cat("Clustering data matrix using SIMLR\n")
	  return(cal_simlr(expr, k_cluster))
  }

  if (method == "kmedoids") {
      if(verbose)
	      cat("Clustering data matrix using K-medoids\n")
      return(cal_kmedoids(expr, k_cluster))
  }

  if (method == "kmeans") {
      return(cal_kmeans(expr, k_cluster))
  }

  if (method == "dbscan") {
      DBSCAN_groups = cal_dbscan(expr, min_pts)
	  return(DBSCAN_groups)
  }

})




getTPM = function(mat){

  geneLength = mat[,4]
  mat = mat[,5:dim(mat)[2]]
  mat = mat/geneLength                 # #R语言用矩阵除以一个向量的时候是按列数，第一列的第一个除以向量的第一个，第一列的第二个除以向量的第二个。。。。。
  mat = t(t(mat)*1000000/colSums(mat))
  return(mat)
}

ARI = function(predict_idnex, true_index){
  return(adjustedRandIndex(predict_idnex, true_index))
}

cal_dbscan = function(data,eps="auto",min_pts,tol=0.01){
  library(dbscan)
  dist = dist(t(data))
  #automatic determination of eps (the "elbow" in the kNNdistplot)
  if(eps=="auto"){
  kNNdist = sort(kNNdist(dist,min_pts))
  i = seq(1,length(kNNdist),as.integer(0.01*length(kNNdist)))
  slope_prev = 100
  for(indx in i){
	slope = kNNdist[indx]/indx
	if(slope_prev>=slope-tol*slope){
		slope_prev = slope
	} else {
	    elbow = indx
		break
	}
  }
  eps = kNNdist[elbow]
  print(paste("Epsilon: ",eps))
  } else if(!is.numeric(eps)){
	stop("Please provide a value for eps or set it to \"auto\"")} else {eps=eps}

  kNNdistplot(dist,k=min_pts)
  #abline(h=eps,col="red")
  res = dbscan(dist,eps = eps,minPts = min_pts)
  return(res$cluster)
}

cal_simlr = function(data,
                     k_cluster,
					 cores_ratio=0.5){
  library("SIMLR")
  res = SIMLR(X=data,
              c=k_cluster,
			  cores.ratio = cores_ratio)
  return(as.numeric(res$y[["cluster"]]))
}

cal_seurat = function(data,
                      dims = 8,
                      resolution = 0.5){

  pbmc <- CreateSeuratObject(counts = data, min.cells = 0,min.features = 0,
                             project = "seurat.clust")
  #pbmc <- NormalizeData(pbmc, verbose = FALSE)
  pbmc <- FindVariableFeatures(pbmc,
                               verbose = FALSE)
  pbmc <- ScaleData(pbmc,
                    verbose = FALSE)
  pbmc <- RunPCA(pbmc, seed.use = 42,
                 verbose = FALSE)

  pbmc <- FindNeighbors(pbmc, reduction = 'pca',
                        dims = 1:dims,
                        verbose = FALSE)
  pbmc <- FindClusters(pbmc, random.seed = 42,
                       resolution = resolution,
                       verbose = FALSE)
  return(as.numeric(Idents(pbmc)))
}

cal_tsne_kmeans = function(data,
                           n_cluster,
                           tsne_dims = 3,
                           pca_dims = 20){
  library(Rtsne)
  tsne_output <- Rtsne(t(data),dims = tsne_dims,
                       initial_dims = pca_dims,
                       pca_center = F,
                       theta=0.1)
  tsne_kmeansOUTPUT <- kmeans(tsne_output$Y, n_cluster)
  return(tsne_kmeansOUTPUT$cluster)
}

cal_kmeans = function(data,
                       n_cluster){
   #print(dim(data))
   kmeans_out = kmeans(t(data), n_cluster)
   #print(unique(kmeans_out$cluster))
   return(kmeans_out$cluster)
}

cal_kmedoids = function(data,
                        n_cluster){
  library(cluster)
  res = pam(data, n_cluster)
	return(res$cluster)
}

cal_umap_kmeans = function(data,
                           n_cluster,
                           dims=3){
  library(umap)
  umap_output <- umap(t(data), dims = dims)
  umap_kmeansOUTPUT <- kmeans(umap_output$layout,n_cluster)
  #table(umap_kmeansOUTPUT$cluster)
  #adjustedRandIndex(umap_kmeansOUTPUT$cluster, cell_type)
  return(umap_kmeansOUTPUT$cluster)
}

cal_hclust = function(data,
                      n_cluster,
                      method = "ward.D"){
  d <- dist(t(data))
  hc <- hclust(d, method = method)
  label <- cutree(hc, k = n_cluster)
  #table(label)
  return(label)
}

cal_mclust = function(data,
                      n_cluster){
  m_clust <- Mclust(t(data), G=n_cluster)
  return(m_clust$classification)
}

cal_SC3 = function(data,
                   n_cluster,
                   n_cores=2){
  exp_cell_exprs <- SingleCellExperiment(assays = list(tmp = norm(data)))
  exp_cell_exprs <- sc3(exp_cell_exprs, ks = n_cluster,
                        biology = FALSE, n_cores = 1, rand_seed = 1)
}

norm <- function(data){
  data = (data-min(data))/(max(data)-min(data))
  return(data)
}

