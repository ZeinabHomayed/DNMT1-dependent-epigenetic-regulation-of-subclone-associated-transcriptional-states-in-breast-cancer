##Load dependencies
library(Seurat)
library(magrittr)
library(pracma)
#edist without having to dowload scperturbR package
edist <- function(seurat_object, groupby = "perturbation", reduction = "pca", verbose = TRUE) {
  labels <- seurat_object[[]][[groupby]]
  groups <- unique(labels)
  emb <- Seurat::Embeddings(seurat_object, reduction = reduction)
  
  df <- setNames(data.frame(matrix(ncol = length(groups),
                                   nrow = length(groups)),
                            row.names = groups),
                 groups)
  if (verbose) {
    print("Computing E-test statistics for each group.")
  }
  completed_groups <- c()
  for (groupx in groups) {
    for (groupy in groups){
      if (groupy %in% completed_groups) {
        next
      }  # skip if already computed
      x <- as.matrix(emb)[labels == groupx, ]
      y <- as.matrix(emb)[labels == groupy, ]
      
      N <- nrow(x)
      M <- nrow(y)
      
      dist_xy <- rdist::cdist(x,y)
      dist_x <- rdist::pdist(x)
      dist_y <- rdist::pdist(y)
      
      ed <- 2 * (sum(dist_xy) / (N * M)) - (sum(dist_x) / (N * N)) - (sum(dist_y) / (M * M))
      
      df[groupx, groupy] <- df[groupy, groupx] <- ed
    }
    completed_groups <- c(completed_groups, groupx)
  }
  return(df)
}

#Load data
seu_list <- readRDS("~/list of processed seurat objects with subclone info batch2.rds")



## All tITH Goyal's analysis with subsampled tumors
#Tumor subsampled 400 cells
seu_list <- lapply(seu_list, function(x){
  return(x[unlist(lapply(x, function(y){
    ncol(y)>=400
  }))])
})

tITH_400c <- lapply(seu_list, function(x){lapply(x, function(y){
  #Subsample to keep randomly 400 cells
  y <- subset(y, cells=sample(Cells(y), min(400, ncol(y))))
  ##PC analysis
  y <- SCTransform(y, verbose = F)
  y <- RunPCA(y, npcs=min(100, (ncol(y)-1)), verbose = F)
  #Generate the same analysis with randomized matrix
  r_matrix <- t(apply(y@assays$RNA@counts, 1, sample))
  colnames(r_matrix) <- colnames(y)
  r_y <- CreateSeuratObject(counts = r_matrix)
  r_y <- SCTransform(r_y, verbose = F)
  r_y <- RunPCA(r_y, npcs = min(100, (ncol(y)-1)), verbose = F)
  var_explained <- y@reductions$pca@stdev^2/sum(apply(y@assays$SCT@scale.data, 1, var))
  r_var_explained <- r_y@reductions$pca@stdev^2/sum(apply(r_y@assays$SCT@scale.data, 1, var))
  pc_needed <- min(which((var_explained-r_var_explained)<0))-1 #Here is the max PC needed
  
  ##Clustering at several resolution
  y <- RunPCA(y, npcs=min(50, (ncol(y)-1)), verbose = F)
  y <- FindNeighbors(y, verbose = F)
  res <- seq(0.0125, 2.5, 0.0125)
  nb_clust <- c()
  dist <- list()
  for (i in 1:length(res)) {
    y <- FindClusters(y, resolution = res[i], verbose = F)
    nb_clust <- c(nb_clust, y[[]][[paste0("SCT_snn_res.", res[i])]] %>% unique %>% length)
    dist_mat <- edist(y, groupby = paste0("SCT_snn_res.", res[i]), verbose=F)
    dist[[i]] <- dist_mat[lower.tri(dist_mat)]
  }
  names(dist) <- nb_clust
  clust_res <- data.frame(resolution = res, nb_clusters = nb_clust)
  #Compute AUC of nb of clusters according to resolution
  clust_res_auc <- clust_res[clust_res$resolution>=0.1 & clust_res$resolution<=2,]
  auc <- trapz(clust_res_auc$resolution, clust_res_auc$nb_clusters)
  #Compute mean euclidean distance between clusters at several resolutions
  mean_dist <- c()
  for (i in c(6, 8, 10)) {
    mean_dist <- c(mean_dist, mean(dist[[as.character(i)]]))
  }
  mean_dist_df <- data.frame(nb_clusters=c(6, 8, 10), mean_dist=mean_dist)
  return(list(nb_clust_by_res=clust_res, AUC=auc, 
              pairwise_dist=dist, mean_e_dist=mean_dist_df,
              pc_needed=pc_needed))
})})



