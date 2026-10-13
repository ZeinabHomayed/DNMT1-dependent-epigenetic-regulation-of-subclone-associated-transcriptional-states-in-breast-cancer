library(Seurat)
library(readxl)
library(HGNChelper)
library(dplyr)
library(magrittr)
source("gene_sets_prepare.R"); source("sctype_score_.R")

#Prepare list of positive and negative markers of immune populations
gs_list <- gene_sets_prepare("ScTypeDB_simple.xlsx", "Immune system")

#Load data seurat objects
gao1 <- readRDS(~/gao/Gao1_tot_filter.rds")
gao2 <- readRDS("~/gao/Gao2_tot_filter.rds")
gao <- merge(gao1, gao2)
qian <- readRDS(~/qian/Qian_tot_filter.rds")
wu <- readRDS("~/Wu_tot_filter.rds")
pal <- readRDS("~/Pal_wo_LN_filter.rds")
griffiths <- readRDS("~/inferCNV/griffiths/Griffiths10X_untreated_filter.rds")
bassez <- readRDS("~/bassez/Bassez_untreated_filter.rds")
xu <- readRDS("~/xu/Xu_tot_filter_man_annot.rds")

seu_list <- list(pal=pal, qian=qian, wu=wu, gao=gao,
                 xu=xu, griffiths=griffiths, bassez=bassez)


#Classic pre-processing
seu_list <- lapply(seu_list, function(x){
  x <- SCTransform(x, variable.features.n = 10000, conserve.memory=T)
  x <- RunPCA(x)
  x <- FindNeighbors(x)
  x <- FindClusters(x, resolution = 0.8)
})

#SCtype computation
for (i in 1:length(seu_list)) {
  es.max <- sctype_score(seu_list[[i]]@assays$SCT@scale.data, scaled = T,
                         gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
    
  cL_resutls = do.call("rbind", lapply(unique(seu_list[[i]]@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(seu_list[[i]]@meta.data[seu_list[[i]]@meta.data$seurat_clusters==cl, ])]),
                     decreasing = !0) #for each cluster, sort cell type signature by score from best to worst
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seu_list[[i]]@meta.data$seurat_clusters==cl)),
         10) #Create dataframe with 10 cell type signatures with best scores for the cluster with signature name, score, number of cells in the cluster and cluster ID
  })) #rbind the dataframes for every cluster
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores) #Select best scores for every cluster
  
  # set low-confident (low ScType score) clusters to "unknown" (if score is inf to the forth of the number of cells in the cluster)
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  
  seu_list[[i]]@meta.data$customclassif = ""
  for(k in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==k,]; 
    seu_list[[i]]@meta.data$customclassif[seu_list[[i]]@meta.data$seurat_clusters == k] = as.character(cl_type$type[1])
  } #Attribute celltype to the seurat object
}

meta_list <- lapply(seu_list, function(x){x@meta.data})
