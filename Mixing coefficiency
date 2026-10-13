##Load dependencies
library(Seurat)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(magrittr)

#Load data
seu_list <- readRDS("list of processed seurat objects with subclone info batch2.rds")

#Keep only tumors with several subclones
seu_oligo_list <- lapply(seu_list, function(x){
  return(x[unlist(lapply(x, function(y){y$subclone[1]!="No subclones"}))])
})

#Compute average spearman correlation compared to random and mixing coefficients of tumor subclones
homogeneity_score_list <- lapply(seu_oligo_list, function(x){lapply(x, function(y){
  #Subset seurat object to have the same number of cells by subclone
  min_nb <- min(table(y$subclone))
  subclone_meta <- data.frame(cell_id=Cells(y), subclone=y$subclone)
  keep_cells <- subclone_meta %>% arrange(sample(1:nrow(subclone_meta))) %>%
    group_by(subclone) %>% slice_sample(n=min_nb) %>% pull(cell_id)
  
  y_sub <- subset(y, cells=keep_cells)
  y_sub <- SCTransform(y_sub)
  y_sub <- RunPCA(y_sub, npcs=10)
  y_sub <- FindNeighbors(y_sub, k.param=10)
  
  
  clone <- unique(y_sub$subclone)
  #Compute correlation matrix of pc embeddings
  pc_emb <- Embeddings(y_sub)
  cor_mat <- cor(t(pc_emb))
  #Extract nearest-neighbors graph
  nn_graph <- y_sub@graphs[["SCT_nn"]]
  
  #Compute scores
  scores <- data.frame()
  for (i in 1:length(clone)) {
    ##Compute spearman corr clone homogeneity
    cells_clone <- Cells(y_sub)[y_sub$subclone == clone[i]]
    cor_mat_sub <- cor_mat[cells_clone, cells_clone]
    mean_cor <- mean(cor_mat_sub[lower.tri(cor_mat_sub)])
    #Same with randomly chosen cells
    r_cells <- sample(Cells(y_sub), length(cells_clone))
    cor_mat_sub <- cor_mat[r_cells, r_cells]
    r_mean_cor <- mean(cor_mat_sub[lower.tri(cor_mat_sub)])
    distance_homogeneity_score <- mean_cor-r_mean_cor
    
    ##Compute nearest neighbors clone homogeneity (normalized mixing coefficient)
    sub_nn_graph <- nn_graph[cells_clone,cells_clone]
    connectivity_homogeneity_score <- median((10-rowSums(sub_nn_graph))/rowSums(sub_nn_graph))/(length(clone)-1)
    nb_self_nn <- median(rowSums(sub_nn_graph))
    
    scores <- rbind(scores, data.frame(subclone = clone[i],
                                       cor_score = distance_homogeneity_score,
                                       nn_score = connectivity_homogeneity_score,
                                       nb_self_nn = nb_self_nn, nb_subclones=length(clone)))
  }
  return(scores)
})})
