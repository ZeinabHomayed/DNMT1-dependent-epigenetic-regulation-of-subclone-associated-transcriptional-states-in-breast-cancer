
#
library(magrittr)
library(numbat)
library(dplyr)
library(glue)
library(logger)
library(rtracklayer)
library(data.table)
library(Seurat)
library(gtools)
library(ggplot2)
library(ggpubr)
library(AUCell)
library(presto)
library(rstatix)

#Open numbat scripts
dir <- "~/numbat_scripts"
r_list <- list.files(dir)
for (i in 1:length(r_list)) {
  source(paste0(dir, "/", r_list[i]))
}

#Load data
seu <- readRDS("~/seurat_object.rds")
gtf_final <- readRDS("~/gtf_file_to_inferCNA_with_chr_X_and_Y.rds")
agg_exp_normal <- readRDS("~/aggregate expression of pal 2021 embo normal cells for inferCNA ref.rds")


##Construct gtf file
centro_pos <- data.frame(CHROM = 1:24, centro_pos = c(123.4, 93.9, 90.9, 50.0, 48.8, 59.8, 60.1, 45.2,
                                                      43.0, 39.8, 53.4, 35.5, 17.7, 17.2, 19.0, 36.8,
                                                      25.1, 18.5, 26.2, 28.1, 12.0, 15.0, 61.0, 10.4)*1e6)
gtf_arm <- merge(gtf_final, centro_pos, by = "CHROM", all = T)
#Add arm annotation to each gene
gtf_arm[["arm"]] <- ifelse(gtf_arm$gene_start < gtf_arm$centro_pos, "p", "q")
gtf_arm[["chr_arm"]] <- paste0(gtf_arm$CHROM, gtf_arm$arm)


##Construct reference file
ref_normal <- agg_exp_normal$RNA
ref_normal <- as.matrix(numbat:::scale_counts(as(ref_normal, "dgCMatrix")))


##Infer CNAs
exp_roll <- exp_hclust(seu@assays$RNA@counts, lambdas_ref = ref_normal, gtf = gtf_final)
saveRDS(exp_roll, "exp_hclust output from SUM159_VF.rds")

#Assign subclones
assign_subclone <- function(seu, exp_roll, gtf_chr_arm){
  #####
  ##Only for female patients
  non_y_genes <- gtf_chr_arm[gtf_chr_arm$CHROM != 24, "gene"]
  cna_mat <- exp_roll$gexp_roll_wide[,colnames(exp_roll$gexp_roll_wide) %in% unlist(as.vector(non_y_genes))]
  
  #####
  #Apply mean by chromosome arm for each cell
  loc_group <- gtf_chr_arm[match(colnames(cna_mat), gtf_chr_arm$gene), c("gene", "chr_arm")]
  cell_profile_by_chr_arm <- split(data.frame(t(cna_mat)), loc_group$chr_arm)
  cell_by_chr_arm <- lapply(cell_profile_by_chr_arm, colMeans)
  cell_by_chr_arm <- cell_by_chr_arm[mixedsort(names(cell_by_chr_arm))]
  cell_by_chr_arm <- t(do.call("rbind", cell_by_chr_arm))
  
  
  #####
  ##############Clustering
  #Cluster cells from hclust and merge clusters according to Tirosh technique
  clusters <- cutree(exp_roll$hc, k =8)
  
  #Split cna matrix by cluster
  cna_by_group <- split(data.frame(cna_mat), clusters)
  
  #Merge clusters with less than 10 cells
  while (any(lapply(cna_by_group, nrow)<10)){
    small_cluster <- lapply(cna_by_group, nrow) %>% unlist %>% which.min %>% names #Select smallest cluster
    cluster_profile <- do.call("rbind", lapply(cna_by_group, colMeans)) #Compute (or recompute) mean profile by cluster
    dist <- as.matrix(parallelDist::parDist(cluster_profile, method = "euclidean")) #Compute euclidean distance matrix between clusters
    new_loc <- names(sort(dist[small_cluster,])[2]) #Select the closest cluster
    cna_by_group[[new_loc]] <- rbind(cna_by_group[[new_loc]], cna_by_group[[small_cluster]]) #Merge the small cluster with the closest cluster
    cna_by_group <- cna_by_group[names(cna_by_group) %in% small_cluster == FALSE] #Remove former small cluster (now merged)
  }
  
  
  #Merge clusters with correlation higher than 0.95
  cluster_profile <- do.call("rbind", lapply(cna_by_group, colMeans)) #Compute (or recompute) mean profile by cluster
  cor_clusters <- cor(t(cluster_profile))
  cor_clusters[upper.tri(cor_clusters, diag = T)] <- 0
  
  while (any(cor_clusters>.95)) { #Repeat as long as there is a pair of clusters with correlation higher than 0.95
    cluster_names <- colnames(cor_clusters)
    cluster_pair <- as.vector(which(cor_clusters == max(cor_clusters), arr.ind = T)) #We search for the clusters pair with the highest correlation
    new_loc <- cluster_names[min(cluster_pair)]
    old_loc <- cluster_names[max(cluster_pair)]
    cna_by_group[[as.character(new_loc)]] <- rbind(cna_by_group[[as.character(new_loc)]], cna_by_group[[as.character(old_loc)]]) #Merge the highest correlated clusters
    cna_by_group <- cna_by_group[names(cna_by_group) %in% as.character(old_loc) == FALSE] #Remove former cluster unmerged
    cluster_profile <- do.call("rbind", lapply(cna_by_group, colMeans)) #Recompute mean of expression for the new clusters
    cor_clusters <- cor(t(cluster_profile)) #Correlation between mean profiles of new clusters
    cor_clusters[upper.tri(cor_clusters, diag = T)] <- 0
  }
  
  if (length(cna_by_group) == 1){
    return("No subclones")
  }else{
    #####
    ################Search for CNA with big spacing between clusters' median
    #Divide the cell_by_chr_arm matrix into precedently computed clusters
    new_clusters <- data.frame(cell = gsub("-", ".", do.call("c", lapply(cna_by_group, rownames))),
                               cluster = rep(names(cna_by_group), times = lapply(cna_by_group, nrow))) #Create an annotation file the new clusters
    new_clusters <- new_clusters[match(rownames(cell_by_chr_arm), new_clusters$cell),]
    cell_chr_by_cluster <- split(data.frame(cell_by_chr_arm), new_clusters$cluster)
    
    #Find CNAs of interest with diff between group with the lowest value and group with the highest
    median_chr_cluster <- do.call("rbind", lapply(cell_chr_by_cluster, function(x){apply(x, 2, median)}))
    cna_diff <- apply(median_chr_cluster, 2, function(x){abs(min(x)-max(x))>0.4})
    cna_of_interest <- names(cna_diff[cna_diff == T])
    
    if (length(cna_of_interest) == 0){
      return("No subclones")
    }else{
      #For CNAs of interest, find cutoff values by performing cutpointr on highest and lowest groups
      library(cutpointr)
      cp_list <- list()
      for (i in 1:length(cna_of_interest)) {
        cna <- cna_of_interest[i]
        max_group <- rownames(median_chr_cluster[which.max(median_chr_cluster[,cna]),,drop=F])
        min_group <- rownames(median_chr_cluster[which.min(median_chr_cluster[,cna]),,drop=F])
        data_for_cutpoint <- data.frame(value = c(cell_chr_by_cluster[[max_group]][,cna], cell_chr_by_cluster[[min_group]][,cna]),
                                        group = rep(c("max", "min"),
                                                    times = c(nrow(cell_chr_by_cluster[[max_group]]), nrow(cell_chr_by_cluster[[min_group]]))))
        cp_list[[i]] <- cutpointr(data_for_cutpoint, value, group)
      }
      names(cp_list) <- cna_of_interest
      
      #Second cutoff taking AUC of 0.95 to separate min and max group
      cp_list <- cp_list[unlist(lapply(cp_list, function(x){x$AUC>0.9}))]
      
      if (length(cp_list) == 0){
        return("No subclones")
      }else{
        #Assign chromosomal indel to each cluster
        cna_by_cluster_list <- list()
        for (i in 1:length(cp_list)) {
          cutoff_cna <- cp_list[[i]]$optimal_cutpoint #extract cutoff precedently computed
          chromosome_arm <- names(cp_list)[i]
          if (cutoff_cna > 0){ #if cutoff sup to 0, we consider the situation as neutral/amp situation (to improve in the next version)
            cna_by_cluster_list[[i]] <- do.call("c", lapply(cell_chr_by_cluster, function(x){
              if (sum(x[,chromosome_arm]>cutoff_cna)>nrow(x)/2){ #if more than half of the cells present an amp, the cluster is considered to be amp
                return("amp")
              } else {
                return("neutral")
              }
            }))
          } else { #if cutoff inf to 0, we consider the situation as neutral/del situation
            cna_by_cluster_list[[i]] <- do.call("c", lapply(cell_chr_by_cluster, function(x){
              if (sum(x[,chromosome_arm]<cutoff_cna)>nrow(x)/2){ #if more than half of the cells present a del, the cluster is considered to be del
                return("del")
              } else {
                return("neutral")
              }
            }))
          }
        }
        
        names(cna_by_cluster_list) <- names(cp_list)
        cna_by_cluster <- data.frame(do.call("cbind", cna_by_cluster_list), cluster = names(cell_chr_by_cluster))
        
        
        #####
        ###############Pattern assignment and final subclone determination
        #Merge clusters with the same set of indels
        cna_by_cluster[["subclone"]] <- 0
        patterns <- distinct(cna_by_cluster[, 1:(ncol(cna_by_cluster)-2), drop = F])
        
        for (i in 1:nrow(patterns)) { ###this doesn't work when several indels
          for (j in 1:nrow(cna_by_cluster)) {
            if (all(cna_by_cluster[j,1:(ncol(cna_by_cluster)-2)] == patterns[i,])){
              cna_by_cluster[j, "subclone"] <- i
            }
          }
        }
        
        #Assign subclone identity to cells
        final_subclones <- merge(new_clusters, cna_by_cluster, by = "cluster", all = T)
        final_subclones <- final_subclones[match(Cells(seu), gsub("\\.", "-", final_subclones$cell)),]
        
        return(final_subclones)
      }
    }
  }
}
subclones <- assign_subclone(seu, exp_roll = exp_roll, gtf_chr_arm = gtf_arm)
seu[["subclone"]] <- subclones$subclone

#Process and plot
seu <- SCTransform(seu)
seu <- RunPCA(seu)
seu <- RunUMAP(seu, dims = 1:10)
umap_clone_no_ccreg <- UMAPPlot(seu, group.by = "subclone") + ggtitle("No cell cycle regression")
seu <- SCTransform(seu, vars.to.regress = c("S.Score", "G2M.Score"))
seu <- RunPCA(seu)
seu <- RunUMAP(seu, dims = 1:10)
umap_clone_ccreg <- UMAPPlot(seu, group.by = "subclone") + ggtitle("Cell cycle regression")
cd_clone <- DotPlot(seu, features = c("THY1", "MCAM", "ITGB4", "MME"), group.by = "subclone", scale = F) + ggtitle("Surface markers")
ggarrange(ggarrange(umap_clone_no_ccreg, umap_clone_ccreg, common.legend = T), cd_clone, widths = c(1.2, 0.8))
ggsave("subclone analysis SUM159.png",
       width = 12, height = 3.5,
       dpi = 300)

plot_exp_roll3 = function(gexp_roll_wide, subclone_label, gtf, lim = 0.8, n_sample = 300, reverse = TRUE) {
  require(ggplot2)
  require(patchwork)
  
  gexp_norm_long = gexp_roll_wide %>% 
    as.data.frame() %>%
    tibble::rownames_to_column('cell') %>%
    reshape2::melt(id.var = 'cell', variable.name = 'gene', value.name = 'exp_rollmean') %>%
    left_join(gtf, by = 'gene') %>%
    arrange(CHROM, gene_start) %>%
    mutate(gene_index = as.integer(factor(gene, unique(gene))))
  
  cells = unique(gexp_norm_long$cell)
  
  cell_sample = sample(cells, min(n_sample, length(cells)), replace = FALSE)
  
  if (length(subclone_label)>1){
    subclone_label$cell <- gsub("\\.", "-", subclone_label$cell)
    cell_order = subclone_label %>% arrange(sample(cell)) %>% arrange(subclone) %>% pull(cell)
  }else{
    cell_order = sample(cells)
  }
  
  chrom_labeller <- function(chr){
    chr[chr %in% c(18, 21)] = ''
    return(chr)
  }
  
  p_heatmap = gexp_norm_long %>%
    filter(cell %in% cell_sample) %>%
    mutate(cell = factor(cell, cell_order)) %>%
    arrange(cell) %>%
    mutate(cell_index = as.integer(cell)) %>%
    ggplot() +
    geom_tile(aes(x = gene_index, y = cell, fill = exp_rollmean)) +
    # geom_rect(aes(xmin = gene_start, xmax = gene_end, ymin = cell_index - 0.5, ymax = cell_index + 0.5, fill = exp_rollmean)) +
    scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white', midpoint = 0, limits = c(-lim,lim), oob = scales::squish) +
    theme_void() +
    scale_x_discrete(expand = expansion(0)) +
    scale_y_discrete(expand = expansion(0)) +
    theme(
      axis.text.y = element_blank(),
      legend.position = 'top',
      panel.spacing = unit(0, 'mm'),
      panel.border = element_rect(size = 0.5, color = 'wheat4', fill = NA),
      strip.text.y = element_blank(),
      # axis.title.y = element_text(angle = 90),
      # axis.title.x = element_text()
    ) +
    facet_grid(~CHROM, scales = 'free', space = 'free', labeller = labeller(CHROM = chrom_labeller)) +
    guides(fill = guide_colorbar(title = 'Expression\nmagnitude')) +
    xlab('Gene index') +
    ylab('Cell')
  
  if (length(subclone_label)>1){
    subclone_label$subclone <- as.factor(subclone_label$subclone)
    subclone_label$cell <- gsub("\\.", "-", subclone_label$cell)
    
    p_label = subclone_label %>%
      filter(cell %in% cell_sample) %>%
      mutate(cell = factor(cell, cell_order)) %>%
      arrange(cell) %>%
      ggplot() +
      geom_tile(aes(x = 1, y = cell, fill = subclone)) +
      theme_void() +
      theme(
        axis.text = element_blank(),
        legend.position = 'top',
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.2, "cm")
      ) +
      guides(fill=guide_legend(ncol=2))
  }
  
  if (length(subclone_label)>1) {
    (p_label | p_heatmap) + plot_layout(widths = c(1,10))
  } else {
    (p_heatmap) + plot_layout(widths = c(10))
  }
}
p <- plot_exp_roll3(gexp_roll_wide = exp_roll$gexp_roll_wide, subclone_label = subclones, gtf = gtf_final, n_sample = 600)



