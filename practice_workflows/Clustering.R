library(bluster)
library(ggplot2)
library(dplyr)
library(plyr)
library(SC3)
library(SingleCellExperiment)

# Some Important Variables
Donor <- sce_Gn$Donor; Sample <- sce_Gn$Sample; PCA_Outliers <- sce_Gn$Outliers_PCA
Sample_Donor <- paste0(sce_Gn$Donor, '_', sce_Gn$Sample)


# Data Frame for plotting 
explan_var <- data.frame(Donor, Sample, Sample_Donor, PCA_Outliers)


# Louvain Clustering - k=10 - TRY WITH K = sqrt(Number_Cells) - TO DO::VISUALISE GRAPH ----
jaccard_graph <- buildSNNGraph(sce_Gn, use.dimred='Horns_PCA', type='jaccard', k=sqrt(ncol(sce_Gn)))
louv_clust <- factor(igraph::cluster_louvain(jaccard_graph)$membership)


# TSNE Plot of Louvain Clustering ---- 
clust_df_L <- cbind(tsne_df, louv_clust, sc3_clust)
clust_L_gg <- ggplot(clust_df_L) + geom_point(aes(x=X1, y=X2, color=louv_clust), size=1, shape=16, alpha=0.5) + 
  theme_classic() +
  labs(title=paste0('TSNE: Visualising Clusters from Louvain Community Detection\n'), 
       x='TSNE_1\n',
       y='TSNE_2\n',
       subtitle = paste0('KNN: ', round(sqrt(ncol(sce_Gn))), ', Edges: Jaccard\n'),
       colour='Clusters') +
  theme(
    legend.title =element_text(size=10),
    legend.text=element_text(size=9),
    legend.direction = 'horizontal',
    legend.box = 'horizontal',
    legend.justification = c('bottom'),
    legend.box.just = 'bottom',
    legend.position='bottom'
  ) +
  guides(shape=FALSE, alpha=FALSE, color=guide_legend(nrow=1)) + 
  geom_vline(xintercept=0, linetype='dashed', colour='#777777', size=0.2) +
  geom_hline(yintercept=0, linetype='dashed', colour='#777777', size=0.2)


# UMAP Plot of Louvain Clustering TO DO : CLEAN UP CODE ---- 
df_L_UMAP <- cbind(umap_df, louv_clust, explan_var, sc3_clust)
UM_l <- list()

j <- 1
for (i in colnames(df_L_UMAP[,-2:-1])) {
  UM_l[[j]] <- ggplot(df_L_UMAP) + 
    geom_point(aes_string(x='X1', y='X2', color=i), size=1, shape=16, alpha=0.7) +
    theme_classic() +
    labs(title=paste0('UMAP: Visualising Clusters from Louvain: ', i), 
         x='UMAP_1\n',
         y='UMAP_2\n',
         subtitle =  paste0('KNN: ', round(sqrt(ncol(sce_Gn))), ', Edges: Jaccard\n'),
         colour=i) +
    theme(
      legend.title =element_text(size=10),
      legend.text=element_text(size=9),
      plot.title = element_text(face='plain'),
      plot.subtitle =  element_text(size=9)
    ) +
    guides(shape=FALSE, alpha=FALSE) + 
    geom_vline(xintercept=0, linetype='dashed', colour='#777777', size=0.2) +
    geom_hline(yintercept=0, linetype='dashed', colour='#777777', size=0.2)
    print(UM_l[[j]])
    j <- j + 1
}


# TSNE of Louvain & SC3 Clustering 
tsne_sc3 <- tsne_func(clust_df_L, colour = 'sc3_clust',  title='TSNE:SC3 Clustering', PCS=horns$n, 
                      perplexity=round(perplexity), iterations = maxiter)

# Combined Plot of UMAP & TSNE Visualisation of Clusters ----
clust_leg <- get_legend(clust_L_gg)
title_clust <- cowplot::ggdraw() + cowplot::draw_label('TSNE & UMAP : Visualising Louvain Clusters')
tsn_umap_gr <- cowplot::plot_grid(clust_L_gg + theme(plot.title = element_blank(),
                                                     legend.position = 'none'), 
                                  clust_L_ggUMP + theme(plot.title = element_blank(),
                                                        legend.position='none'))
tsne_UMAP_clus <- cowplot::plot_grid(tsn_umap_gr, clust_leg, ncol=1, rel_heights=c(2, 0.2)) # SAVE






