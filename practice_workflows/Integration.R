library(batchelor)
library(Seurat)
library(scran)
library(scater)
library(SC3)
library(SingleCellExperiment)
library(bluster)
library(ggplot2)


# Integrate Across Donors - MNN ---- 
don_correct <-  fastMNN(sce_Gn, batch=sce_Gn$Donor, d = horns$n, k=20, subset.row=hvg.3000, 
                        BSPARAM=BiocSingular::IrlbaParam())
batch <- don_correct$batch

# UMAP & TSNE on Corrected Dim Red 
set.seed(1)
UMAP_int <- calculateUMAP(don_correct, dimred='corrected', n_epochs=1500, n_neighbors=perplexity,  min_dist=0.1, init='random')
set.seed(1)
TSNE_1_int <- calculateTSNE(don_correct, ntop=FALSE, perplexity=perplexity, max_iter=maxiter, dimred='corrected')

# Louvain Clustering - k=10 - TRY WITH K = sqrt(Number_Cells) - VISUALISE GRAPH ----
jac_int <- buildSNNGraph(don_correct, use.dimred='corrected', type='jaccard', k=sqrt(ncol(don_correct)))
louv_int <- factor(igraph::cluster_louvain(jac_int)$membership)

# SC3 Clustering ON CORRECTED ----
# Run SC3
sce_sc3 <- sce_Gn
rowData(sce_sc3)$feature_symbol <- rownames(sce_sc3)
logcounts(sce_sc3) <- as.matrix(logcounts(sce_sc3))
sce_sc3 <- sc3(sce_sc3, ks=8,  n_cores=1, biology=TRUE, gene_filter = FALSE) # Same clustes for Comp ----
sc3_clust <- colData(sce_sc3)$sc3_8_clusters
saveRDS(sce_sc3, 'SC3_clustering')

# Extract Some SC3 Plots
sc3_calc_consens(sce_sc3)
# sc3_plot_consensus(sce_sc3, k=8) # TO DO ----
sc3_markers_pl <- sc3_plot_markers(sce_sc3, k=8) 
sc3_markers_pl
ggsave('Plots/Marker_gene_plots.png', sc3_markers_pl,height = 20, width=7)

# Pairwise Rand : SC3 & Louv FIX TO COMPARE CORRECTED ----
louv_sc3_rand <- pairwiseRand(sc3_clust, louv_int, mode='index')
# % of cells in each cluster in Louvain & SC3
louv_sc3 <- table(Louvain=louv_int, SC3=sc3_clust)
louv_sc3 <- as.data.frame(louv_sc3/rowSums(louv_sc3))

# Heatmap of % of cells in each cluster from louvain & SC3 ---- FIX COLOUR / TIDY VISUALISATION
louv_sc3_heat <- ggplot(louv_sc3) +
  geom_tile(aes(x=Louvain, y=SC3, fill=Freq), color='#777777') + 
  scale_fill_gradientn(colours=viridis::viridis(100)) +
  theme_classic() +
  labs(x='Louvain Clustering',
       y='SC3',
       fill='Cell \nProportion:',
       subtitle = 'Proportion of Cells\nbelonging to each Cluster',
       caption = paste0('Rand Index: ', round(louv_sc3_rand,3))) + 
  theme(plot.caption = element_text(hjust=0, size=11, face='bold'))

ggsave('Plots/Clusters/sil.png', silhoette_widths_sc3) # TO FIXXXX!!!!!

# Silhouette widths on Louv & SC3 
silhoette_widths_sc3 <- sc3_plot_silhouette(sce_sc3, k=8) # SC3
sil_louv <- as.data.frame(approxSilhouette(reducedDim(sce_Gn, 'Horns_PCA'), cluster=louv_int)) # Louv
sil_louv$closest <- factor(ifelse(sil_louv$width > 0, louv_int, sil_louv$other))
sil_louv$cluster <- factor(louv_int)

sil_louvgg <- ggplot(sil_louv, (aes(x=cluster, y=width, colour=closest))) +  
        ggbeeswarm::geom_quasirandom(method='smiley', width = 0.4, alpha=0.8, size=1.2) + 
        theme_classic() +
        labs(subtitle='Silhouette Widths: Louvain Cluster Separation',
             x='Cluster',
             y='Silhouette Width',
             colour='Nearest\nCluster') +
        theme() +
          geom_hline(yintercept=0, linetype='dashed', colour='#777777', size=0.2) +
        guides(shape=FALSE, alpha=FALSE, 
        color=guide_legend(override.aes = list(size=3, alpha=1)))

# Plot Both on Same Plot :: TO DO :: ----
ggsave('Plots/silhoette.png', device = 'png')

# Data Frame For Plotting TSNE & UMAP
int_df <- data.frame(UMAP=UMAP_int, TSNE=TSNE_1_int, louv_clusters=louv_int, batch,
                     sc3=sc3_clust)

# Visualise UMAP
umap_1_int <- ggplot(int_df) + 
  geom_point(aes(x=UMAP.1, y=UMAP.2, color=louv_clusters), size=1, shape=16, alpha=0.7) +
  theme_classic() +
  labs(title=paste0('UMAP: Integrated across Donors'), 
       x='UMAP_1\n',
       y='UMAP_2\n',
       subtitle = paste0('KNN: ', round(sqrt(ncol(don_correct))), ', Edges: Jaccard\n'),
       colour='Cluster') +

  theme(
    legend.title =element_text(size=10),
    legend.position=c(0.9, 0.2),
    legend.text=element_text(size=9),
    plot.title = element_text(face='plain'),
    plot.subtitle =  element_text(size=9)
  ) +
  guides(shape=FALSE, alpha=FALSE) + 
  geom_vline(xintercept=0, linetype='dashed', colour='#777777', size=0.2) +
  geom_hline(yintercept=0, linetype='dashed', colour='#777777', size=0.2)


# Visualise TSNE : TO DO:: LOOP OVER BATCHES ----
tsne_1_int_sc3 <- ggplot(int_df) +
  geom_point(aes(x=TSNE.1, y=TSNE.2, color=sc3), size=1, shape=16, alpha=0.5) + 
  theme_classic() +
  labs(subtitle='Single-Cell Consensus (SC3) Clustering', 
       x='TSNE_1\n',
       y='TSNE_2\n',
       colour='Cluster:',
       caption = 'Based on K-means') +
  theme(
    plot.caption = element_text(hjust=0),
    legend.title =element_text(size=10),
    legend.text=element_text(size=9),
    legend.direction = 'horizontal',
    legend.box = 'horizontal',
    legend.justification = c('bottom'),
    legend.box.just = 'bottom',
    legend.position='bottom'
  ) +
  guides(shape=FALSE, alpha=FALSE,
         color=guide_legend(override.aes = list(size=3, alpha=1.5), nrow=1)) + 
  geom_vline(xintercept=0, linetype='dashed', colour='#777777', size=0.2) +
  geom_hline(yintercept=0, linetype='dashed', colour='#777777', size=0.2)

tsne_1_int_louv <- ggplot(int_df) +
  geom_point(aes(x=TSNE.1, y=TSNE.2, color=louv_int), size=1, shape=16, alpha=0.5) + 
  theme_classic() +
  labs(subtitle='Single-Cell Consensus (SC3) Clustering', 
       x='TSNE_1\n',
       y='TSNE_2\n',
       colour='Cluster:',
       caption = 'Based on K-means') +
  theme(
    plot.caption = element_text(hjust=0),
    legend.title =element_text(size=10),
    legend.text=element_text(size=9),
    legend.direction = 'horizontal',
    legend.box = 'horizontal',
    legend.justification = c('bottom'),
    legend.box.just = 'bottom',
    legend.position='bottom'
  ) +
  guides(shape=FALSE, alpha=FALSE,
         color=guide_legend(override.aes = list(size=3, alpha=1.5), nrow=1)) + 
  geom_vline(xintercept=0, linetype='dashed', colour='#777777', size=0.2) +

# Combined TSNE Plots of SC3 and Louvain Clustering ----
title_clust_both <- cowplot::ggdraw() + cowplot::draw_label('TSNE: Comparing Seurat* and SC3 Clustering')
subtitle_int <- cowplot::ggdraw() + cowplot::draw_label(paste0('Number of PCs:', horns$n, ', perplexity: ', 
                                                               round(perplexity), ', Iterations: ', maxiter, '\n'), size=10)
tsne_clust_gr <- cowplot::plot_grid(tsne_1_int_louv,
                                    tsne_1_int_sc3, nrow=2, ncol=1)
tsne_clust_com <- cowplot::plot_grid(title_clust_both, subtitle_int, tsne_clust_gr, int_leg_clu, 
                                     ncol=1, rel_heights = c(0.1, 0.1, 2, 0.1), rel_widths=c(0.1,0.1,2,0.1))
ggsave('Plots/Clusters/SC3_Seurat.png', tsne_clust_com, width=10, height=10)

tsne_1_int_louv <- tsne_1_int_louv + coord_fixed(1) + theme(legend.position = 'none',
                                         plot.title = element_blank(),
                                         plot.subtitle  =element_blank(),
                                         plot.caption=element_blank())
tsne_1_int_sc3 <- tsne_1_int_sc3 + coord_fixed(1) + theme(legend.position = 'none',
                                        plot.title=element_blank(),
                                        plot.subtitle =  element_blank(),
                                        plot.caption=element_blank())

# Combined Plot - Evaluating Clustering ---- CAREFUL OF LEGEND
int_leg_clu <- get_legend(tsne_1_int_sc3)
title_clust_int <- cowplot::ggdraw() + cowplot::draw_label('Evaluating Seurat* and SC3 Clustering')
clust_com <- cowplot::plot_grid(scree + theme(legend.position=c('bottom')), 
                                sil_louvgg + theme(legend.position = 'none'),
                                louv_sc3_heat,
                                sc3_plot_silhouette(sce_sc3, k=8), nrow=2, ncol=2, 
                                rel_widths = c(0.7, 1.0))
tsne_clus <- cowplot::plot_grid(title_clust_int, clust_com, int_leg_clu, ncol=1, rel_heights=c(0.1,5, 0.2)) # SAVE

ggsave('Plots/Clusters/sil_widths.png', tsne_clus) # SAVED





  
