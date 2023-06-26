library(scran)
library(scater)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)
library(cowplot)



ggsave("Plots/sc3_markers.png", sc3_markers) # FIX FIX ---- 

# Added Names of 'POTENTIAL' Cells After 
sc3_markers <- c(Acinar='REG1A', Acinar='PRSS3P2', Islet='REG1B', Thyroid='TTR', Trypsin='PRSS1',
                Islet='REG3A',Pancreatitis='SPINK1', 
                Beta='INS', Beta='INS-IGF2',
                Alpha='GCG', 
                Exocrine='SPP1',
                Ductal='KRT19', 
                Pseudogene='PGM5P2', 
                Pancreatic_poly_peptide='PPY',
                Delta='UNC5B', Delta='FFAR4', Delta='SST',  Delta='GABRG2',
                Epsilon='GHRL',  Epsilon='ESAM'
               )

Delta <- c('SST', 'UNC5B', 'GABRB3', 'GABRG2', 'CASR', 'FFAR4', 'GPR120',
           'KCNJ2')
epsilon <- c('GHRL', 'ESAM')
# TO DO :: TEST THESE !!!!
# CHRM3 = Gamma / RGS4 = AKPHA / SST (delta??)
# DELTA: UNC5B, GABRB3, GABRG2, CASR, FFAR4/GPR120, and KCNJ2)
# GAMMA : SERTM1, SPOCK1, ABCC9, and SLITRK6
# NPY1R,OPRK1, PTGER4, and ASGR1, and processing enzymes such
# 'FTH1', 'FTL'

# Get Labels 
# TRIAL WITH LOUVAIN CLUSTERING CHECK ---- 
df_marker <- as.data.frame(TSNE_1_int) %>% 
                    mutate(clusters=louv_int) %>% 
                    mutate(ori_clusters=louv_clust)
df_clus_lab <- df_marker %>% 
                group_by(clusters)  %>% 
                summarise(x=mean(V1), y=mean(V2))
# UPDATED LABELS
df_clus_lab$clusters <- plyr::mapvalues(df_clus_lab$clusters,
                                        from=c(1,2,3,4,5,6,7,8),
                                        to=c('PGM5P2', 'Epsilon', 'Delta\nPP', 'SPP1', 
                                             'Acinar\nREG3A', 'Alpha', 'KRT19', 'Beta'))
# Cluster Numbers
df_clus_lab$numbers <- df_marker %>% 
  group_by(clusters)  %>% 
  summarise(x=mean(V1), y=mean(V2))

# TSNE for each Marker Gene ----
tsne_marker <- purrr::imap(sc3_markers, ~.y, function(marker){
  df_marker %>% 
    dplyr::mutate(gene=as.numeric(logcounts(sce_Gn[marker,]))) %>%
    ggplot(., aes(V1, V2)) + 
    theme_classic() + 
    labs(x='TSNE_1',
         y='TSNE_2',
         title=.y) + 
    geom_point(aes_string(color='gene'), alpha=0.8) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue") +
    theme(plot.title = element_text(face='plain', size=9)) +
    geom_text(data=df_clus_lab$numbers, 
              aes(label=clusters, x, y)) 
}) %>%
      cowplot::plot_grid(plotlist=.) 
marker_title <- cowplot::ggdraw() + cowplot::draw_label('Identifying Marker Genes')
tsne_marker_grid_test <- cowplot::plot_grid(marker_title, tsne_marker, ncol=1, rel_heights =c(0.1,2))

ggsave('Plots/Marker_gene_plots/marker_gene_tsne.png', tsne_marker_grid, width = 10, height = 10)



# FINAL UMAP PLOT
df_umap <- as.data.frame(int_df) %>% 
  mutate(clusters=louv_int) %>% 
  dplyr::select(UMAP.1, UMAP.2, clusters)
# Get coordinated for labels
cell_umap_lab <- df_umap %>% 
  group_by(clusters)  %>% 
  summarise(x=mean(UMAP.1), y=mean(UMAP.2))
# Replace Clusters with Cell Annotation
cell_umap_lab$clusters <- plyr::mapvalues(cell_umap_lab$clusters,
                                        from=c(1,2,3,4,5,6,7,8),
                                        to=c('PGM5P2', 'Epsilon', 'Delta\nPP', 'SPP1', 
                                             'Acinar\nREG3A', 'Alpha', 'KRT19', 'Beta'))

# FINAL PLOT ----
UMAP_FIN <- ggplot(df_umap) +
  geom_point(aes(x=UMAP.1, y=UMAP.2, colour=I(clusters)), size=0.8, shape=16, alpha=1/3) + 
  theme_classic() +
  coord_fixed(1) +
  labs(title='UMAP: Cell Type Annotation',
      subtitle='Seurat* Clustering on Integrated Data', 
       x='UMAP_1\n',
       y='UMAP_2\n',
       colour='Cluster:') +
  theme(
    plot.caption = element_text(hjust=0),
    legend.title =element_text(size=10),
    legend.text=element_text(size=9),
    legend.position=c('top', 'right'),
    plot.title = element_text(hjust=0.5),
    plot.subtitle = element_text(hjust=0.5)
  ) +
  guides(shape=FALSE, alpha=FALSE,
         color=guide_legend(override.aes = list(size=3, alpha=1.5))) + 
  geom_vline(xintercept=0, linetype='dashed', colour='#777777', size=0.2) +
  geom_hline(yintercept=0, linetype='dashed', colour='#777777', size=0.2)  + 
  geom_text_repel(data=cell_umap_lab, aes(label=clusters,x,y),
                box.padding = 0.4,max.overlaps=Inf, size=3, nudge_x = 0.2, nudge_y = -0.1)

ggsave('Plots/Clusters/UMAP_FINAL.png', UMAP_FIN)










