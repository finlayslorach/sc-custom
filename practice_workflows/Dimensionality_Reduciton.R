library(PCAtools)
library(scran)
library(scater)
library(ggplot2)
library(dplyr)
library(gridExtra)


# Top 3000 HVG - REPEAT WITH HIGHER LOWER SIGNIFICANCE THRESHOLD ----
hvg.3000 <- rowData(sce_Gn)$sct.variable
horns <- PCAtools::parallelPCA(logcounts(sce_Gn)[hvg.3000,], threshold = 0.05, 
                               BSPARAM=BiocSingular::IrlbaParam())

# 12 PCs identified 
reducedDim(sce_Gn, 'Horns_PCA') <- horns$original$rotated[,1:20]

# Prepare Plot
# Mean of variance for each PC from each iteration 
permuted <- rowMeans(horns$permuted)
perm_var_df <- data.frame(permuted, original_var=horns$original$variance)
# Add Column for Plotting PCs (Numeric)
perm_var_df <- perm_var_df %>% 
            top_n(50) %>% 
            mutate(number=seq(1:50))

# Plot Of Variance of Permuted Matrix vs Original ----
scree <-ggplot(perm_var_df) + 
        geom_point(aes(x=number,y=original_var, group=1, shape='.', color='Explained by PCS'))+
        geom_line(aes(x=number,y=original_var))+
        geom_point(aes(x=number,y=permuted, group=1, shape='.', color='Explained by chance')) +
        geom_line(aes(x=number,y=permuted)) +
        theme_classic() +
        labs(x='Number of PCS',
             y='Variance Explained',
             subtitle = paste0('Scree plot: Permuted vs Original\nOptimum Number of PCs:', horns$n)) +
           theme(
            legend.title = element_blank(),
            legend.position=c(0.9, 0.9),
            legend.margin = margin(4, 4, 4, 4)
          ) 
          
scree <- scree + guides(shape=FALSE, alpha=FALSE, 
                 color=guide_legend(override.aes = list(size=3, alpha=1))) + 
      geom_vline(xintercept=horns$n, linetype='dashed', 
                 colour='#777777', size=1)
# SAVE

# Hyperparamters 
maxiter <- 1500
perplexity <- sqrt(ncol(sce_Gn))

# Run TSNE on Top 12 PCs - REPEAT AT MORE ITERATIONS? ----
set.seed(1)
TSNE_1 <- calculateTSNE(sce_Gn, ntop=FALSE, perplexity=perplexity, max_iter=maxiter, dimred='Horns_PCA')
tsne_df <- data.frame(TSNE_1, phases=assignments$scores$phase)



# Plot TSNE - CC & Split CC ----
tsne_1_CC <- tsne_func(tsne_df, colour='phases', title='Visualising Optimum Number of PCS', 
                          PCS=horns$n, perplexity=round(perplexity), iterations = maxiter)
tsne_CC_split <- tsne_func(tsne_df, 'phases' , 'TSNE: Cell Cycle Variation', PCS=horns$n, 
                           perplexity=round(perplexity), iterations=maxiter) + facet_wrap(~phases) + 
                            theme(legend.position = 'none')

# Save Graphs 
ggsave(tsne_1_CC, filename='Plots/Dim_red_plots/TSNE_CC_single.png') # SAVED
ggsave(tsne_1_CC, filename='Plots/Dim_red_plots/TSNE_CC_faceted.png') #SAVED


# Run UMAP on Top Features - Using Same Paramters as Used for tSNE 
set.seed(1)
UMAP_1 <- calculateUMAP(sce_Gn, dimred='Horns_PCA', n_epochs=1500, n_neighbors=perplexity,  min_dist=0.1, init='random')
umap_df <- data.frame(UMAP_1, phases=assignments$scores$phase)

# Plot UMAP ----
umap_1_CC <- ggplot(umap_df) + 
          geom_point(aes(x=X1, y=X2, color=phases), size=0.8, shape=16, alpha=0.7) +
          labs(title=paste0('UMAP: Cell Cycle Variation'), 
               x='UMAP_1\n',
               y='UMAP_2\n',
               subtitle = paste0('Number of PCs: ', horns$n, ', Number of Neighbours: ', round(perplexity), ', \nIterations: ', maxiter, ', Min Distance: 0.1\n'),
               colour='Cell Phase') +
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

# UMAP Faceted by Cell Cycle Variation ----
umap_1_CC_split <- ggplot(umap_df) + 
  geom_point(aes(x=X1, y=X2, color=phases), size=0.8, shape=16, alpha=0.7) +
  facet_wrap(~phases) +
  theme_classic() + 
  labs(title=paste0('UMAP: Cell Cycle Variation'), 
       x='UMAP_1\n',
       y='UMAP_2\n',
       subtitle = paste0('Number of PCs: ', horns$n, ', Number of Neighbours: ', round(perplexity), ', \nIterations: ', maxiter, ', Min Distance: 0.1\n'),
       colour='Cell Phase') +
  theme(
    legend.position = 'none',
    legend.text=element_text(size=9),
    plot.title = element_text(face='plain'),
    plot.subtitle =  element_text(size=9)
  ) +
  guides(shape=FALSE, alpha=FALSE) + 
  geom_vline(xintercept=0, linetype='dashed', colour='#777777', size=0.2) +
  geom_hline(yintercept=0, linetype='dashed', colour='#777777', size=0.2)
# SAVE


# TSNE & UMAP Plot Comparison
# Get Legend from Normalization FIX ----
nt <- theme(legend.position='none', plot.subtitle = element_blank())
title <- cowplot::ggdraw() + cowplot::draw_label('Visualising PCS by TSNE and UMAP')
p_grid <- cowplot::plot_grid(umap_1_CC + nt, tsne_1_CC + nt, umap_1_CC_split + nt, tsne_1_CC_split+nt, 
                             rel_heights = c(0.15,0.2))
comb_TSNE_UMAPplot <- cowplot::plot_grid(title, p_grid, legend_cc, ncol=1, rel_heights=c(0.1, 2, 0.1)) # SAVE 


saveRDS(sce_Gn, 'SCE_DIM_RED')
        
        
        


            
