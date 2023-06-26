library(seurat)
library(SingleCellExperiment)
library(scran)
library(scater)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(tibble)

# plot pearson residuels variance to identify all HVG
gene_meta <- as.data.frame(rowData(sce_Gn))

# rownames to column for plotting 
gene_meta <- rownames_to_column(gene_meta, 'Gene')


# data frame for labels 
hvg.10 <- gene_meta %>% 
  arrange(desc(sct.residual_variance)) %>% 
  slice_head(n=10)

# test <- gene_meta %>% filter(sct.variable==TRUE) : No Ribosomal Genes are HVG
# sum(grepl("^RP[SL]\\d+", rownames(test)))

# Get Threshold for HVG detection
threshold <- gene_meta %>% filter(sct.variable==FALSE) %>% top_n(sct.residual_variance, n=1) %>%
              dplyr::select(sct.residual_variance)

# PLOT of HVG ----
HVG_plot <- ggplot(gene_meta) + 
  geom_point(aes(x=sct.gmean, y=sct.residual_variance, 
                 color=dplyr::case_when(sct.residual_variance >= min(hvg.10$sct.residual_variance) ~ '#1b9e77',
                                        sct.variable ~ '#d95f02',
                                        TRUE ~ "#7570b3")), size=1.5, alpha=1/3) + 
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  theme_classic() + 
  labs(title='Top 10 Highly Variable Genes after SCTransform Normalization\n', 
       color='Highly Variable Genes',
       x='log(Mean Count)\n',y='log(Residual Variance)') +  
  guides(shape=FALSE, alpha=FALSE, colour=guide_legend(override.aes = list(size=3, alpha=1))) +
  theme(
    plot.title = element_text(face='plain', hjust=0.0),
    plot.caption = element_text(hjust=0.5),
    legend.position = 'none') 



# Add Labels and Line for 
HVG_plot <- HVG_plot + geom_text_repel(data=hvg.10, aes(sct.gmean, sct.residual_variance, label=Gene),
                          box.padding = 0.5,max.overlaps=Inf, size=3) + 
                          geom_hline(yintercept=threshold$sct.residual_variance, linetype='dashed', 
                                     colour='#777777', size=1)  
ggsave('Plots/Marker_gene_plots/HVG_plot.png', HVG_plot)

# save object
saveRDS(sce_Gn, 'HVG_Detection')
