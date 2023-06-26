library(Seurat)
library(plyr)
library(scran)
library(scater)
library(plyr)
library(dplyr)
library(ggplot2)
sceG

# check expression distribution of FILTERED cells 
hist_exp <- hist(log10(sceG$Total_Reads), xlab='# of reads before normalization', breaks=100)

# create temp copy for testing LIB NORM ---- 
sce_Gn <- sceG
sce_Gn <- computeLibraryFactors(sce_Gn)
sce_Gn <- logNormCounts(sce_Gn)

# Check Donor variation 
sce_Gn <- runPCA(sce_Gn, name='quickPCA', BSPARAM=BiocSingular::IrlbaParam())
PCA_quick <- plotReducedDim(sce_Gn, dimred="QUICKPCA", colour_by="Donor") +
  geom_hline(yintercept=0, linetype='dashed', colour='#777777') +
  geom_vline(xintercept=0, linetype='dashed', colour='#777777') +
  labs(subtitle='Donor Variation', 
       x='PC1 (26%)',
       Y='PC2 (12%)',
       caption = 'Library size Normalization') +  
  coord_fixed(1) +
  theme(
    plot.caption = element_text(hjust=0.0),
    legend.title = element_blank(),
    legend.justification = c('right', 'top'),
    legend.box.background = element_rect(color='black', size=0.1),
    legend.box.just = 'right',
    legend.margin = margin(6, 6, 6, 6),
    legend.text = element_text(size = 8)) 
PCA_quick <- PCA_quick +  guides(color=guide_legend(override.aes = list(size=3, alpha=1)))

# check whether to regress on cell cycle variation 
# check all genes w/o QC
cyclin.genes <- grep("^CCN[ABDE][0-9]$", rowData(sce_Gn)$symbol)
cyclin.genes <- rownames(sce_Gn)[cyclin.genes]

# pre-trained marker gene pairs to classify phase boundaries
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))

# ensembl names to compare to hs.pairs
cell_phase <- annotations %>% 
      dplyr::select(gene_id, gene_name) %>% 
      dplyr::rename(ensembl=gene_id, symbol=gene_name) %>%
      join(as.data.frame(rowData(sce_Gn)), by='symbol') 

# remove NA & duplicates - TO DO::REMOVE ----
temp_cc <- sce_Gn
cell_phase <- cell_phase[!is.na(cell_phase$originalName) & !duplicated(cell_phase$originalName),] 
rowData(temp_cc) <- join(as.data.frame(rowData(temp_cc)), cell_phase, by='symbol')

# CC phases using traning set ----
assignments <- cyclone(temp_cc, hs.pairs, gene.names=rowData(temp_cc)$ensembl, verbose=FALSE)


# Run PCA on cyclin genes TO DO MAKE PLOT NICER ----
PCA_CC_check <- plotReducedDim(sce_Gn, dimred="QUICKPCA", colour_by=I(assignments$phases)) +
  geom_hline(yintercept=0, linetype='dashed', colour='#777777') +
  geom_vline(xintercept=0, linetype='dashed', colour='#777777') +
  labs(title='Gene variation from Cell Cycle', 
       x='PC1',
       y='PC2',
       caption = 'Library Size Normalization of Counts') +  
  coord_fixed(1) +
  theme(
    plot.title = element_text(face='plain', hjust=0.5),
    plot.caption = element_text(hjust=0.0),
    legend.position = c(.5,.9),
    legend.box = 'horizontal',
    legend.direction = 'horizontal',
    legend.title = element_blank(),
    legend.justification = c('right', 'bottom'),
    legend.box.just = 'right',
    legend.margin = margin(6, 6, 6, 6),
    legend.text = element_text(size = 8))

# large differences in CC variation so regress out
# Seurat Normalization regularized binomial distribution
seurat_G <- as.Seurat(sce_Gn)
seurat_G <- SCTransform(seurat_G, new.assay.name = 'SCTransform')
# convert back to SCE 
sce_Gn <- as.SingleCellExperiment(seurat_G, assay='SCTransform')


# check read count distirbution after normalization 
total_reads <- colSums(counts(sce_Gn))
hist_exp1 <- hist(log10(total_reads), xlab='log[# of total reads after normalization]', breaks=150)


# Compare to Deconv Normalization ----
# Normalize Endogenous Genes
set.seed(0101)
quickclusters <- quickCluster(sceG, min.size=100)
sce_deconv <- computeSumFactors(sceG, cluster=quickclusters)
summary(sizeFactors(sce_deconv))

# Normalize Spike-Ins
sce_deconv <- computeSpikeFactors(sceG, spikes ='ERCC')
summary(sizeFactors(sce_deconv, 'ERCC'))

# Log Transformation
sce_deconv <- logNormCounts(sce_deconv, use.altexps='ERCC')

# Cell Phases Accounting for Unknown Classification
rownames(assignments$scores) <- colnames(sce_deconv)

phase <- rep("S", ncol(sce_deconv))
phase[assignments$scores$G1 > 0.5] <- "G1"
phase[assignments$scores$G2M > 0.5] <- "G2M"
phase[assignments$scores$G1 > 0.5 & assignments$scores$G2M > 0.5] <- "Unknown"

# Add Phase Column 
assignments$scores$phase <- phase


# PCA Plots of Diff Normalization Strategies FIX ----
# temp_cc == LS Normalization
norm_methods <- list('Deconvolution'=sce_deconv, 'SCTransform'=sce_Gn, 
                     'Library Size'=temp_cc)

output <- list()

j <- 1 # NEED TO FIX?? ----
for(i in norm_methods){
  pca <- runPCA(i, name='PCA',  BSPARAM=BiocSingular::IrlbaParam())
  var_1 <- attr(reducedDim(pca, 'PCA'), "percentVar" )[1]
  var_2 <- attr(reducedDim(pca, 'PCA'), "percentVar" )[2]
  pca <- reducedDim(pca, 'PCA')[,1:2]
  df <- data.frame(pca, phases=assignments$scores)
  output[[j]] <- ggplot(df) + 
    geom_point(aes(x=PC1,y=PC2, color=phases.phase), size=1, alpha=1/2) + 
    geom_hline(yintercept=0, linetype='dashed', colour='#777777') +
    geom_vline(xintercept=0, linetype='dashed', colour='#777777') +
    theme_classic() + 
    coord_fixed(1) +
    labs(x=paste0('PC1 (', round(var_1), '%)'), 
         y=paste0('PC2 (', round(var_2), '%)'), color="Cell Phases") + 
    theme(
      plot.title = element_text(face='plain', hjust=0.5),
      plot.caption = element_text(hjust=0.5),
      legend.position = c(.7, .2),
      legend.box = 'horizontal',
      legend.justification = c('right', 'bottom'),
      legend.box.just = 'right',
      legend.margin = margin(6, 6, 6, 6),
      legend.text = element_text(size = 8)) +
      guides(shape=FALSE, alpha=FALSE) 
  print(output[[j]])
  j <- j+1
}


# Grid for 1x Donor Plot and 3x Norm PCA
norm_tit <- cowplot::ggdraw() + 
  cowplot::draw_label('PCA: Donor and Cell Phase Variation ')
norm_don_gr <- cowplot::plot_grid(PCA_quick, output[[1]] + 
                                    labs(subtitle='Method: Deconvolution') +
                                theme(legend.position = 'none') ,
                              output[[2]]  + 
                                labs(subtitle='Method: SCTransform') +
                                theme(legend.position = 'none'), 
                              output[[3]]  +
                                labs(subtitle='Method: Library Size') +
                                theme(legend.position='none'), rel_widths = c(2, 2, 3, 2))
# Combined Grid Plot 
com_norm_dongr <- cowplot::plot_grid(norm_tit, norm_don_gr + guides(color=guide_legend(override.aes = list(size=3, alpha=1))), 
                                  legend_cc, ncol=1, 
                                  rel_heights=c(0.1, 2, 0.1))  
ggsave('Plots/Norm_plots/norm_don.png', com_norm_dongr)
# PCA Norm Strategies : TO FIX ----
purrr::map(norm_methods, function(norm){
  norm <- runPCA(norm, name='PCA',  BSPARAM=BiocSingular::IrlbaParam())
  pca <- reducedDim(norm, 'PCA')[,1:2]
  df <- data.frame(pca, phases=assignments$scores)
    ggplot(df, aes(PC1, PC2, color=phases.phase), size=1, alpha=1/2) + 
    theme_classic() + 
    coord_fixed(1) +
    labs(x='PC1',
         y='PC2',
         colour='Cell Phase',
         subtitle=paste0('Method:')) + 
    theme(plot.title = element_text(face='plain', hjust=0.5),
          plot.caption = element_text(hjust=0.5),
          legend.position = 'none') +
          guides(shape=FALSE, alpha=FALSE) 
}) %>%
  cowplot::plot_grid(plotlist=.) 
norm_title <- cowplot::ggdraw() + cowplot::draw_label('Choosing Normalization Method')
norm_grid <- cowplot::plot_grid(norm_title, norm_PCA, legend_cc, ncol=1, 
                                       rel_heights =c(0.1,2, 0.1))



# PCA plot split by cell cycle phases 
quickpca <- reducedDim(sce_Gn, 'QUICKPCA')[,1:2]
df_split <- data.frame(quickpca, phases=assignments$scores)
CC_split <- ggplot(df_split) + 
  geom_point(aes(x=quickpca_1,y=quickpca_2, color=phases.phase, shape='.',alpha=0.1)) +
  facet_wrap(~phases.phase) +
  theme_classic() + 
  labs(title='PCA plot: Variation From Cell phase',
      color="Cell Phases",
      x="PC_1",y="PC_2") + theme(
        legend.box = 'horizontal',
        legend.justification = c('bottom'),
        legend.box.just = 'bottom',
        legend.title = element_blank(),
        legend.position='bottom'
      ) +
  guides(shape=FALSE, alpha=FALSE, color=guide_legend(override.aes = list(size=3, alpha=1)))
# Get Cell Cycle Legend
legend_cc <- get_legend(CC_split)

# save norm objects 
saveRDS(seurat_G, 'Seurat_Grun_Norm')
# save SCE normalized using seurat
saveRDS(sce_Gn, 'SCE_G_Norm')
saveRDS(sce_G, 'SCE_G_Deconv_Norm')




