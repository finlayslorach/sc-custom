######### LOAD PACKAGES ##########
load_packages<- function(){library(SingleCellExperiment)
  library(dplyr)
  library(scran)
  library(scater)
  library(Seurat)
  library(magrittr)
  library(BiocGenerics)
  library(AnnotationHub)
  library(EnsDb.Hsapiens.v86)
  library(DESeq2)
  library(edgeR)
  library(Seurat)
  library(forcats)
  library(ggplot2)
  library(reshape2)
  library(ggrepel)
  library(batchelor)
  library(stringr)
  library(RColorBrewer)
  library(SC3)
  library(dendextend)
  library(dynamicTreeCut)
  library(limma)
  library(SeuratWrappers)
  library(Nebulosa)
  library(cluster, quietly = TRUE)
  library(harmony)
  library(clustree)
  library(destiny)
  library(plotly)
  library(org.Hs.eg.db)
  library(edgeR)
  library(ComplexHeatmap)
  library(circlize)
  library(VennDiagram)
}

####'GET GENOME ANNOTATIONS ####
load_annotations <- function(){
  ah <- AnnotationHub::AnnotationHub()
  edb <- ah[["AH83216"]]
  annotations <- ensembldb::genes(edb, return.type="data.frame")  %>%
    dplyr::select(gene_id, gene_name, gene_biotype, description, symbol)
  return(annotations)
}

#####'LOAD DATA #####
X <- read.csv('/hpc/scratch/hdd2/fs541623/scRNAseq/Run1_3_Cisplatin_080721/Secondary_analysis/Copy of UMI_matrix_150821.csv')
root <- '/hpc/scratch/hdd2/fs541623/scRNAseq/Run1_3_Cisplatin_080721/Secondary_analysis'
load_packages()
theme_set(theme_bw())
cc <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"),brewer.pal(9,"Set3"))

#########' FORMAT DATA ##########
Y <- X %>% dplyr::select(-c(Gene, gene.id))
rownames(Y) <- X$gene.id
sce_3 <- SingleCellExperiment(list(counts=Y), rowData=DataFrame(symbol=X$Gene)) %>% 
  set_rownames(uniquifyFeatureNames(rownames(.), rowData(.)$symbol))
saveRDS(sce_3, paste0(root, '/new_SCE_objects/FINAL_WORKFLOW/raw_sce.rds'))
## get sample labels ##
sce_3$sample <- str_extract(colnames(sce_3), '\\d+h')
sce_3$sample[is.na(sce_3$sample)] <- '0h'
sce_3$sample <- factor(sce_3$sample, levels=c('0h', '2h', '8h','16h','24h', '48h','72h'))
## get day labels ##
days <- str_extract(colnames(Y), '.*_day\\d+')
days[is.na(days)] <- 'Other'
sce_3$day <- days

########### QC #############################
counts(sce_3) <- data.matrix(counts(sce_3))
sce_3 <- addPerCellQC(sce_3, 
  subsets=list(mito=grepl('^MT-', rowData(sce_3)$symbol), 
  ribo=grepl('^RP', rowData(sce_3)$symbol)))

#### QC plots ####
plotColData(sce_3, x="sample", y="sum") + scale_y_log10() 
ggsave(paste0(root, '/new_visualisations/FINAL_WORKFLOW/QC_sum_logscale.png'))
discard <- sce_3$sum < 600 | sce_3$subsets_mito_percent < 1 | sce_3$sum > 20000 | sce_3$detected < 500
sce_3$discard <- discard
filtered_sce_3 <-sce_3[,!(discard)]
filtered_sce_3$lowgenes <- filtered_sce_3$detected < 500
plotColData(filtered_sce_3, x="sample", y="sum", colour_by = 'lowgenes') + scale_y_log10() 
ggsave(paste0(root, '/new_visualisations/FINAL_WORKFLOW/QC_sum_logscale_afterfiltering.png'))


############## NORMALIZATION ################### 
set.seed(123)
s <- CreateSeuratObject(counts(filtered_sce_3)) %>%
  PercentageFeatureSet(.,pattern='MT-', col.name='mito.percent') %>%
  AddMetaData(.$mito.percent > 20, col.name='high.mito') %>% 
  SCTransform(., verbose=F, return.only.var.genes = F, 
  variable.features.n = 1000)%>%
  RunPCA(assay='SCT', features=VariableFeatures(.))

### Write HVG to file ### 
annotations <- load_annotations()
top_hvg <- HVFInfo(s) %>% 
  mutate(., bc = rownames(.)) %>% 
  arrange(desc(residual_variance)) %>% 
  top_n(1000, residual_variance)  %>% 
  rename(gene_name=bc) 
hvg.anno <- annotations %>% 
  dplyr::filter(gene_name %in% top_hvg$gene_name) %>% 
  left_join(., top_hvg, by='gene_name') %>% 
  arrange(desc(residual_variance))
write.csv(hvg.anno, paste0(root, '/new_tables/FINAL_WORKFLOW/hvglist.csv'))

### Write PCA loadings to file ###
PCAloadings1 <- Loadings(s) %>% 
  as.data.frame(.) %>% dplyr::select(PC_1, PC_2) %>% 
  tibble::rownames_to_column(var='gene_name') %>% 
  left_join(., annotations, by='gene_name') %>% 
  arrange(-PC_1)
PCAloadings2 <- Loadings(s) %>% 
  as.data.frame(.) %>% dplyr::select(PC_1, PC_2) %>% 
  tibble::rownames_to_column(var='gene_name') %>% 
  left_join(., annotations, by='gene_name') %>% 
  arrange(-PC_2)
write.csv(PCAloadings1, paste0(root, '/new_tables/FINAL_WORKFLOW/PC1_loading_ordered.csv'))
write.csv(PCAloadings2, paste0(root, '/new_tables/FINAL_WORKFLOW/PC2_loading_ordered.csv'))


############# DIM REDUCTION ###################
horns <- PCAtools::parallelPCA(GetAssayData(s, assay='SCT', slot='scale.data')[VariableFeatures(s),])
horns$n ## 23 
ElbowPlot(s) ## 10
ggsave(paste0(root, '/new_visualisations/FINAL_WORKFLOW/Elbow_plot.png'))
s$sample <- filtered_sce_3$sample
s$day <- filtered_sce_3$day

########### PCA #############
to.plot <- as.data.frame(Embeddings(s[['pca']])) %>% 
  mutate(sample=s$sample) 
ggplot(to.plot) + 
  geom_point(aes(PC_1, PC_2, color=sample), size=2, alpha=0.8)+
  coord_fixed(1) +
  geom_vline(xintercept=0, linetype='dashed', colour='#777777', size=0.2) +
  geom_hline(yintercept=0, linetype='dashed', colour='#777777', size=0.2) +
  scale_colour_manual(values=c(cc[1], cc[2], cc[3], cc[4], cc[5], cc[6], cc[7]), 
  guide = guide_legend(override.aes = list(size = 5), title='Timepoint'))
ggsave(paste0(root, '/new_visualisations/FINAL_WORKFLOW/PCA_plot.png'))
## Run PCA on subset ##
s.s <- subset(s, subset=day!='Other')
s.s <- s.s %>% RunPCA() 
to.plot <- as.data.frame(Embeddings(s.s[['pca']])) %>% 
  mutate(day=s.s$day) 
ggplot(to.plot) + 
  geom_point(aes(PC_1, PC_2, color=day), size=2)+
  coord_fixed(1) +
  geom_vline(xintercept=0, linetype='dashed', colour='#777777', size=0.2) +
  geom_hline(yintercept=0, linetype='dashed', colour='#777777', size=0.2) +
  scale_colour_manual(values=c(cc[1], cc[2], cc[3], cc[4], cc[5], cc[6], cc[7]), 
  guide = guide_legend(override.aes = list(size = 5), title='Timepoint'))
ggsave(paste0(root, '/new_visualisations/FINAL_WORKFLOW/PCA_plot_batcheffect.png'))

########## MDS ###############
d <- dist(t(GetAssayData(s, slot = "scale.data")))
mds <- cmdscale(d = d, k = 8)
colnames(mds) <- paste0("MDS_", 1:2)
s[["mds"]] <- CreateDimReducObject(embeddings = mds, key = "MDS_", assay = DefaultAssay(s))
to.plot <- as.data.frame(Embeddings(s[['mds']])) %>% 
  mutate(sample=s$sample) 
ggplot(to.plot) + 
  geom_point(aes(MDS_1, MDS_2, color=sample), size=2)+
  coord_fixed(1) +
  geom_vline(xintercept=0, linetype='dashed', colour='#777777', size=0.2) +
  geom_hline(yintercept=0, linetype='dashed', colour='#777777', size=0.2) +
  scale_colour_manual(values=c(cc[1], cc[2], cc[3], cc[4], cc[5], cc[6], cc[7]), 
  guide = guide_legend(override.aes = list(size = 5), title='Timepoint'))
ggsave(paste0(root, '/new_visualisations/FINAL_WORKFLOW/mds_plot.png'))

## Run MDS on subset ##
d <- dist(t(GetAssayData(s.s, slot = "scale.data")))
mds <- cmdscale(d = d, k = 8)
colnames(mds) <- paste0("MDS_", 1:2)
s.s[["mds"]] <- CreateDimReducObject(embeddings = mds, key = "MDS_", assay = DefaultAssay(s.s))
to.plot <- as.data.frame(Embeddings(s.s[['mds']])) %>% 
  mutate(day=s.s$day) 
ggplot(to.plot) + 
  geom_point(aes(MDS_1, MDS_2, color=day), size=2)+
  coord_fixed(1) +
  geom_vline(xintercept=0, linetype='dashed', colour='#777777', size=0.2) +
  geom_hline(yintercept=0, linetype='dashed', colour='#777777', size=0.2) +
  scale_colour_manual(values=c(cc[1], cc[2], cc[3], cc[4], cc[5], cc[6], cc[7]), 
  guide = guide_legend(override.aes = list(size = 5), title='Timepoint'))
ggsave(paste0(root, '/new_visualisations/FINAL_WORKFLOW/mds_plot_batcheffect.png'))

### 3D plots ###
df <- as.data.frame(Embeddings(s, reduction='mds')) %>% mutate(sample=s$sample)
plot_ly(data = df, x=df$MDS_1, y=df$MDS_2, z=df$MDS_3, type="scatter3d", mode="markers", color=df$sample, size = 2)

### SAVE OBJECT ###
saveRDS(s, paste0(root, '/new_SCE_objects/FINAL_WORKFLOW/dim_reduction_obj.rds'))
