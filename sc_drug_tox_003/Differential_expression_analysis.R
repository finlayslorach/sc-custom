################## LOAD DATA ############################
root <- '/hpc/scratch/hdd2/fs541623/scRNAseq/Run1_3_Cisplatin_080721/Secondary_analysis'
s <- readRDS(paste0(root, '/new_SCE_objects/FINAL_WORKFLOW/dim_reduction_obj.rds'))
sce_3 <- readRDS(paste0(root, '/new_SCE_objects/FINAL_WORKFLOW/raw_sce.rds'))
counts_corrected <- GetAssayData(s, assay='SCT', slot='counts')
theme_set(theme_bw())

## Get genome annotations ##
load_annotations <- function(){
  ah <- AnnotationHub::AnnotationHub()
  edb <- ah[["AH83216"]]
  annotations <- ensembldb::genes(edb, return.type="data.frame")  %>%
  dplyr::select(gene_id, gene_name, gene_biotype, description, symbol)
  return(annotations)
  }

############### DEG ACROSS TIME ########################
de_res <- sctransform::diff_mean_test(
  y=GetAssayData(s, assay='SCT', slot='counts'), group_labels=s$sample, only_pos = F, mean_th=0.25)
de_res <- readRDS(paste0(root, '/new_SCE_objects/FINAL_WORKFLOW/deg_cellthr2.rds'))

### get annotations ###
de_res.anno <- annotations %>% 
  dplyr::filter(gene_name %in% de_res$gene) %>% 
  dplyr::rename(gene='gene_name') %>%
  left_join(., de_res, by='gene') %>% 
  arrange(group1, pval_adj) 
write.csv(de_res.anno, paste0(root, '/new_tables/FINAL_WORKFLOW/deg_cellsthr2annotated.csv'))

#### DEG across controls ####
times <- c('0h','2h','8h','16h','24h','48h','72h')
DEG.between_contr <- 
  lapply(times[2:7], function(x) sctransform::diff_mean_test(
  y=counts_corrected, group_labels=s$sample, only_pos = F, 
  compare=c(x, '0h'), mean_th=0.25))
de_res_contr <- do.call(rbind, DEG.between_contr)
de_res_contr$group1 <- factor(de_res_contr$group1)
saveRDS(de_res_contr, paste0(root, '/new_SCE_objects/FINAL_WORKFLOW/deg_cellthr2_control.rds'))

### get annotations ###
annotations <- load_annotations()
de_res.anno <- annotations %>% 
  dplyr::filter(gene_name %in% de_res_contr$gene) %>% 
  dplyr::rename(gene='gene_name') %>%
  left_join(., de_res_contr, by='gene') %>% 
  arrange(group1, pval_adj) 
write.csv(de_res.anno, paste0(root, '/new_tables/FINAL_WORKFLOW/deg_cellsthr2annotated_control.csv'))

### get gene ids ###
de_res_anno$gene.ids <- mapIds(org.Hs.eg.db, keys =de_res_anno$gene, column='ENSEMBL', keytype='SYMBOL')
de_res_anno <- de_res_anno %>% dplyr::select(gene, group1, pval_adj,log2FC, gene.ids)
de_res_anno <- de_res_anno %>% split.data.frame(factor(.$group1)) %>% purrr::reduce(merge, by='gene')
write.csv(de_res_anno, paste0(root, '/new_tables/FINAL_WORKFLOW/deg_cellsthr2_control_IPAformatted.csv'))


####### use vs time ######
### get gene ids and format for IPA  ###
de_res$gene.ids <- mapIds(org.Hs.eg.db, keys =de_res$gene, column='ENSEMBL', keytype='SYMBOL')
de_res <- de_res %>% dplyr::select(gene, group1, pval_adj,log2FC, gene.ids)
de_path_rest <- de_res %>% split.data.frame(factor(.$group1)) %>% purrr::reduce(merge, by='gene')
write.csv(de_path_rest, paste0(root, '/new_tables/FINAL_WORKFLOW/deg_cellsthr2_IPA.csv'))

### VOLCANO PLOT ###
de_res <- readRDS(paste0(root, '/new_SCE_objects/FINAL_WORKFLOW/deg_cellthr2.rds'))
de_res_filt <- de_res %>% dplyr::filter(group1 %in% c('2h','8h','16h','24h','48h','72h'))
markers_vol <- group_by(de_res, group1) %>% 
  dplyr::filter(rank(pval_adj, ties.method='first') <= 6) %>%  
  dplyr::select(group1, gene, mean1, mean2, log2FC, zscore, emp_pval_adj, pval_adj)  %>% 
  dplyr::filter(group1 %in% c('2h','8h','16h','24h','48h','72h'))
### Plot most significant genes ###
ggplot(de_res_filt, aes(x=pmin(log2FC,10), y=pmin(-log10(pval_adj), 10))) + 
  geom_point(aes(color = pval_adj < 0.05))  + 
  geom_point(data=markers_vol, color='deeppink') +
  geom_text_repel(data = markers_vol, mapping = aes(label = gene), max.overlaps=Inf, force=3) + 
  theme(legend.position = "bottom") + 
  facet_wrap(~group1, ncol = 3) + ylab("-log10(pval_adj)") + 
  geom_vline(xintercept=0, linetype='dashed', colour='#777777', size=0.2) +
  geom_hline(yintercept=0, linetype='dashed', colour='#777777', size=0.2) +
  xlab("log2 mean fold-change")
ggsave(paste0(root, '/new_visualisations/FINAL_WORKFLOW/volcano_top6_log10pvaladj.png'), height=10, width=10)

############ HEATMAP ############################
marker_genes <- group_by(de_res, group1) %>% dplyr::filter(rank(pval_adj, ties.method = "first") <= 15)  
s <- GetResidual(s, features=marker_genes$gene, verbose=F)
sig.genes <- de_res %>% dplyr::filter(pval_adj < 0.05) %>% dplyr::filter(group1!='0h')
## plot most upregulated genes ##
s.s <- subset(s, subset=sample!='0h')
s.s <- GetResidual(s.s, features=sig.genes$gene, verbose=F)
DoHeatmap(s.s, group.by='sample', features=sig.genes$gene, slot = "scale.data") + 
scale_fill_gradient2(low = "blue", mid = "white", high = "#FF2800")
ggsave(paste0(root, '/new_visualisations/FINAL_WORKFLOW/heatmap_upreg_pval_all.png'), width=12, height=15)

## complex heatmap ##
s$pseudotime <- sce$slingPseudotime_1
annotations <- FetchData(s, vars=c('sample', 'pseudotime', 'mito.percent'))
colnames(annotations) <- c('Time', 'psuedotime', 'mitochondrial.percent')
ha.anno <- HeatmapAnnotation(df=annotations, 
    col = list(Time=c('0h'=cc[1], '2h'=cc[2],'8h'=cc[3],'16h'=cc[4],'24h'=cc[5],'48h'=cc[6],'72h'=cc[7]),
    psuedotime=colorRamp2(c(0,0.75, 1.5), c("green", "white", "blue")),
    mitochondrial.percent=colorRamp2(c(0, 15, 30), c("white", "red", "darkred"))
    ))
markers <- read.csv(paste0(root, '/new_tables/FINAL_WORKFLOW/Copy of Copy of deg_cellsthr2annotated_170821.csv'))
markers <- markers[!(markers %in% c(dna.gene, plat.genes))][1:20]
mark.toplot <- markers[markers$LABEL == '*',] %>% pull(symbol)
plat.genes <- markers %>% dplyr::filter(LABEL == 'Platinum') %>% pull(symbol)
dna.gene <- markers %>% dplyr::filter(LABEL == 'DNA') %>% pull(symbol)
counts.mat <- counts[sig.genes$gene,]
markers <- c(plat.genes, dna.gene, unique(mark.toplot))
plot.new()
png(paste0(root, '/new_visualisations/FINAL_WORKFLOW/DE_complex_heatmap1.png'), height=680, width=480)
Heatmap(counts.mat,  
    column_split = factor(s$sample),
    column_title = 'Differentially expressed genes across time',
    name='Expression',
    col = colorRamp2(c(-2, 0, 2), c("blue", "white", "#FF2800")),
    cluster_columns = F,
    show_column_dend = FALSE,
    cluster_column_slices = TRUE,
    column_title_gp = gpar(fontsize = 12),
    column_gap = unit(0.5, "mm"),
    width=unit(8,'cm'),
    height=unit(15, 'cm'),
    cluster_rows = F,
    show_row_dend = T,
    column_title_rot = 0,
    top_annotation = ha.anno,
    show_column_names = FALSE,
    show_row_names = F) +
rowAnnotation(link=anno_mark(at = which(rownames(counts.mat) %in% markers),
    labels = (rownames(counts.mat)[rownames(counts.mat) %in% markers]), 
    labels_gp = gpar(fontsize = 9, col='black'), padding = unit(0.5, "mm")))


##### SAVE OBJECT #####
saveRDS(s, paste0(root, '/new_SCE_objects/FINAL_WORKFLOW/DE_obj.rds'))
