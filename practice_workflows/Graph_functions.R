# TSNE Function
tsne_func <- function(df, colour){
  plot <- ggplot(df) + geom_point(aes_string(x='V1', y='V2', color=colour), size=0.8, shape=16, alpha=0.5) + 
    theme_classic() +
    labs(title='TSNE', 
         x='TSNE_1\n',
         y='TSNE_2\n',
         colour=colour) +
    coord_fixed(1) +
    theme(
      legend.title =element_text(size=10),
      legend.position=c(0.9, 0.2),
      legend.text=element_text(size=9)
    ) +
    guides(shape=FALSE, alpha=FALSE) + 
    geom_vline(xintercept=0, linetype='dashed', colour='#777777', size=0.2) +
    geom_hline(yintercept=0, linetype='dashed', colour='#777777', size=0.2)
  return(plot)
}


# Walds test comparison
get_dds_resultsD17vsD7 <- function(x, contrast, ref, n, A, B){
  # Subset Beta Cells Metadata
  clusters <- levels(metadata$cluster_ID)
  cluster_metadata <- metadata[which(metadata$cluster_ID == clusters[x]), ]
  rownames(cluster_metadata) <- cluster_metadata$Donor_Sample_id
  
  # Subset Gene Counts to Beta Cells
  counts <- as.matrix(agg_sce[[clusters[x]]])
  
  # Get Cluster Counts for Cells in Sample
  cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
  
  
  dds <- DESeqDataSetFromMatrix(cluster_counts, colData = cluster_metadata,
                                design= ~ Donor)
  
  # Distances between Rows for Data Vis
  rld <- rlog(dds, blind=TRUE)
  SampleDists <-  as.matrix(dist(t(assay(rld))))
  mds <- data.frame(cmdscale(SampleDists))
  mds <- cbind(mds, colData(rld))
  ggplot(mds) + geom_point(aes(X1,X2,color=Donor), size=4) + 
    labs(title=paste0('MDS Plot: ', clusters[5],'Cells'),
                                                                  xlab='MDS_1',
                                                                  ylab='MDS_2')
  ggsave(paste0("results/", clusters[5], "_specific_MDSplot.png"))
  
  
  # Hierachal Clustering 
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  # Heatmap 
  pheatmap::pheatmap(rld_cor, annotation=cluster_metadata[, c('Donor'), drop=F])
  # Differential Expression Analaysis 
  dds$Donor <- relevel(dds$Donor, ref = 'D7')
  dds <- DESeq(dds)
  
  # Dispersion Estimates 
  plotDispEsts(dds)
  
  # Results table for D3 vs D10
  resultsNames(dds)
  res <- results(dds, name='Donor_D17_vs_D7')
  res <- lfcShrink(dds, coef = contrast, type = 'apeglm', res = res)
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble()
  
  # Filter Table based on Significance Threshold 
  padj_cutoff <- 0.05
  sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
    dplyr::arrange(padj)
  
  # Get Top 20 Genes 
  norm_counts <- counts(dds, normalized=TRUE)
  top20_sig_genes <- sig_res %>%
    dplyr::arrange(padj) %>% 
    dplyr::pull(gene) %>%
    head(n=10)
  top20_sig_norm <- data.frame(normalized_counts) %>% 
    rownames_to_column(var='gene') %>% 
    dplyr::filter(gene %in% top20_sig_genes)
  
  # Nice Data frame of Sample_id and Counts for Plotting
  gathered_top20_sig <- top20_sig_norm %>%
    gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value =
             "normalized_counts")
  
  # Join with Metadata 
  agg <- agg %>% 
    mutate(Donor_Sample_id=stringr::str_replace_all(Donor_Sample_id, ' |,', '.'))
  gathered_top20_sig <- inner_join(agg[, c('Donor', "Sample", "Donor_Sample_id" )], 
                                   gathered_top20_sig, by = 
                                     c("Donor_Sample_id"                                                      
                                       ="samplename"))
  
 
  gathered_top20_sig %>% filter(Donor %in% c('D7', 'D17')) %>% 
    ggplot() +
    geom_point(aes(x = gene, 
                   y = normalized_counts+1, 
                   color = Donor), 
               position=position_jitter(w=0.1,h=0)) +
    scale_y_log10() +
    theme_classic() + 
    xlab("Genes") +
    ylab("log10 Normalized Counts") +
    ggtitle(paste0("Top 20 Significant DE Genes in ", clusters[5], "Cells")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0("results/", clusters[5], "top20.png"))
  
  
  res_table_thres <- res_tbl %>% 
    mutate(threshold = padj < 0.05) 
  labels_DE <- res_table_thres %>% 
    arrange(desc(log2FoldChange)) %>% 
    slice_head(n=5) 
  
  ## Volcano plot
  ggplot(res_table_thres) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
    ggtitle(paste0("Volcano plot of DE Genes in ",  clusters[5], "cells")) +
    xlab("log2 fold change") + 
    theme_classic() + 
    ylab("-log10 adjusted p-value") +
    scale_y_continuous(limits = c(0,50)) +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25))) + 
    ggrepel::geom_text_repel(data=labels_DE, aes(label=gene,log2FoldChange,-log10(padj)),
                             box.padding = 0.4,max.overlaps=Inf, size=3, nudge_x = 0.2, nudge_y = -0.1)
  ggsave(paste0("results/", clusters[5], "volcano.png"))
  
}


# Volcano Plot 
volcano_plot <- function(df, x, y, name){
  ggplot(res_table_thres) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  ggtitle(paste0("Volcano plot of DE Genes in ", name, "cells")) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 
}
  # Run the script on all clusters comparing stim condition relative to control condition
ggrepel::geom_text_repel(data=labels_DE, aes(label=gene,log2FoldChange,-log10(padj)),
                         box.padding = 0.4,max.overlaps=Inf, size=3, nudge_x = 0.2, nudge_y = -0.1)
