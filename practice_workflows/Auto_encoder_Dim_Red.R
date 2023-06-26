setwd('C:/Users/fs541623/OneDrive - GSK/Desktop/scRNASeq/GrunPancreas_analysis_120121')
saveRDS(sce_Gn, 'SingleCellExpressionObject.RData')

library(SingleCellExperiment)
counts_1 <- as.data.frame(logcounts(sce_Gn)) %>% 
          rbind(Cell_Type = sce_Gn$Cluster_id)

write.csv(t(counts_1), 'SingleCellExpressionCounts.csv', col.names=FALSE)


