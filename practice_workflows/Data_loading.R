library(SingleCellExperiment)
library(dplyr)

# Retrieve Grun pancreas (2016) data 
sce_G <- scRNAseq::GrunPancreasData(ensembl = TRUE)
rownames(sce_G) <- uniquifyFeatureNames(rownames(sce_G), rowData(sce_G)$symbol)

# query AH for ref genome to get Gene annotation TO USE!! ---- 
ah <- AnnotationHub::AnnotationHub()
edb <- ah[["AH83216"]]
annotations <- ensembldb::genes(edb, return.type="data.frame")  %>%
          dplyr::select(gene_id, gene_name, gene_biotype, description)




