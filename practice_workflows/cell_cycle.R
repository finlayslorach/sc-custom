cyclin.genes <- grep("^CCN[ABDE][0-9]$", rowData(sce)$symbol)
cyclin.genes <- rownames(sce)[cyclin.genes]
