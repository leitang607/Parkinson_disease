suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(dplyr))
library(ggplot2)
Convert('./stored each_subcluster h5 files/dans_umap.h5ad',dest = "h5seurat", overwrite = TRUE)
ad_dans <- LoadH5Seurat('./stored each_subcluster h5 files/dans_umap.h5seurat',assays="RNA")
for (i in unique(ad_dans$label_9)){
    ad_temp <- subset(ad_dans,label_9==i)
    Idents(ad_temp) <-'exp_condition1'
    all_markers <- FindAllMarkers(ad_temp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    write.csv(all_markers,file=paste0('./de_genes/',i,'.csv'))
}
