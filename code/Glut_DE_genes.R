suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(dplyr))
library(ggplot2)
Convert('../data/gluts_umap.h5ad',dest = "h5seurat", overwrite = TRUE)
ad_exc <- LoadH5Seurat('../data/gluts_umap../data/gluts_umap.h5seurat',assays="RNA")
vulnerable_sort <- c()
for (i in ad_exc$label_5){
    if (i %in% c('Ex_MEIS2','Ex_MEIS1')){
        vulnerable_sort <- c(vulnerable_sort,'vulnerable')
    }else {
        if (i %in% c('Ex_LMX1A')){
            vulnerable_sort <- c(vulnerable_sort,'resistant')
        }else{vulnerable_sort <- c(vulnerable_sort,'other')}
    }
}
ad_exc$vulnerable_cluster <- vulnerable_sort
ad_exc <- subset(ad_exc,vulnerable_cluster %in% c('vulnerable','resistant'))
Idents(ad_exc)<-'vulnerable_cluster'
all_markers <- FindAllMarkers(ad_exc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all_markers,file='./de_genes/exc_vulnerable_resistant_cluster_markers.csv')
