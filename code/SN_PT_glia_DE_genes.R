suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(dplyr))
library(ggplot2)
ad <- readRDS('./SN_umap.Rdata')
Idents(ad) <- 'label_5'
ad <- subset(ad,downsample=12000)
ad <- NormalizeData(ad, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(ad)
ad <- ScaleData(ad, features = all.genes)
for (i in unique(ad$label_5)){
    ad_temp = subset(ad,label_5==i)
    Idents(ad_temp) <- 'exp_condition1'
    print(paste0('current cell cluster is ',i))
    print(table(ad_temp$exp_condition1))
    print('------------')
    de_gene <- FindMarkers(ad_temp, ident.1 = "PD", ident.2 = "control",logfc.threshold = 0.25, test.use = "wilcox")
    write.csv(de_gene,file=paste0('./de_genes/snc_',i,'_degene.csv'))
}

ad_put <- readRDS('./PT_umap.Rdata')
ad_put <- NormalizeData(ad_put, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(ad_put)
ad_put <- ScaleData(ad_put, features = all.genes)
for (i in c('PVALB/TH IN','SST/NPY IN')){
    ad_temp = subset(ad_put,label_5==i)
    Idents(ad_temp) <- 'exp_condition1'
    print(paste0('current cell cluster is ',i))
    print(table(ad_temp$exp_condition1))
    print('------------')
    de_gene <- FindMarkers(ad_temp, ident.1 = "PD", ident.2 = "control",logfc.threshold = 0.25, test.use = "wilcox")
    write.csv(de_gene,file=paste0('./de_genes/putamen_',gsub('/','-',i),'_degene.csv'))
}
for (i in c('Astrocytes','Microglia','OPC','eSPN','dSPN','iSPN','Immune','CHAT IN','Endothelial','Vasc','Excitatory','MOL')){
    ad_temp = subset(ad_put,label_5==i)
    Idents(ad_temp) <- 'exp_condition1'
    print(paste0('current cell cluster is ',i))
    print(table(ad_temp$exp_condition1))
    print('------------')
    de_gene <- FindMarkers(ad_temp, ident.1 = "PD", ident.2 = "control",logfc.threshold = 0.25, test.use = "wilcox")
    write.csv(de_gene,file=paste0('./de_genes/putamen_',gsub('/','-',i),'_degene.csv'))
}

