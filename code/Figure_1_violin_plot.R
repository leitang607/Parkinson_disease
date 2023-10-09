library(Seurat)
library(ggplot2)
ad <- readRDS("./sn_umap.Rdata")
Idents(ad) <- "label_subtype"
re_levels <- c('Astrocytes','Gluts','Inhibitory','GABA', # nolint
                'HdNs','DaNs','OPC','Oligo','Microglia','Immune','Endothelial','Vasc','Ependymal') # nolint
ad$label_subtype <- factor(ad$label_subtype,levels = re_levels)
Idents(ad) <- 'label_subtype' # nolint
violin_color <- c('#6c00bf','#00835A','#FFB307','#BF480D','#FF6600','#BF8219','#74A0FF','#0000ff','#D64849','#f252c5','#A19922','#807B30','#FFA388') # nolint
ftr_Astrocytes <- c('AQP4', 'SLC4A4')
ftr_Microglia <- c('MS4A7', 'CSF1R')
ftr_Immune <- c('CD44','CD48')
ftr_OPC <- c('VCAN','PDGFRA')
ftr_GABA <- c('SPHKAP','GRIK1')
ftr_Gluts <- c('SLC17A6','FAM19A1')
ftr_Inhibitory <- c('GAD1', 'GAD2')
ftr_HdNs <- c('HDC','HCRTR2')
ftr_DaNs <- c('SLC6A3','TH')
ftr_Endothelial <- c('FLT1','SLC2A1')
ftr_Vasc <- c('LAMA2','DCN')
ftr_Ependymal <- c('FOXJ1', 'SLC5A5')
ftr_Oligo <- c('MOBP','MOG')
violin_features <- c(ftr_Astrocytes,ftr_Gluts,ftr_Inhibitory,ftr_GABA,ftr_HdNs,ftr_DaNs,ftr_OPC,ftr_Oligo,ftr_Microglia,ftr_Immune,ftr_Endothelial,ftr_Vasc,ftr_Ependymal)# nolint
violin_figure <- VlnPlot(ad,features = violin_features,cols = violin_color,stack = TRUE,fill.by = "ident",flip=TRUE)+ scale_y_discrete(expand = c(0,0))
ggsave(violin_figure,width = 5,height = 5,filename = './Figures1/snc_violin_plot.pdf')# nolint

# For putamen
ad <- readRDS('./putamen_figout.Rdata')# nolint
Idents(ad) <- "label_subtype" 
re_levels <- c('Astrocytes','Gluts','eSPN','iSPN','dSPN','PVALB/TH IN','SST/NPY IN','CHAT IN','OPC','NFOL','Oligo','Microglia','Immune','Endothelial','Vasc') # nolint
ad$label_subtype <- factor(ad$label_subtype,levels = re_levels)
Idents(ad) <- 'label_subtype' # nolint
violin_color <- c('#6c00bf','#00835A','#1f77b4','#279e68','#98df8a','#c49c94','#ffbb78','#ff7f0e','#74A0FF','#0000ff','#D64849','#f252c5','#A19922','#807B30')# nolint
#plot violin PLot
ftr_Astrocytes <- c( 'AQP4','ALDOC')  # nolint
ftr_Microglia <- c('MS4A7', 'CSF1R')  # nolint
ftr_Immune <- c('CD44','CD48')  # nolint
ftr_OPC <- c('VCAN','PDGFRA')  # nolint
ftr_eSPN <- c('PPP1R1B','OTOF','CACNG5','TAC1')  # nolint
ftr_PVALB_TH_IN <-c('LHX6','IL1RAPL2') # nolint
ftr_iSPN <- c('DRD2', 'FIG4')   # nolint
ftr_dSPN <- c('MEIS1', 'RELN')    # nolint
ftr_SST_NPY_IN <- c('NPY','NOS1')  # nolint
ftr_Gluts <- c('SLC17A7', 'SATB2')  # nolint
ftr_Endothelial <- c('FLT1')  # nolint
ftr_Vasc <- c('DCN')  # nolint
ftr_Oligo <- c('MOBP','MOG') # nolint
ftr_CHAT_IN <- c('CHAT','NEB')  # nolint
violin_features <- c(ftr_Astrocytes,ftr_Gluts,ftr_eSPN,ftr_iSPN,ftr_dSPN,ftr_PVALB_TH_IN,ftr_SST_NPY_IN,ftr_CHAT_IN,ftr_OPC,ftr_Oligo,ftr_Microglia,ftr_Immune,ftr_Endothelial,ftr_Vasc) # nolint
putamen_violin <- VlnPlot(ad,features = violin_features,stack = TRUE,cols = violin_color,fill.by = "ident",flip=TRUE)+ scale_y_discrete(expand = c(0,0))
ggsave(putamen_violin,width = 5,height = 5,filename = './Figures1/putamen_violin.pdf')     # nolint                                                                                     
