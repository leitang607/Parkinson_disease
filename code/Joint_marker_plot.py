import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import *
rcParams['pdf.fonttype']=42
rcParams['ps.fonttype']=42

ad_hu = sc.read('../data/homo_DaNs_umap.h5')
ad_mk = sc.read('../data/DaNs_umap.h5')
marker_dict = dict()
marker_dict['SOX6+_RYR3']=['RYR3','SOX6','PEX5L','PLD5','FMN1']
marker_dict['SOX6+_CUX2']=['CUX2','OLFM3','DPP10','PSD3','COL11A1']
marker_dict['SOX6-_EBF1']=['EBF1','CDH13','LRRTM4','MAML3','UNC5D']
marker_dict['SOX6-_CALCR']=['RBFOX1','CALCR','IL1RAPL2','CDH6','NTNG1']
marker_dict['SOX6-_CBLB']=['CBLB','NEGR1','KCNAB1','GRIA4','GRID2','PTPRD']
marker_dict['SOX6-_COBLL1']=['COBLL1','TMEM132D','PCDH9','CDH18','TMTC1','GRIA4']
marker_dict['SOX6-_SORCS3']=['SORCS3','RELN','ANKFN1','CACNA2D1','KCNH1']
final_order = ['SOX6-_SORCS3', 'SOX6-_CBLB', 'SOX6-_COBLL1', 'SOX6-_CALCR', 'SOX6-_EBF1', 'SOX6+_RYR3', 'SOX6+_CUX2']
var_names=[]
for i in final_order:
    var_names = var_names+marker_dict.get(i)
fig,ax = plt.subplots(figsize=(10,3),dpi=600)
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.size']=7
plt.rcParams['font.weight']='normal'
sc.pl.DotPlot(ad_mk,groupby='label_9',var_names=var_names,categories_order=final_order,use_raw=False,ax=ax).style(cmap='GnBu',dot_edge_lw=0,largest_dot=160).show()
fig.tight_layout()
fig.savefig('./Figures/joint_marker_mk_DaNs_dotplot.pdf')

hu_order = ['SOX6_AGTR1','SOX6_PART1','SOX6_DDT','CALB1_CRYM_CCDC68','CALB1_PPP1R17','CALB1_CALCR','CALB1_RBP4','SOX6_GFRA2','CALB1_TRHR','CALB1_GEM']
hu_order.reverse()
fig,ax = plt.subplots(figsize=(10,3),dpi=600)
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.size']=7
plt.rcParams['font.weight']='normal'
sc.pl.DotPlot(ad_hu,groupby='label_subclass',var_names=var_names,categories_order=hu_order,use_raw=False,ax=ax).style(cmap='GnBu',dot_edge_lw=0,largest_dot=160).show()
fig.tight_layout()
fig.savefig('./Figures/joint_marker_hu_DaNs_dotplot.pdf')
