import pandas as pd
import scanpy as sc
import numpy as np
import seaborn as sns
from matplotlib import *
import matplotlib.pyplot as plt
rcParams['pdf.fonttype']=42
rcParams['ps.fonttype']=42
colors_gluts = plt.cm.Greens(np.linspace(0.2, 1, 6))
a = sns.color_palette(colors_gluts,6)
color_order = [colors.to_hex(i) for i in colors_gluts]

ad_glu = sc.read('../data/gluts_umap.h5')
order = ['Glu_MEIS2','Glu_MEIS1','Glu_CLMN','Glu_DPYD','Glu_PCSK5','Glu_LMX1A']
color_dict = dict(zip(order,color_order))
ad_glu.uns['label_subtype_colors'] = [color_dict.get(i)  for i in np.unique(ad_glu.obs['label_subtype'])]
fig,ax = plt.subplots(figsize=(5,4),dpi=600)
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.size']=7
plt.rcParams['font.weight']='normal'
sc.pl.umap(ad_glu,color='label_subtype',ax=ax,legend_loc='',frameon=False,title='',s=40)
fig.savefig('./figure/gluts_umap.png')

fig,ax = plt.subplots(figsize=(4,3),dpi=600)
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.size']=7
plt.rcParams['font.weight']='normal'
categories_order = ['Glu_MEIS2','Glu_MEIS1','Glu_CLMN','Glu_DPYD','Glu_PCSK5','Glu_LMX1A']
sc.pl.DotPlot(ad_glu,groupby='label_5',var_names=['MEIS2','MEIS1','CLMN','DPYD','PCSK5',\
                 'LMX1A','SORCS3','FOXP2'],categories_order=categories_order,use_raw=False,ax=ax).style(cmap='GnBu',largest_dot=250,dot_edge_lw=0).show()
fig.tight_layout()
fig.savefig('./Figures_pip_ver2/Figures2/exc_dotplot.pdf')



