import scanpy as sc
import pandas as pd
import numpy as np
from matplotlib import *
import matplotlib.pyplot as plt
#Define a nice colour map for gene expression
colors2 = plt.cm.Reds(np.linspace(0, 0.7, 100))
colors3 = plt.cm.Blues_r(np.linspace(0,1,100))
colorsComb = np.vstack([colors3, colors2])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)


ad_exc = sc.read('./exc_figout.h5')
def plot_umap_genes(umap_position,exp_matrix,sample_area,gene,vmax_para,vmin_para):
    fig,ax = plt.subplots(figsize=(2.1,2),dpi=600)
    ax2 = fig.add_axes([0.75, 0.1, 0.15, 0.2])
    ax.scatter(umap_position['X'],umap_position['Y'],c=exp_matrix[gene],cmap=mymap,linewidths=0,vmax=vmax_para,vmin=vmin_para,s=0.2)
    norm = colors.Normalize(vmin=vmin_para, vmax=vmax_para, clip=True)
    cb = fig.colorbar(cm.ScalarMappable(norm=norm,cmap=mymap),aspect=8,ax=ax2)
    cb.set_ticks([vmin_para,vmax_para])
    cb.ax.tick_params(direction = 'out',width =0.4)
    cb.outline.set_color('black')
    cb.outline.set_linewidth(0.5)
    ax.axis('off')
    ax2.axis('off')
    fig.tight_layout()
    fig.savefig(f'./Figures1/{sample_area}_{gene}.png')
