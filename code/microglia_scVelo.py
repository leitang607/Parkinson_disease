import scanpy as sc
import scvelo as scv
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

vcl = scv.read("/data/dopamine_all.loom")
ad_mic = sc.read('./Figures_pipline/microglia_figout.h5')
ad = scv.utils.merge(ad_mic,vcl,id_length=16)
scv.pl.proportions(ad, groupby = 'label_5')
scv.pp.filter_and_normalize(ad, min_shared_counts=20)
scv.pp.moments(ad, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(ad,n_jobs=12) 
scv.tl.velocity(ad,mode='dynamical')
scv.tl.velocity_graph(ad)
scv.tl.terminal_states(ad)
scv.pl.velocity_embedding_stream(ad, basis='umap', color = 'label_5',legend_loc='right',alpha=0.65,s=11)
fig,ax = plt.subplots(figsize=(3,3),dpi=300)
scv.pl.velocity_embedding_stream(ad, basis='umap', color = 'label_5',ax=ax,legend_loc='right',alpha=0.65,s=11,save='Figures_pipline/Figures4/microglia_scvelo.pdf')
df = ad.var
df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]
kwargs = dict(xscale='log', fontsize=16)
with scv.GridSpec(ncols=3) as pl:
    pl.hist(df['fit_alpha'], xlabel='transcription rate', **kwargs)
    pl.hist(df['fit_beta'] * df['fit_scaling'], xlabel='splicing rate', xticks=[.1, .4, 1], **kwargs)
    pl.hist(df['fit_gamma'], xlabel='degradation rate', xticks=[.1, .4, 1], **kwargs)
scv.get_df(ad, 'fit*', dropna=True).head()
from matplotlib import *
#Define colour map
colors2 = plt.cm.Reds(np.linspace(0, 1, 100))
colors3 = plt.cm.Blues_r(np.linspace(0,0.9,100))
colorsComb = np.vstack([colors3, colors2])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)
scv.tl.latent_time(ad)
scv.pl.scatter(ad, color='latent_time', color_map=mymap,save='Figures_pipline/Figures4/microglia_latent_time.pdf')

top_genes = ad.var['fit_likelihood'].sort_values(ascending=False).index[:250]
scv.pl.heatmap(ad, var_names=top_genes, sortby='latent_time', col_color='label_5', n_convolve=200,save='Figures_pipline/Figures4/scvelo_heatmap_microglia.png')

# draw marker gene in umap
vmax = 5.5
vmin = 2
gene = 'VAV3'
fig,ax = plt.subplots(figsize=(2,2),dpi=300)
ax2 = fig.add_axes([0.7, 0.15, 0.2, 0.2])
norm = colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=mymap),aspect=8,ax=ax2)
cb.set_ticks([vmin,vmax])
ax2.axis('off')
sc.pl.umap(ad_mic,color=gene,cmap=mymap,vmax=vmax,vmin=vmin,s=8,frameon=False,colorbar_loc=None,ax=ax)
fig.savefig(f'./Figures_pipline/Figures4/microglia_{gene}_umap.pdf')

vmax = 5
vmin = 0
gene = 'GPNMB'
fig,ax = plt.subplots(figsize=(2,2),dpi=300)
ax2 = fig.add_axes([0.7, 0.15, 0.2, 0.2])
norm = colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=mymap),aspect=8,ax=ax2)
cb.set_ticks([vmin,vmax])
ax2.axis('off')
sc.pl.umap(ad_mic,color=gene,cmap=mymap,vmax=vmax,vmin=vmin,s=8,frameon=False,colorbar_loc=None,ax=ax)
fig.savefig(f'./Figures_pipline/Figures4/microglia_{gene}_umap.pdf')

vmax = 5
vmin = 1
gene = 'TNFAIP8L3'
fig,ax = plt.subplots(figsize=(2,2),dpi=300)
ax2 = fig.add_axes([0.7, 0.15, 0.2, 0.2])
norm = colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=mymap),aspect=8,ax=ax2)
cb.set_ticks([vmin,vmax])
ax2.axis('off')
sc.pl.umap(ad_mic,color=gene,cmap=mymap,vmax=vmax,vmin=vmin,s=8,frameon=False,colorbar_loc=None,ax=ax)
fig.savefig(f'./Figures_pipline/Figures4/microglia_{gene}_umap.pdf')

vmax = 3
vmin = 0
gene = 'MKI67'
fig,ax = plt.subplots(figsize=(2,2),dpi=300)
ax2 = fig.add_axes([0.7, 0.15, 0.2, 0.2])
norm = colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=mymap),aspect=8,ax=ax2)
cb.set_ticks([vmin,vmax])
ax2.axis('off')
sc.pl.umap(ad_mic,color=gene,cmap=mymap,vmax=vmax,vmin=vmin,s=8,frameon=False,colorbar_loc=None,ax=ax)
fig.savefig(f'./Figures_pipline/Figures4/microglia_{gene}_umap.pdf')
