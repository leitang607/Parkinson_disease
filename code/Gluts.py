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
def milo_cell_rate(adata,k=25,d=25,prop=0.1,rep='X_aucell'):
    sc.pp.neighbors(adata, use_rep=rep, n_neighbors=k, n_pcs=d,key_added='milo')
    milo.make_nhoods(adata, neighbors_key='milo', prop=prop)
    nhood_size = adata.obsm['nhoods'].toarray().sum(0)
    plt.hist(nhood_size, bins=20);
    plt.xlabel('# cells in neighbourhood');
    plt.ylabel('# neighbouthoods');
    milo.count_nhoods(adata, sample_col="sampleid")
    adata.uns['nhood_adata']
    mean_n_cells = adata.uns['nhood_adata'].X.toarray().mean(1)
    plt.plot(nhood_size, mean_n_cells, '.');
    plt.xlabel('# cells in nhood');
    plt.ylabel('Mean # cells per sample in nhood')
    milo.DA_nhoods(adata, design="~exp_condition1", model_contrasts='exp_condition1PD-exp_condition1control')
    milo_results_salmonella = adata.uns["nhood_adata"].obs.copy()
    def plot_milo_diagnostics(adata):
        nhood_adata = adata.uns['nhood_adata'].copy()
        alpha = 0.1 
    
        with matplotlib.rc_context({'figure.figsize':[12,12]}):
            plt.subplot(2,2,1);
            plt.hist(nhood_adata.obs['PValue'], bins=20);
            plt.xlabel('Uncorrected P-value');
            plt.subplot(2,2,2);
            plt.scatter(nhood_adata.obs['PValue'], nhood_adata.obs['SpatialFDR'], s=3);
            plt.xlabel('Uncorrected P-value');
            plt.ylabel('SpatialFDR');
            plt.subplot(2,2,3);
            plt.scatter(nhood_adata.obs['logFC'], -np.log10(nhood_adata.obs['SpatialFDR']), s=3);
            plt.axhline(y=-np.log10(alpha), color='red', linewidth=1, label=f'{int(alpha*100)} % SpatialFDR');
            plt.legend();
            plt.xlabel('log-Fold Change');
            plt.ylabel('- log10(SpatialFDR)');
            plt.tight_layout()
            df = nhood_adata.obs
            emp_null = df[df['SpatialFDR'] >= alpha]['logFC'].mean()
            df['Sig'] = df['SpatialFDR'] < alpha
    
            plt.subplot(2,2,4);
            sns.scatterplot(data=df, x="logCPM", y="logFC", hue='Sig')
            plt.axhline(y=0, color='grey', linewidth=1)
            plt.axhline(y=emp_null, color='purple', linewidth=1);
            plt.xlabel('Mean log-counts');
            plt.ylabel('log-Fold Change');
            plt.show()
    plot_milo_diagnostics(adata)
    milopy.utils.build_nhood_graph(adata)
    with matplotlib.rc_context({'figure.figsize':[5,5]}):
        milopy.plot.plot_nhood_graph(adata, alpha=0.1, min_size=5, plot_edges=False)
        sc.pl.umap(adata, color='label_5', legend_loc='on data')
    milopy.utils.annotate_nhoods(adata, anno_col='label_5')
    nhood_adata = adata.uns['nhood_adata'].copy()
    nhood_adata.obs.loc[nhood_adata.obs['nhood_annotation_frac'] < 0.75, 'nhood_annotation'] = 'Mixed'
    adata.uns['nhood_adata'] = nhood_adata.copy()
    milopy.plot.plot_DA_beeswarm(adata)
    plt.show()


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
sc.pl.DotPlot(ad_glu,groupby='label_subtype',var_names=['MEIS2','MEIS1','CLMN','DPYD','PCSK5',\
                 'LMX1A','SORCS3','FOXP2'],categories_order=categories_order,use_raw=False,ax=ax).style(cmap='GnBu',largest_dot=250,dot_edge_lw=0).show()
fig.tight_layout()
fig.savefig('./figures/gluts_dotplot.pdf')

ad_glu = milo_cell_rate(ad_glu)
milo_outs = ad_glu.uns['nhood_adata'].obs
print(milo_outs.groupby('nhood_annotation').mean())
milo_outs['relative_to_mixed_logFC'] = milo_outs['logFC']-milo_outs.groupby('nhood_annotation').mean().loc['Mixed','logFC']
milo_outs['is_significant']= (milo_outs['PValue']<0.05).astype(str)
order = ['Ex_MEIS2','Ex_MEIS1','Ex_CLMN','Ex_DPYD','Ex_PCSK5','Ex_LMX1A']
fig,ax = plt.subplots(figsize=(7,5),dpi=300,constrained_layout=True)
sns.violinplot(data=milo_outs, y='nhood_annotation', x="relative_to_mixed_logFC",
                       size=190, inner=None, orient="h",order=order,
                       linewidth=0,
                       scale="width",ax=ax)
sns.stripplot(data=milo_outs, y='nhood_annotation', x="relative_to_mixed_logFC", size=2,order=order,
                  orient="h",hue='is_significant', palette=['grey', 'black'],alpha=0.6,ax=ax)
plt.axvline(x=0, ymin=0, ymax=1, color="black", linestyle="--",linewidth=0.5)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_linewidth(0.5)
#ax.yaxis.set_major_formatter(plt.NullFormatter())
#ax.yaxis.set_minor_formatter(plt.NullFormatter())
ax.spines['left'].set_visible(False)
fig.savefig('Figures/Gluts_abundance.pdf')

