import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import milopy
import milopy.core as milo
from matplotlib import *
def milo_cell_rate(adata,k=50,d=30,prop=0.3,rep='X_pca'):
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
            plt.legend(title=f'< {int(alpha*100)} % SpatialFDR')
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
adata = sc.read('../All_single_neuclear_cell_SNC_putamen/Figures_pip_ver2/snc_figout.h5')
milo_cell_rate(adata,k=20,d=10,prop=0.3,rep='X_harmonypca')
milo_results = adata.uns["nhood_adata"].obs.copy()
cell_cluster_order = ['DaNs','GABA','HdNs','Excitatory','Inhibitory','Immune','OPC','Endothelial','Astrocytes','MOL','Vasc','Microglia']
color_dict={
    'DaNs': '#BF8219','Excitatory': '#00835A','GABA': '#BF480D','HdNs': '#FF6600','Inhibitory': '#FFB307','Endothelial': '#A19922',
    'Ependymal': '#FFA388','Microglia': '#D64849','Astrocytes': '#6c00bf','OPC': '#74A0FF','NFOL': '#69A8E6','MOL': '#0000ff','Vasc': '#807B30',
    'SMC/VLMC':'#9A9a9a','doublets': '#d6d6d6','low_quality':'#669D6A','Immune': '#f252c5','CHAT IN':'#ff7f0e','PVALB/TH IN':'#c49c94',
    'SST/NPY IN':'#ffbb78','dSPN':'#98df8a','eSPN':'#1f77b4','iSPN':'#279e68'
}
color_order = [color_dict.get(i) for i in cell_cluster_order]
import seaborn as sns
fig,ax = plt.subplots(figsize=(5,4),dpi=300,constrained_layout=True)
sns.violinplot(data=milo_results.loc[milo_results['nhood_annotation']!='Mixed',], y='nhood_annotation', x="logFC",
                       size=190, inner=None, orient="h",palette=color_order,order=cell_cluster_order,
                       linewidth=0,
                       scale="width",ax=ax)
sns.stripplot(data=milo_results.loc[milo_results['nhood_annotation']!='Mixed',], y='nhood_annotation', x="logFC", size=0.5,order=cell_cluster_order,
                  orient="h",hue='is_signif', palette=['grey', 'black'],alpha=0.5,ax=ax)
plt.axvline(x=0, ymin=0, ymax=1, color="black", linestyle="--",linewidth=0.8)
ax.legend_.remove()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
fig.savefig('./sn_maintype_milo_cellrate.pdf')
