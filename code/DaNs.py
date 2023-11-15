import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import *
rcParams['pdf.fonttype']=42
rcParams['ps.fonttype']=42
import seaborn as sns
import SCCAF
import milopy
import milopy.core as milo

def plot_PD_control_proportion(ad,sample_area):
    cell_rate = pd.crosstab(ad.obs['Condition'],ad.obs['label_subtype'])
    all_cells = ad.obs['Condition'].value_counts()
    plot_cell_rate = pd.DataFrame(index=cell_rate.index,columns=cell_rate.columns)
    for i in cell_rate.columns:
        plot_cell_rate[i]=cell_rate[i]/all_cells
    plot_cell_rate = plot_cell_rate/plot_cell_rate.sum(0)
    plot_cell_rate = plot_cell_rate.loc[:,['SOX6-_SORCS3', 'SOX6-_COBLL1', 'SOX6-_CBLB', 'SOX6-_CALCR', 'SOX6-_EBF1', 'SOX6+_CUX2', 'SOX6+_RYR3']]
    plot_cell_rate = plot_cell_rate.loc[['control','MPTP'],]
    plt.rcParams['font.size']=7
    plt.rcParams['font.weight']='normal'
    sample_color={'control':'#2E72A4','MPTP':'#BC2823'}
    fig,ax = plt.subplots(figsize=(3.5,3),dpi=300,constrained_layout=True)
    counter=1
    for i in plot_cell_rate.columns:
        width=0
        for j in plot_cell_rate.index:
            ax.barh(counter,plot_cell_rate.loc[j,i],left=width,color=sample_color.get(j))
            width+=plot_cell_rate.loc[j,i]
        counter+=1
        ax.set_yticks(range(1,len(plot_cell_rate.columns)+1))
        ax.set_yticklabels(plot_cell_rate.columns.tolist())
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_xticks([0,0.5,1])
        ax.xaxis.set_major_formatter(plt.NullFormatter())
        ax.tick_params(width=0)
        ax.set_xticklabels(['0','50','100'])
        ax.set_ylim([0.5,len(plot_cell_rate.columns)+0.5])
        ax.xaxis.set_ticks_position('top') 
        ax.set_title('% Cells',fontsize=14)
    plt.show()
    fig.savefig('./Figures/DaNs_cell_rete.pdf')
    return plot_cell_rate

def milo_cell_rate(adata,k=10,d=5,prop=0.3,rep='X_harmonypca'):
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
    milo.DA_nhoods(adata, design="~Condition", model_contrasts='ConditionMPTP-Conditioncontrol')
    milo_results_salmonella = adata.uns["nhood_adata"].obs.copy()
    def plot_milo_diagnostics(adata):
        nhood_adata = adata.uns['nhood_adata'].copy()
        alpha = 0.5 ## significance threshold  
        with matplotlib.rc_context({'figure.figsize':[12,12]}):    
            ## Check P-value histogram
            plt.subplot(2,2,1);
            plt.hist(nhood_adata.obs['PValue'], bins=20);
            plt.xlabel('Uncorrected P-value');   
            ## Visualize extent of multiple-testing correction
            plt.subplot(2,2,2);
            plt.scatter(nhood_adata.obs['PValue'], nhood_adata.obs['SpatialFDR'], s=3);
            plt.xlabel('Uncorrected P-value');
            plt.ylabel('SpatialFDR');    
            ## Visualize volcano plot
            plt.subplot(2,2,3);
            plt.scatter(nhood_adata.obs['logFC'], -np.log10(nhood_adata.obs['SpatialFDR']), s=3);
            plt.axhline(y=-np.log10(alpha), color='red', linewidth=1, label=f'{int(alpha*100)} % SpatialFDR');
            plt.legend();
            plt.xlabel('log-Fold Change');
            plt.ylabel('- log10(SpatialFDR)');
            plt.tight_layout()    
            ## Visualize MA plot
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
        sc.pl.umap(adata, color='label_7', legend_loc='on data')
    milopy.utils.annotate_nhoods(adata, anno_col='label_7')
    # Define as mixed if fraction of cells with same label < 0.75
    nhood_adata = adata.uns['nhood_adata'].copy()
    nhood_adata.obs.loc[nhood_adata.obs['nhood_annotation_frac'] < 0.75, 'nhood_annotation'] = 'Mixed'
    adata.uns['nhood_adata'] = nhood_adata.copy()
    milopy.plot.plot_DA_beeswarm(adata)
    plt.show()
    return adata

ad_dans = sc.read('../data/DaNs_umap.h5')
colors2 = plt.cm.Oranges(np.linspace(0.2, 1, 7))
sns.palplot(sns.color_palette(colors2,7))
color_list = [colors.to_hex(i) for i in colors2]
order = ['SOX6+_RYR3','SOX6+_CUX2','SOX6-_EBF1','SOX6-_CBLB','SOX6-_CALCR','SOX6-_COBLL1','SOX6-_SORCS3']
color_order = color_list
color_dict = dict(zip(order,color_order))
ad_dans.uns['label_subtype_colors'] = [color_dict.get(i) for i in ad_dans.obs['label_subtype'].cat.categories]
fig,ax = plt.subplots(figsize=(4,3.5),dpi=600)
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.size']=7
plt.rcParams['font.weight']='normal'
sc.pl.umap(ad_dans,color='label_subtype',ax=ax,legend_loc='',frameon=False,title='',s=120)
fig.savefig('./figures/DaNs_umap.pdf')
categories_order = ['SOX6+_RYR3','SOX6+_CUX2','SOX6-_EBF1','SOX6-_CALCR','SOX6-_CBLB','SOX6-_COBLL1','SOX6-_SORCS3']
fig,ax = plt.subplots(figsize=(4.1,3),dpi=600)
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.size']=7
plt.rcParams['font.weight']='normal'
sc.pl.DotPlot(ad_dans,groupby='label_subtype',var_names=['SOX6','KCNJ6','RYR3','CUX2','EBF1','CALCR','CBLB','COBLL1','SORCS3','FOXP2'],categories_order=categories_order,use_raw=False,ax=ax).style(cmap='GnBu',dot_edge_lw=0).show()
fig.tight_layout()
fig.savefig('./figures/DaNs_dotplot.pdf')
plot_cell_rate = plot_PD_control_proportion(ad_dans,'sn')
#we observed merged human and monkey datasets and choose their common hvg for further analysis 
ad_homo_da = sc.read('../data/homo_DaNs_umap.h5')
ad_mk_da = sc.read('../data/DaNs_umap.h5')
ad_homo_da = ad_homo_da[:,ad_homo_da.var_names.isin(ad_mk_da.var_names)]
ad_mk_da = ad_mk_da[:,ad_homo_da.var_names]
ad_homo_da.obs['org'] = 'human'
ad_mk_da.obs['org']='monkey'
sc.pp.highly_variable_genes(ad_homo_da, flavor="seurat_v3", n_top_genes=1500,layer='counts')
monkey_hvgs = ad_mk_da.var_names[ad_mk_da.var['highly_variable']]
human_hvg = ad_homo_da.var_names[ad_homo_da.var['highly_variable']]
sele_gene = monkey_hvgs[monkey_hvgs.isin(human_hvg)].tolist()
ad_mk_da = ad_mk_da[:,sele_gene]
ad_homo_da = ad_homo_da[:,sele_gene]
ad = ad_homo_da.concatenate(ad_mk_da)
sc.tl.pca(ad, use_highly_variable = False)
sc.pp.neighbors(ad)
sc.tl.umap(ad, maxiter=200, min_dist=1, spread=2)
ho = hm.run_harmony(ad.obsm['X_pca'], ad.obs, ['org'])
ad.obsm['X_harmonypca'] = ho.Z_corr.T
sc.pp.neighbors(ad, use_rep='X_harmonypca')
ad.obsm['X_umapraw'] = ad.obsm['X_umap']
sc.tl.umap(ad)
ad.obsm['X_umapharmony'] = ad.obsm['X_umap']
ad_mk = ad[ad.obs['org']=='monkey',]
ad_mk = milo_cell_rate(ad_mk)
milo_outs = ad_mk.uns['nhood_adata'].obs
print(milo_outs.groupby('nhood_annotation').mean())
milo_outs['relative_to_mixed_logFC'] = milo_outs['logFC']-milo_outs.groupby('nhood_annotation').mean().loc['Mixed','logFC']
milo_outs['is_significant']= (milo_outs['PValue']<0.05).astype(str)
order = ['SOX6+_RYR3','SOX6+_CUX2','SOX6-_EBF1','SOX6-_CALCR','SOX6-_CBLB','SOX6-_COBLL1','SOX6-_SORCS3']
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
fig.savefig('Figures/dans_abundance.pdf')




