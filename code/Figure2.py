import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import *
rcParams['pdf.fonttype']=42
rcParams['ps.fonttype']=42
import seaborn as sns

def plot_PD_control_proportion(ad,sample_area):
    cell_rate = pd.crosstab(ad.obs['exp_condition1'],ad.obs['label_9'])
    all_cells = ad.obs['exp_condition1'].value_counts()
    plot_cell_rate = pd.DataFrame(index=cell_rate.index,columns=cell_rate.columns)
    for i in cell_rate.columns:
        plot_cell_rate[i]=cell_rate[i]/all_cells
    plot_cell_rate = plot_cell_rate/plot_cell_rate.sum(0)
    plot_cell_rate = plot_cell_rate.loc[:,['SOX6-_SORCS3', 'SOX6-_COBLL1', 'SOX6-_CBLB', 'SOX6-_CALCR', 'SOX6-_EBF1', 'SOX6+_CUX2', 'SOX6+_RYR3']]
    plot_cell_rate = plot_cell_rate.loc[['control','PD'],]
    plt.rcParams['font.size']=7
    plt.rcParams['font.weight']='normal'
    sample_color={'control':'#2E72A4','PD':'#BC2823'}
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
    fig.savefig('./Figures_pip_ver3/Figure2/DaNs_cell_rete.pdf')
    return plot_cell_rate

ad_dans = sc.read('./Figures_pip_ver2/DaNs_figout.h5')
colors2 = plt.cm.Oranges(np.linspace(0.2, 1, 7))
sns.palplot(sns.color_palette(colors2,7))
color_list = [colors.to_hex(i) for i in colors2]
order = ['SOX6+_RYR3','SOX6+_CUX2','SOX6-_EBF1','SOX6-_CBLB','SOX6-_CALCR','SOX6-_COBLL1','SOX6-_SORCS3']
color_order = color_list
color_dict = dict(zip(order,color_order))
ad_dans.uns['label_9_colors'] = [color_dict.get(i) for i in ad_dans.obs['label_9'].cat.categories]
fig,ax = plt.subplots(figsize=(4,3.5),dpi=600)
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.size']=7
plt.rcParams['font.weight']='normal'
sc.pl.umap(ad_dans,color='label_9',ax=ax,legend_loc='',frameon=False,title='',s=120)
fig.savefig(f'./Figures_pip_ver3/Figure2/DaNs_umap.pdf')
categories_order = ['SOX6+_RYR3','SOX6+_CUX2','SOX6-_EBF1','SOX6-_CALCR','SOX6-_CBLB','SOX6-_COBLL1','SOX6-_SORCS3']
fig,ax = plt.subplots(figsize=(4.1,3),dpi=600)
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.size']=7
plt.rcParams['font.weight']='normal'
sc.pl.DotPlot(ad_dans,groupby='label_9',var_names=['SOX6','KCNJ6','RYR3','CUX2','EBF1','CALCR','CBLB','COBLL1','SORCS3','FOXP2'],categories_order=categories_order,use_raw=False,ax=ax).style(cmap='GnBu',dot_edge_lw=0).show()
fig.tight_layout()
fig.savefig('./Figures_pip_ver3/Figure2/DaNs_dotplot.pdf')
plot_cell_rate = plot_PD_control_proportion(ad_dans,'snc')
