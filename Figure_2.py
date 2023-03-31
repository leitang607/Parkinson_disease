import re
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import *
import matplotlib.pyplot as plt
#Define a nice colour map for gene expression
colors2 = plt.cm.Reds(np.linspace(0, 0.7, 100))
colors3 = plt.cm.Greys_r(np.linspace(0.6,1,100))
colorsComb = np.vstack([colors3, colors2])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)
rcParams['pdf.fonttype']=42 # in order to save fonttype in AI
rcParams['ps.fonttype']=42
# plot PD control proportion
def plot_PD_control_proportion(ad,sample_area):
    cell_rate = pd.crosstab(ad.obs['exp_condition1'],ad.obs['label_5'])
    all_cells = ad.obs['exp_condition1'].value_counts()
    plot_cell_rate = pd.DataFrame(index=cell_rate.index,columns=cell_rate.columns)
    for i in cell_rate.columns:
        plot_cell_rate[i]=cell_rate[i]/all_cells
    plot_cell_rate = plot_cell_rate/plot_cell_rate.sum(0)
    plot_cell_rate = plot_cell_rate.sort_values(by='control',axis=1)
    plot_cell_rate = plot_cell_rate.loc[['control','PD'],]
    plt.rcParams['font.sans-serif']='Arial'
    plt.rcParams['font.size']=7
    plt.rcParams['font.weight']='normal'
    sample_color={'control':'#86C166','PD':'#EF7A82'}
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
        fig.savefig(f'./Figures2/cell_rate_{sample_area}.pdf')
        fig.savefig(f'./Figures2/cell_rate_{sample_area}.png')

#load data
ad_dans = sc.read('./DaNs_figout.h5')
colors_dic = {'SOX6+_LIN7A': '#60A726',
 'SOX6-_PLCB4': '#1f77b4',
 'SOX6+_PLEKHG1': '#C0417A',
 'SOX6-_SORCS3': '#D86900',
 'SOX6-_ARPP21': '#9467bd',
 'SOX6-_CBLB': '#ff7f0e'}
ad_dans.uns['label_5_colors']=[colors_dic.get(i) for i in ad_dans.obs['label_5'].cat.categories]
# plot umap
fig,ax = plt.subplots(figsize=(5,4),dpi=600)
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.size']=7
plt.rcParams['font.weight']='normal'
sc.pl.umap(ad_dans,color='label_5',ax=ax,legend_loc='',frameon=False,title='',s=70)
fig.savefig(f'./Figures2/DaNs_umap.png')
#plot dotplot
fig,ax = plt.subplots(figsize=(4,4),dpi=600)
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.size']=7
plt.rcParams['font.weight']='normal'
categories_order = ['SOX6+_PLEKHG1','SOX6+_LIN7A','SOX6-_PLCB4','SOX6-_ARPP21','SOX6-_CBLB','SOX6-_SORCS3']
sc.pl.dotplot(ad_dans,groupby='label_5',var_names=['SOX6','PLEKHG1','LIN7A',\
 'PLCB4','ARPP21','CBLB','SORCS3','FOXP2'],categories_order=categories_order,use_raw=False,swap_axes=True,ax=ax)
fig.tight_layout()
fig.savefig('./Figures2/DaNs_dotplot.pdf')
#plot PD proportiion
plot_PD_control_proportion(ad_dans,'DaNs')
#plot regulon heatmap
plot_data = pd.read_csv('./DaNs_regulon_heatmap.csv',index_col=0)
plot_data = plot_data.T
regulon_order = ['Regulon(PBX1)','Regulon(NRF1)','Regulon(ZNF91)','Regulon(FOXP1)','Regulon(MEF2A)','Regulon(RORA)','Regulon(ZEB1)', 'Regulon(SREBF2)','Regulon(MEF2D)',\
                'Regulon(LIN28B)','Regulon(ATF1)','Regulon(ZMAT4)','Regulon(FOXO3)','Regulon(CREB5)','Regulon(FOXP2)', 'Regulon(SMAD9)','Regulon(BCL11A)',\
                'Regulon(MLX)', 'Regulon(TFAP2A)', 'Regulon(NFATC2)']
cluster_order = ['SOX6+_PLEKHG1', 'SOX6+_LIN7A', 'SOX6-_PLCB4', 'SOX6-_CBLB', 'SOX6-_ARPP21', 'SOX6-_TSPAN18', 'SOX6-_ZFHX2']
plot_data = plot_data.loc[regulon_order,cluster_order]
regulon_heatmap = sns.clustermap(plot_data,standard_scale=0,col_cluster=False,row_cluster=False,cmap=mymap,figsize=(4,9))
fig = regulon_heatmap.fig
fig.savefig('./Figures2/DaNs_regulon_heatmap.pdf')

#For exc
ad_exc = sc.read('./exc_figout.h5')
color_dict = {'Ex_CLMN': '#B0D0E5',
 'Ex_DPYD': '#AFDF64',
 'Ex_LMX1A': '#F86160',
 'Ex_MEIS1': '#F79023',
 'Ex_MEIS2': '#A64E2C',
 'Ex_PCSK5': '#873DB1'}
ad_exc.uns['label_5_colors']=[color_dict.get(i) for i in ad_exc.obs['label_5'].cat.categories]
# plot umap
fig,ax = plt.subplots(figsize=(5,4),dpi=600)
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.size']=7
plt.rcParams['font.weight']='normal'
sc.pl.umap(ad_exc,color='label_5',ax=ax,legend_loc='',frameon=False,title='',s=30)
fig.savefig(f'./Figures2/exc_umap.png')
#plot dotplot
fig,ax = plt.subplots(figsize=(4,4),dpi=600)
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.size']=7
plt.rcParams['font.weight']='normal'
categories_order = ['Ex_PCSK5','Ex_MEIS1','Ex_MEIS2','Ex_LMX1A','Ex_DPYD','Ex_CLMN']
sc.pl.dotplot(ad_exc,groupby='label_5',var_names=['PCSK5','KHDRBS2','MEIS1','MEIS2','LMX1A','DPYD','GRID2','CLMN'],categories_order=categories_order,use_raw=False,swap_axes=True,ax=ax)
fig.tight_layout()
fig.savefig('./Figures2/exc_dotplot.pdf')
#plot PD proportiion
plot_PD_control_proportion(ad_exc,'exc')
#plot regulon heatmap
#plot regulon heatmap
plot_data = pd.read_csv('./exc_regulon_heatmap.csv',index_col=0)
plot_data = plot_data.T
regulon_heatmap = sns.clustermap(plot_data,standard_scale=0,col_cluster=False,row_cluster=False,cmap=mymap,figsize=(4,9))
fig = regulon_heatmap.fig
fig.savefig('./Figures2/exc_regulon_heatmap.pdf')
 # For interneuron
ad_inh = sc.read('./inh_figout.h5')
color_dict = {'Int_NOVA1': '#e377c2',
 'Int_DPYD': '#ff7f0e',
 'Int_GRIK1': '#69A891',
 'Int_MEIS2': '#d62728',
 'Int_EBF3': '#4E9F40'}
ad_inh.uns['label_5_colors']=[color_dict.get(i) for i in ad_inh.obs['label_5'].cat.categories]
#plot umap
fig,ax = plt.subplots(figsize=(5,4),dpi=600)
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.size']=7
plt.rcParams['font.weight']='normal'
sc.pl.umap(ad_inh,color='label_5',ax=ax,legend_loc='',frameon=False,title='',s=30)
fig.savefig(f'./Figures2/inh_umap.png')
#plot dotplot
fig,ax = plt.subplots(figsize=(4,4),dpi=600)
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.size']=7
plt.rcParams['font.weight']='normal'
categories_order = ['Int_NOVA1','Int_DPYD','Int_GRIK1','Int_MEIS2','Int_EBF3']
sc.pl.dotplot(ad_inh,groupby='label_5',var_names=['NOVA1','DPYD','GRIK1','MEIS2','EBF3'],categories_order=categories_order,use_raw=False,swap_axes=True,ax=ax)
fig.tight_layout()
fig.savefig('./Figures2/inh_dotplot.pdf')
#plot PD proportiion
plot_PD_control_proportion(ad_inh,'inh')