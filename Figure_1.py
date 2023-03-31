from matplotlib.figure import Figure
import scanpy as sc
import pandas as pd
import numpy as np
import os
from matplotlib import *
import matplotlib.pyplot as plt
#Define a nice colour map for gene expression
colors2 = plt.cm.Reds(np.linspace(0, 0.7, 100))
colors3 = plt.cm.Greys_r(np.linspace(0.6,1,100))
colorsComb = np.vstack([colors3, colors2])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)
rcParams['pdf.fonttype']=42 # in order to save fonttype in AI
rcParams['ps.fonttype']=42

# define parameters
color_dict={
    'DaNs': '#BF8219','Excitatory': '#00835A','GABA': '#BF480D','HdNs': '#FF6600','Inhibitory': '#FFB307','Endothelial': '#A19922',
    'Ependymal': '#FFA388','Microglia': '#D64849','Astrocytes': '#6c00bf','OPC': '#74A0FF','NFOL': '#69A8E6','MOL': '#0000ff','Vasc': '#807B30',
    'SMC/VLMC':'#9A9a9a','doublets': '#d6d6d6','low_quality':'#669D6A','Immune': '#f252c5','CHAT IN':'#ff7f0e','PVALB/TH IN':'#c49c94',
    'SST/NPY IN':'#ffbb78','dSPN':'#98df8a','eSPN':'#1f77b4','iSPN':'#279e68'
}
snc_vmax_vmin_dict = {
    'AQP4':(3.5,0),'TH':(1,0),'GAD2':(2,0),'GAD1':(2,0),'SLC17A6':(1,0),'MOG':(3,0),'PDGFRA':(2.5,0),'DCN':(2.5,0),'FLT1':(4,0),'HDC':(3,0),'C1QA':(3,0),
    'CSF1R':(3,0),'LSP1':(2,0),'CD200R1':(1,0),'CD3G':(2,0),'SLC6A3':(2,0),'CSPG4':(1.5,0),'VCAN':(4,0),'BCAS1':(3,0)
}
putamen_vmax_vmin_dict = {
    'AQP4':(3,0),'GAD2':(2.5,0),'LHX6':(2,0),'CHAT':(2,0),'SLC17A7':(1,0),'DRD2':(2.5,0),'DRD1':(0.5,0),'OTOF':(2,0),'FOXP2':(3,0),'VCAN':(3.5,0),'PDGFRA':(2.5,0),
   'GAD2':(2.5,0),'MOG':(3.5,0),'BCAS1':(3.5,0),'CSF1R':(3,0),'CD3G':(2,0),'CD200R1':(1,0),'LSP1':(2,0),'DCN':(3,0),'FLT1':(3.5,0),'NOS1':(4.5,0),'CNTNAP4':(3,0),
}
# define plot function
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

#plot umap
def plot_umap(ad,sample_area):
    fig,ax = plt.subplots(figsize=(5,5),dpi=600)
    plt.rcParams['font.sans-serif']='Arial'
    plt.rcParams['font.size']=7
    plt.rcParams['font.weight']='normal'
    sc.pl.umap(ad,color='label_5',ax=ax,legend_loc='',frameon=False,title='',s=2.5)
    fig.savefig(f'./Figures1/{sample_area}_umap.png')

# plot cell_No. and cell_gene counts
def plot_cell_proportion(ad_snc,sample_area):
    fig = plt.figure(figsize=(8.5,4.5),dpi=300)
    grid = plt.GridSpec(1,9,hspace=0.2,wspace=0.1)
    plt.gcf().subplots_adjust(bottom=0.2)
    #ax_dot = fig.add_subplot(grid[0,0:3])
    ax_cell_number = fig.add_subplot(grid[0,4:5])
    ax_genes = fig.add_subplot(grid[0,5:7])
    ad_cell_rate = fig.add_subplot(grid[0,7:9])
    cell_cluster_order = ad_snc.obs['label_5'].value_counts().index
    expression = pd.DataFrame(ad_snc.layers['counts'].A,index = ad_snc.obs_names,columns = ad_snc.var_names)
    expression['label_5']=ad_snc.obs['label_5']
    filter_exp = expression.groupby('label_5').apply(lambda x:(x>1).sum(0))
    detected_gene = (filter_exp>1).sum(1)
    cell_type = detected_gene.index.tolist()
    cell_count=ad_snc.obs.groupby('label_5').size()
    cell_count = cell_count[cell_cluster_order]
    size = ad_snc.obs.groupby('label_5').size()
    props = size*100/size.sum()
    props = props[cell_cluster_order]
    # plot cell counts
    plt.rcParams['font.sans-serif']='Arial'
    plt.rcParams['font.size']=7
    plt.rcParams['font.weight']='normal'
    with plt.style.context('seaborn-paper'):
        ax_cell_number.set_yticks(np.arange(1,len(cell_type)+1))
        ax_cell_number.set_yticklabels(cell_cluster_order)
        ax_cell_number.set_ylim(0.5,len(cell_type)+1)
        ax_cell_number.set_xticks([0.25])
        ax_cell_number.set_xlim(0,0.3)
        ax_cell_number.xaxis.set_major_formatter(plt.NullFormatter())
        ax_cell_number.xaxis.set_major_locator(plt.NullLocator())
        ax_cell_number.set_xlabel('No. of cells')
        ax_cell_number.spines['right'].set_visible(False)
        ax_cell_number.spines['top'].set_visible(False)
        counter =0
        for i in cell_count:
            ax_cell_number.text(0.1,counter+0.85,str(i))
            counter+=1
    #plot cell genes
    plt.rcParams['font.sans-serif']='Arial'
    plt.rcParams['font.size']=7
    plt.rcParams['font.weight']='normal'
    with plt.style.context('seaborn-paper'):
        ax_genes.set_yticks(np.arange(1,len(cell_type)+1))
        ax_genes.set_yticklabels(cell_cluster_order)
        ax_genes.set_ylim(0.5,len(cell_type)+1)
        ax_genes.yaxis.set_major_formatter(plt.NullFormatter())
        #ax.yaxis.set_major_locator(plt.NullLocator())
        ax_genes.set_xticks([0,5000,10000,15000])
        ax_genes.set_xlabel('Detected genes per\n   cell type')
        ax_genes.spines['right'].set_visible(False)
        ax_genes.spines['top'].set_visible(False)
        counter =1
        for i in cell_cluster_order:
            ax_genes.barh(counter,detected_gene[i],color=color_dict.get(i))
            counter+=1
    
    with plt.style.context('seaborn'):
        sum_ex = 0
        counter =0
        for i in props.index:
            ad_cell_rate.bar(1,props[i],width=0.5,bottom=sum_ex,color=color_dict.get(i))
            pro_index = ('% .1f'% props[i])
            if (props[i] > 8):
                ad_cell_rate.text(0.9,sum_ex+props[i]/2,i)
                ad_cell_rate.text(0.9,sum_ex+props[i]/2+3,f'{pro_index}%')
            else:
                ad_cell_rate.text(1.25,counter,i)
                ad_cell_rate.text(1.25,counter+4,f'{pro_index}%')
            sum_ex = sum_ex+props[i]
            counter+=9
        ad_cell_rate.yaxis.set_major_locator(plt.NullLocator())
        ad_cell_rate.xaxis.set_major_locator(plt.NullLocator())
        ad_cell_rate.spines['bottom'].set_visible(False)#设置边框是否可见
        ad_cell_rate.spines['right'].set_visible(False)
        ad_cell_rate.spines['top'].set_visible(False)
        ad_cell_rate.spines['left'].set_visible(False)
    fig.tight_layout()
    fig.savefig(f'./Figures1/{sample_area}_cell_structre.pdf')

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
        fig.savefig(f'./Figures1/cell_rate_{sample_area}.pdf')
        fig.savefig(f'./Figures1/cell_rate_{sample_area}.png')    


if os.path.exists('Figures1'):
    pass
else:
    os.mkdir('Figures1')
ad_snc = sc.read('./snc_figout.h5')
color_label_5 = [color_dict.get(i) for i in ad_snc.obs['label_5'].cat.categories]
ad_snc.uns['label_5_colors']=color_label_5
plot_cell_proportion(ad_snc,'snc')
plot_PD_control_proportion(ad_snc[ad_snc.obs['seq_type']=='10x_nuclear',],'snc')
plot_umap(ad_snc,'snc') # plot umap and save fig
umap_position = pd.DataFrame(ad_snc.obsm['X_umap'],index=ad_snc.obs_names,columns=['X','Y'])
exp_matrix = ad_snc.to_df()
for i in snc_vmax_vmin_dict:  # plot each genes and save fig
    plot_umap_genes(umap_position,exp_matrix,'snc',i,snc_vmax_vmin_dict.get(i)[0],snc_vmax_vmin_dict.get(i)[1])
# plot all indicator
for indicator in ['DaNs','Excitatory','Inhibitory','Microglia','Astrocytes','MOL']:
    fig,ax = plt.subplots(figsize=(2,2),dpi=600)
    ax.scatter(umap_position['X'],umap_position['Y'],marker='o',c=[{key:value for key,value in color_dict.items() if key==indicator}.get(x,'grey') for x in ad_snc.obs['label_5']],s=0.2,linewidths=0)
    ax.axis('off')
    fig.savefig(f'./Figures2/{indicator}_umap_indicator.png')
del ad_snc
del umap_position
del exp_matrix    # free memory

#load putamen dataset and plot...
ad_putamen = sc.read('./putamen_figout.h5')
color_label_5 = [color_dict.get(i) for i in ad_putamen.obs['label_5'].cat.categories]
ad_putamen.uns['label_5_colors']=color_label_5
plot_PD_control_proportion(ad_putamen,'putamen')
plot_umap(ad_putamen,'putamen') # plot umap and save fig
umap_position = pd.DataFrame(ad_putamen.obsm['X_umap'],index=ad_putamen.obs_names,columns=['X','Y'])
exp_matrix = ad_putamen.to_df()
for i in putamen_vmax_vmin_dict:  # plot each genes and save fig
    plot_umap_genes(umap_position,exp_matrix,'putamen',i,putamen_vmax_vmin_dict.get(i)[0],putamen_vmax_vmin_dict.get(i)[1])
plot_cell_proportion(ad_putamen,'putamen')
del ad_putamen
del umap_position
del exp_matrix    # free memory
