import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import *
rcParams['pdf.fonttype']=42
rcParams['ps.fonttype']=42

def plot_type_composition(plot_data,file_path):
    fig,ax = plt.subplots(figsize=(10,5),dpi=200)
    counter=1
    for i in plot_data.index:
        sum_ex = 0
        for j in plot_data.columns:
            ax.barh(counter,width=plot_data.loc[i,j],left=sum_ex,color=color_dict.get(j))
            sum_ex+=plot_data.loc[i,j]
        counter+=1
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    if len(plot_data.index)<5:
        plt.yticks([1,5])
        ax.set_yticklabels(['1','5'])
    else:
        if len(plot_data.index)<=10:
            plt.yticks([1,5,10])
            ax.set_yticklabels(['1','5','10'])
        else:
            plt.yticks([1,5,10,15])
            ax.set_yticklabels(['1','5','10','15'])
    plt.xlim(0,1)
    plt.title('Cell type composition within individual sample')
    fig.savefig(file_path)

ad_sn = sc.read('../data/sn_umap.h5')
nUMI = pd.DataFrame(ad_sn.layers['counts'].sum(1))[0]
fig,ax = plt.subplots(figsize=(3,3),dpi=300)
nbin=300
plt.hist(nUMI,bins=nbin,color='#008000')
plt.xlim(-500,40000)
plt.ylim(0,18000)
mean_gene = int(np.mean(nUMI))
median_gene = int(np.median(nUMI))
plt.plot([mean_gene,mean_gene],[0,16000],color='yellow',label=f'Mean({mean_gene})',linewidth=1) #mean
plt.plot([median_gene,median_gene],[0,16000],color='red',label=f'Median({median_gene})',linewidth=1)#media
plt.ylabel('# Of Cells')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.yticks([0,5000,10000,15000])
ax.set_yticklabels(['0','5','10','15'])
plt.legend(loc='upper right')
plt.text(-500,20000,'*1,000')
plt.title('UMI counts')
fig.savefig('./figures/detected_UMIs.pdf')

n_genes = pd.DataFrame((ad_sn.layers['counts']>0).sum(1))[0]
fig,ax = plt.subplots(figsize=(3,3),dpi=300)
nbin=100
plt.hist(n_genes,bins=nbin,color='#2179B5')
plt.xlim(-500,10000)
plt.ylim(0,9000)
mean_gene = int(np.mean(n_genes))
median_gene = int(np.median(n_genes))
plt.plot([mean_gene,mean_gene],[0,8000],color='yellow',label=f'Mean({mean_gene})',linewidth=1) #mean
plt.plot([median_gene,median_gene],[0,8000],color='red',label=f'Median({median_gene})',linewidth=1)#media
plt.ylabel('# Of Cells')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.yticks([0,2000,4000,6000,8000])
ax.set_yticklabels(['0','2','4','6','8'])
plt.legend(loc='upper right')
plt.text(-500,10000,'*1,000')
plt.title('Detected genes')
fig.savefig('./figures/detected_genes.pdf')

plot_data = pd.DataFrame(columns =ad_sn.obs['label_5'].cat.categories,index=ad_sn.obs['sampleid'].cat.categories )
for i in plot_data.columns:
    each_sample_info = ad_sn[ad_sn.obs['label_5']==i,].obs['sampleid'].value_counts()
    for j in plot_data.index:
        if j in each_sample_info.index:
            plot_data.loc[j,i]=each_sample_info[j]
        else:
            plot_data.loc[j,i]=0
for i in plot_data.index:
    plot_data.loc[i,] = plot_data.loc[i,]/plot_data.sum(1)[i]
color_dict={
    'DaNs': '#BF8219','Gluts': '#00835A','GABA': '#BF480D','HdNs': '#FF6600','Inhibitory': '#FFB307','Endothelial': '#A19922',
    'Ependymal': '#FFA388','Microglia': '#D64849','Astrocytes': '#6c00bf','OPC': '#74A0FF','Oligo': '#0000ff','Vasc': '#807B30',
    'SMC/VLMC':'#9A9a9a','doublets': '#d6d6d6','Immune': '#f252c5','CHAT IN':'#ff7f0e','PVALB/TH IN':'#c49c94',
    'SST/NPY IN':'#ffbb78','dSPN':'#98df8a','eSPN':'#1f77b4','iSPN':'#279e68'
}
single_cell_samples = ['SN_112801', 'SN_112802', 'SN_112803', 'SN_112804']
sn_nuclei = plot_data.loc[~plot_data.index.isin(single_cell_samples),]
sn_cell = plot_data.loc[plot_data.index.isin(single_cell_samples),]
plot_type_composition(sn_nuclei,'./figures/nuclei_type_composition.pdf')
plot_type_composition(sn_cell,'./figures/cell_type_composition.pdf')

each_sample_info = pd.read_csv('../data/each_sample_info.csv',index_col=0)
color_dict = {'SNC':'#98df8a','Putaman':'#279e68','control':'#2E72A4','PD':'#BC2823'}
fig,ax = plt.subplots(6,1,figsize=(4,8),dpi=300)
counter=0
for i in plot_data.columns:
    ax[5].bar(counter,each_sample_info.loc['n_genes',i],color='grey',width=1,edgecolor='white',linewidth=0.1)
    ax[5].set_xlim(-0.8,31.3)
    ax[5].spines['right'].set_visible(False)
    ax[5].spines['top'].set_visible(False)
    ax[5].set_yticks([0,15000])
    ax[5].xaxis.set_major_formatter(plt.NullFormatter())
    ax[5].xaxis.set_major_locator(plt.NullLocator())
    ax[5].set_ylabel('n_genes',rotation='horizontal')
    ax[4].bar(counter,each_sample_info.loc['cell_counts',i],color='grey',width=1,edgecolor='white',linewidth=0.1)
    ax[4].set_xlim(-0.8,31.3)
    ax[4].spines['right'].set_visible(False)
    ax[4].spines['top'].set_visible(False)
    ax[4].set_yticks([0,15000])
    ax[4].xaxis.set_major_formatter(plt.NullFormatter())
    ax[4].xaxis.set_major_locator(plt.NullLocator())
    ax[4].set_ylabel('cell_counts',rotation='horizontal')
    ax[3].bar(counter,each_sample_info.loc['low quaility cells',i],color='grey',width=1,edgecolor='white',linewidth=0.1)
    ax[3].set_xlim(-0.8,31.3)
    ax[3].spines['right'].set_visible(False)
    ax[3].spines['top'].set_visible(False)
    ax[3].set_yticks([0,15000])
    ax[3].xaxis.set_major_formatter(plt.NullFormatter())
    ax[3].xaxis.set_major_locator(plt.NullLocator())
    ax[3].set_ylabel('low quaility cells',rotation='horizontal')
    ax[2].bar(counter,each_sample_info.loc['doublets',i],color='grey',width=1,edgecolor='white',linewidth=0.1)
    ax[2].set_xlim(-0.8,31.3)
    ax[2].spines['right'].set_visible(False)
    ax[2].spines['top'].set_visible(False)
    ax[2].set_yticks([0,15000])
    ax[2].xaxis.set_major_formatter(plt.NullFormatter())
    ax[2].xaxis.set_major_locator(plt.NullLocator())
    ax[2].set_ylabel('doublets',rotation='horizontal')
    ax[1].bar(counter,1,color=color_dict.get(each_sample_info.loc['Condition',i]),width=1,edgecolor='white',linewidth=0)
    ax[1].axis('off')
    ax[1].set_xlim(-0.8,31.3)
    ax[0].bar(counter,1,color=color_dict.get(each_sample_info.loc['Brain region',i]),width=1,edgecolor='white',linewidth=0)
    ax[0].axis('off')
    ax[0].set_xlim(-0.8,31.3)
    counter+=1
fig.savefig('./Figures/sample_info.pdf')







