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
  
ad_snc = sc.read('./scvi_sorted_fx_refined_umap.h5')
nUMI = pd.DataFrame(ad_snc.layers['counts'].sum(1))[0]
fig,ax = plt.subplots(figsize=(3,3),dpi=300)
nbin=200
plt.hist(nUMI,bins=nbin,color='#008000')
plt.xlim(-500,40000)
plt.ylim(0,5000)
mean_gene = int(np.mean(nUMI))
median_gene = int(np.median(nUMI))
plt.plot([mean_gene,mean_gene],[0,16000],color='yellow',label=f'Mean({mean_gene})',linewidth=1) #mean
plt.plot([median_gene,median_gene],[0,16000],color='red',label=f'Median({median_gene})',linewidth=1)#media
plt.ylabel('# Of Cells')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.yticks([0,2500,5000])
ax.set_yticklabels(['0','2.5','5'])
plt.legend(loc='upper right')
plt.text(-500,5000,'*1,000')
plt.title('UMI counts')
fig.savefig('./figures/FASC_sample_info_detected_UMIs.pdf')

n_genes = pd.DataFrame((ad_snc.layers['counts']>0).sum(1))[0]
fig,ax = plt.subplots(figsize=(3,3),dpi=300)
nbin=90
plt.hist(n_genes,bins=nbin,color='#2179B5')
plt.xlim(-500,10000)
plt.ylim(0,1000)
mean_gene = int(np.mean(n_genes))
median_gene = int(np.median(n_genes))
plt.plot([mean_gene,mean_gene],[0,8000],color='yellow',label=f'Mean({mean_gene})',linewidth=1) #mean
plt.plot([median_gene,median_gene],[0,8000],color='red',label=f'Median({median_gene})',linewidth=1)#media
plt.ylabel('# Of Cells')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.yticks([0,500,1000])
ax.set_yticklabels(['0','1','2'])
plt.legend(loc='upper right')
plt.text(-500,1000,'*1,000')
plt.title('Detected genes')
fig.savefig('./figures/FASC_sample_info_detected_genes.pdf')

plot_data = pd.DataFrame(columns =ad_snc.obs['label_main'].cat.categories,index=ad_snc.obs['sampleid'].cat.categories )
for i in plot_data.columns:
    each_sample_info = ad_snc[ad_snc.obs['label_main']==i,].obs['sampleid'].value_counts()
    for j in plot_data.index:
        if j in each_sample_info.index:
            plot_data.loc[j,i]=each_sample_info[j]
        else:
            plot_data.loc[j,i]=0
for i in plot_data.index:
    plot_data.loc[i,] = plot_data.loc[i,]/plot_data.sum(1)[i]
plot_order = ['Astrocytes','DaNs','Endothelial','Gluts','GABA','HdNs','Inhibitory','Oligo','Microglia','OPC','Vasc']
plot_data = plot_data.loc[:,plot_order]
color_dict={
    'DaNs': '#BF8219','Gluts': '#00835A','GABA': '#BF480D','HdNs': '#FF6600','Inhibitory': '#FFB307','Endothelial': '#A19922',
    'Ependymal': '#FFA388','Microglia': '#D64849','Astrocytes': '#6c00bf','OPC': '#74A0FF','NFOL': '#69A8E6','Oligo': '#0000ff','Vasc': '#807B30',
    'SMC/VLMC':'#9A9a9a','doublets': '#d6d6d6','low_quality':'#669D6A','Immune': '#f252c5','CHAT IN':'#ff7f0e','PVALB/TH IN':'#c49c94',
    'SST/NPY IN':'#ffbb78','dSPN':'#98df8a','eSPN':'#1f77b4','iSPN':'#279e68'
}
fig_fasc = plot_type_composition(plot_data,'./figures/fig1_extend_nuclei_type_composition.pdf')

fig,ax = plt.subplots(figsize=(5,5),dpi=600)
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.size']=7
plt.rcParams['font.weight']='normal'
sc.pl.umap(ad_snc,color='sampleid',ax=ax,frameon=False,title='',s=2.5)
fig.savefig(f'./figures/fasc_sampleid_umap.pdf')
