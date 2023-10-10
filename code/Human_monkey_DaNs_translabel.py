
import pandas as pd
import numpy as np
import scanpy as sc
import harmonypy as hm
import os,sys
import seaborn as sns
import matplotlib.pyplot as plt
import SCCAF
from matplotlib import *
rcParams['pdf.fonttype']=42
rcParams['ps.fonttype']=42

def sele_hvg(ad_mk,ad_hu,n_hvg=1500):
    ad_mk.X = ad_mk.layers['counts'].copy()
    sc.pp.normalize_total(ad_mk, target_sum=1e4)
    sc.pp.log1p(ad_mk)
    sc.pp.highly_variable_genes(ad_mk, flavor="seurat_v3", n_top_genes=n_hvg,layer='counts')
    monkey_hvgs = ad_mk.var_names[ad_mk.var['highly_variable']]
    ad_hu.X = ad_hu.layers['counts'].copy()
    sc.pp.normalize_total(ad_hu, target_sum=1e4)
    sc.pp.log1p(ad_hu)
    sc.pp.highly_variable_genes(ad_hu, flavor="seurat_v3", n_top_genes=n_hvg,layer='counts')
    human_hvg = ad_hu.var_names[ad_hu.var['highly_variable']]
    sele_gene_intersection =  monkey_hvgs[monkey_hvgs.isin(human_hvg)].tolist()
    sele_gene_intergration = list(set(monkey_hvgs.tolist()+human_hvg.tolist()))
    return sele_gene_intersection, sele_gene_intergration

def sccaf_predict(ad_mk_sele,ad_hu_sele):
    ad_mk_sele.obs = ad_mk_sele.obs.loc[:,~ad_mk_sele.obs.columns.str.contains('Regulon')]
    ad_mk_sele.obs['donor_id'] = ad_mk_sele.obs['sampleid']
    ad_mk_sele.obs['org']='monkey'
    ad_hu_sele.obs['org'] = 'human'
    ad_mk_sele.obs = ad_mk_sele.obs.loc[:,['Condition','label_subtype','org']]
    ad = ad_hu_sele.concatenate(ad_mk_sele)
    sc.tl.pca(ad, use_highly_variable = False)
    sc.pp.neighbors(ad)
    sc.tl.umap(ad, maxiter=200, min_dist=1, spread=2)
    sc.pl.umap(ad, color=['org'])
    ho = hm.run_harmony(ad.obsm['X_pca'], ad.obs, ['org'])
    ad.obsm['X_harmonypca'] = ho.Z_corr.T
    sc.pp.neighbors(ad, use_rep='X_harmonypca')
    ad.obsm['X_umapraw'] = ad.obsm['X_umap']
    sc.tl.umap(ad)
    ad.obsm['X_umapharmony'] = ad.obsm['X_umap']
    sc.pl.umap(ad, color=['org'])
    ad_homo_temp=ad[ad.obs['org']=='human',]
    ad_mk_temp=ad[ad.obs['org']=='monkey',]
    y_prob, y_pred, y_test, clf, cvsm, acc = SCCAF.SCCAF_assessment(ad_homo_temp.obsm['X_harmonypca'], ad_homo_temp.obs['label_subclass'],n=200)
    aucs = SCCAF.plot_roc(y_prob, y_test, clf, cvsm=cvsm, acc=acc)
    plt.show()
    ad_mk_temp.obs['kamath'] = clf.predict(ad_mk_temp.obsm['X_harmonypca'])
    ad_mk_temp.obs_names= ad_mk_temp.obs_names.str.rsplit('-',1).str[0]
    return ad_mk_temp.obs.loc[:,['label_subtype','kamath']]


if __name__ == '__main__':
    ad_mk = sc.read('../data/DaNs_umap.h5')
    ad_hu = sc.read('../data/homo_DaNs_umap.h5')
    ad_hu = ad_hu[:,ad_hu.var_names.isin(ad_mk.var_names)]
    ad_mk = ad_mk[:,ad_hu.var_names]
    sccaf_predict_summary = pd.DataFrame(index=ad_mk.obs_names)
    for i in range(200,4000,200):
        gene_inter_intersection,gene_intergration = sele_hvg(ad_mk=ad_mk,ad_hu=ad_hu,n_hvg=i)
        if len(gene_inter_intersection)==0:
            print('current intersection_gene = 0 goto intergration')
        else:
            predict_table_intersection = sccaf_predict(ad_mk[:,gene_inter_intersection],ad_hu[:,gene_inter_intersection])
            sccaf_predict_summary[f'hvg_{i}_intersection_length{len(gene_inter_intersection)}']=predict_table_intersection['kamath']
        predict_table_intergration = sccaf_predict(ad_mk[:,gene_intergration],ad_hu[:,gene_intergration])
        sccaf_predict_summary[f'hvg_{i}_intergration_length{len(gene_intergration)}'] = predict_table_intergration['kamath']
    sccaf_predict_summary.to_csv('../data/sccaf_first_predict_result.csv')
