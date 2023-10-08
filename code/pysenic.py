import os, glob, re, pickle
from functools import partial
from collections import OrderedDict
from cytoolz import compose
import operator as op

import pandas as pd
import seaborn as sns
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt

from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
from pyscenic.aucell import aucell

def derive_regulons(motifs, db_names=('hg19-500bp-upstream-10species.mc9nr','hg19-500bp-upstream-7species.mc9nr',
 'hg19-tss-centered-10kb-10species.mc9nr','hg19-tss-centered-10kb-7species.mc9nr','hg19-tss-centered-5kb-10species.mc9nr','hg19-tss-centered-5kb-7species.mc9nr')):
    motifs.columns = motifs.columns.droplevel(0)

    def contains(*elems):
        def f(context):
            return any(elem in context for elem in elems)
        return f

    # For the creation of regulons we only keep the 10-species databases and the activating modules. We also remove the
    # enriched motifs for the modules that were created using the method 'weight>50.0%' (because these modules are not part
    # of the default settings of modules_from_adjacencies anymore.
    motifs = motifs[
        np.fromiter(map(compose(op.not_, contains('weight>50.0%')), motifs.Context), dtype=np.bool) & \
        np.fromiter(map(contains(*db_names), motifs.Context), dtype=np.bool) & \
        np.fromiter(map(contains('activating'), motifs.Context), dtype=np.bool)]
    # We build regulons only using enriched motifs with a NES of 3.0 or higher; we take only directly annotated TFs or TF annotated
    # for an orthologous gene into account; and we only keep regulons with at least 10 genes.
    regulons = list(filter(lambda r: len(r) >= 10, df2regulons(motifs[(motifs['NES'] >= 3.0) 
                                                                      & ((motifs['Annotation'] == 'gene is directly annotated')
                                                                        | (motifs['Annotation'].str.startswith('gene is orthologous to')
                                                                           & motifs['Annotation'].str.endswith('which is directly annotated for motif')))
                                                                     ])))
    
    # Rename regulons, i.e. remove suffix.
    return list(map(lambda r: r.rename(r.transcription_factor), regulons))

os.system(f'pyscenic grn DaNs_exp_CPM.csv /media/tang/raid2/Monkey_Parkinson/regulon_GSE157783/hs_hgnc_curated_tfs.txt\
                            -o DaNs_exp_CPM.adjacencies.tsv --num_workers 48')

file_path = '/media/tang/raid2/Monkey_Parkinson/regulon_GSE157783/'
os.system(f'pyscenic ctx DaNs_exp_CPM.adjacencies.tsv \
        {file_path}hg19-500bp-upstream-10species.mc9nr.feather {file_path}hg19-tss-centered-5kb-10species.mc9nr.feather \
        {file_path}hg19-tss-centered-10kb-10species.mc9nr.feather --annotations_fname {file_path}motifs-v9-nr.hgnc-m0.001-o0.0.tbl\
        --expression_mtx_fname DaNs_exp_CPM.csv --output DaNs_exp_CPM.motifs.csv --num_workers 48')
df_motifs = load_motifs('./DaNs_exp_CPM.motifs.csv')
regulons = derive_regulons(df_motifs)
exp_mtx = pd.read_csv('./DaNs_exp_CPM.csv',index_col=0)
auc_mtx = aucell(exp_mtx, regulons, num_workers=32)
ad = sc.read('./SN_clustered.h5')
ad_dans = ad[ad.obs['label']=='DaNs',]
add_scenic_metadata(ad_dans, auc_mtx, regulons)

