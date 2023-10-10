import pandas as pd
import numpy as np
import scanpy as sc
from matplotlib import *
import matplotlib.pyplot as plt
rcParams['pdf.fonttype']=42
rcParams['ps.fonttype']=42

#mouse striatum
ad_striatum = sc.read('../data/GSE116470_Striatum.h5')
color_order = ['#0F90BF' for i in range(0,10)]
fig,ax = plt.subplots(figsize=(5,0.75),dpi=300)
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.size']=7
plt.rcParams['font.weight']='normal'
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.set_ylim(0,6)
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax.xaxis.set_minor_formatter(plt.NullFormatter())
sc.pl.violin(ad_striatum,groupby='label_main',keys=['Foxp2'],rotation=90,stripplot=False,palette=color_order,ax=ax,linewidth=0.3)

#mouse SN
ad_snc = sc.read('../data/GSE116470_SNC.h5')
color_order = ['#0F90BF' for i in range(0,10)]
fig,ax = plt.subplots(figsize=(5,0.75),dpi=300)
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.size']=7
plt.rcParams['font.weight']='normal'
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.set_ylim(0,6)
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax.xaxis.set_minor_formatter(plt.NullFormatter())
sc.pl.violin(ad_snc,groupby='label_main',keys=['Foxp2'],rotation=90,stripplot=False,palette=color_order,ax=ax,linewidth=0.3)

#monkey SN
ad_snc = sc.read('../data/SN_umap.h5')
color_order = ['#E8320E' for i in range(0,10)]
fig,ax = plt.subplots(figsize=(5,0.75),dpi=300)
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.size']=7
plt.rcParams['font.weight']='normal'
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.set_ylim(0,4)
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax.xaxis.set_minor_formatter(plt.NullFormatter())
sc.pl.violin(ad_snc,groupby='label',keys=['FOXP2'],rotation=90,stripplot=False,palette=color_order,ax=ax,linewidth=0.3)

#monkey PT
ad_put = sc.read('../data/putamen_umap.h5')
color_order = ['#E8320E' for i in range(0,10)]
fig,ax = plt.subplots(figsize=(5,0.75),dpi=300)
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.size']=7
plt.rcParams['font.weight']='normal'
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.set_ylim(0,4)
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax.xaxis.set_minor_formatter(plt.NullFormatter())
sc.pl.violin(ad_put,groupby='label',keys=['FOXP2'],rotation=90,stripplot=False,palette=color_order,ax=ax,linewidth=0.3)

#human SN
ad_snc = sc.read('../data/GSE157783_human_SN_umap.h5')
color_order = ['#F7AF00' for i in range(0,10)]
fig,ax = plt.subplots(figsize=(5,0.75),dpi=300)
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.size']=7
plt.rcParams['font.weight']='normal'
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.set_ylim(0,6)
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax.xaxis.set_minor_formatter(plt.NullFormatter())
sc.pl.violin(ad_snc,groupby='cell_ontology',keys=['FOXP2'],rotation=90,stripplot=False,palette=color_order,ax=ax,linewidth=0.3)

#human PT
ad_pt = sc.read('../data/GSE152058_human_PT_umap.h5')
color_order = ['#F7AF00' for i in range(0,10)]
fig,ax = plt.subplots(figsize=(5,0.75),dpi=300)
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.size']=7
plt.rcParams['font.weight']='normal'
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.set_ylim(0,6)
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax.xaxis.set_minor_formatter(plt.NullFormatter())
sc.pl.violin(ad_pt,groupby='CellType',keys=['FOXP2'],rotation=90,stripplot=False,palette=color_order,ax=ax,linewidth=0.3)
