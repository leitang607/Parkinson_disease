import pandas as pd
import numpy as np
import scanpy as sc
from matplotlib import *
import matplotlib.pyplot as plt
rcParams['pdf.fonttype']=42
rcParams['ps.fonttype']=42

#mouse striatum
ad_striatum = sc.read('/media/tang/raid2/Monkey_Parkinson/GSE116470/Striatum_files/striatum_no_doublets_labeled.h5')
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
ad_snc = sc.read('/media/tang/raid2/Monkey_Parkinson/GSE116470/SNC_files/snc_no_doublets_labeled.h5')
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
ad_snc = sc.read('./Figures_pip_ver2/snc_figout.h5')
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
sc.pl.violin(ad_snc,groupby='label_5',keys=['FOXP2'],rotation=90,stripplot=False,palette=color_order,ax=ax,linewidth=0.3)

#monkey PT
ad_put = sc.read('Figures_pip_ver2/putamen_figout.h5')
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
sc.pl.violin(ad_put,groupby='label_5',keys=['FOXP2'],rotation=90,stripplot=False,palette=color_order,ax=ax,linewidth=0.3)

#human SN
ad_snc = sc.read('/media/tang/raid2/Monkey_Parkinson/GSE157783/gse157783_umap.h5')
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
ad_snc = sc.read('/media/tang/raid2/Monkey_Parkinson/GSE152058/GSE152058_human.h5')
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
sc.pl.violin(ad_snc,groupby='CellType',keys=['FOXP2'],rotation=90,stripplot=False,palette=color_order,ax=ax,linewidth=0.3)
