import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['ps.fonttype']=42
import seaborn as sns

plot_data = pd.read_csv('../data/Kurlan_scores.csv')
color_dict = {'81005':'#ed7d31','102006':'#5b9bd5'}
plot_data['ID'] = plot_data['ID'].astype(str)
fig,ax = plt.subplots(figsize=(3,2))
sns.lineplot(plot_data,x='Day',y='Total',hue='ID',style='ID',markers=True,palette=color_dict,ax=ax,markersize=6,linewidth=1,markeredgewidth=0)
ax.set_ylabel('Kurlan Scores')
ax.set_xlabel('Days')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
fig.savefig('./Figures/Kurlan_scores.pdf')
