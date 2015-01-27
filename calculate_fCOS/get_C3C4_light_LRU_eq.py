import sys, os, os.path
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

C4_plants = ['corn-light', 'sorghum-light', 'amaranthus-light']
C3_plants = ['sage-light', 'hibiscus-light']
#light_experiments = C3_plants.append(C4_plants)

fpath = os.path.join(os.getenv('HOME'), 'work', 'Data',
                     'Stimler_COS_exchange_data.csv')
all_data = pd.read_csv(fpath, sep=',', skiprows=[1])
light_data = all_data.loc[:, ['PAR_umol m-2s-1', 'lru', 'plant']]
light_data = light_data.ix[all_data.plant.str.contains('(light)', na=False)]

# light_data.loc[:, ('C3C4')] = 'C3'
isC4 = light_data.plant.str.contains('(corn|sorghum|amaranthus)')
light_data.loc[isC4, ('C3C4')] = 'C4'
light_data.loc[~isC4, ('C3C4')] = 'C3'

light_data.columns = light_data.columns = [u'PAR', u'LRU', u'plant', u'C3C4']

sns.set_style('ticks')
sns.lmplot("PAR", "LRU", light_data, hue="C3C4", markers=["x", "o"], lowess=True)
g = sns.lmplot("PAR", "LRU", light_data, hue="plant", lowess=True, legend_out=False)
plt.show()
# g = sns.FacetGrid(light_data, hue="plant", size=5, 
#                   palette=sns.color_palette('Dark2'),
#                   legend_out=False)
# g.map(plt.scatter, "PAR", "LRU")
# g.add_legend();
# plt.show()
