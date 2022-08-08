import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ann = pd.read_csv('063022_input_ann.txt',sep='\t')
a = pd.read_csv('063022_output_3D_dist_129_avg.txt',sep='\t')
NUM = ann.shape[0] # number of loci

x=np.zeros((NUM,NUM))
a[['bin1','bin2']] -= 1

for i in range(0,a.shape[0]):
    x[ a['bin1'][i], a['bin2'][i] ] = a['out.mean'][i]
    x[ a['bin2'][i], a['bin1'][i] ] = a['out.mean'][i]

np.fill_diagonal(x, np.nan)

plt.imshow(x, cmap='RdBu')
plt.colorbar()
plt.show()


