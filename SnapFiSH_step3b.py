import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

ann = pd.read_csv('063022_input_ann.txt',sep='\t')
x = pd.read_csv('063022_output_loop_candidate_CAST.txt',sep='\t')
NUM = ann.shape[0]

x['NegLog10FDR'] = -np.log10(x['FDR'])

u = np.zeros((NUM,NUM))
for i in range(0,x.shape[0]):
    u[ x['bin1'][i], x['bin2'][i] ] = x['NegLog10FDR'][i]
    u[ x['bin2'][i], x['bin1'][i] ] = x['NegLog10FDR'][i]

u[u == 0] = 'nan'


plt.imshow(u, cmap='Reds')
plt.colorbar()
plt.show()


