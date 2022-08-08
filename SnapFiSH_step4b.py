import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ann = pd.read_csv('063022_input_ann.txt',sep='\t')
x = pd.read_csv('063022_output_loop_summit_129.txt',sep='\t')
NUM = ann.shape[0]

u = np.zeros((NUM,NUM))
for i in range(0,x.shape[0]):
    u[ x['bin1'][i], x['bin2'][i] ] = 1
    u[ x['bin2'][i], x['bin1'][i] ] = 1
    

plt.imshow(u, cmap='Reds')
plt.colorbar()
plt.show()



