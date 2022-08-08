import numpy as np
import pandas as pd

a = pd.read_csv('063022_output_Ttest_129.txt',sep='\t')
ann = pd.read_csv('063022_input_ann.txt',sep='\t')
NUM = ann.shape[0] # number of loci
BINSIZE = ann['end'][0] - ann['start'][0]

# only consider bin pairs 100Kb ~ 1MB
BOUND_LO = 1e5
BOUND_UP = 1e6

# user defined cutoff
cutoff1 = 1.1
cutoff2 = 1.05

a = a[ (a['bin2']-a['bin1']>=BOUND_LO/BINSIZE) & (a['bin2']-a['bin1']<=BOUND_UP/BINSIZE)
    & (~np.isnan(a['Ctrl'])) & (~np.isnan(a['Ctrl.ll'])) & (~np.isnan(a['Ctrl.h'])) & (~np.isnan(a['Ctrl.v']))]

x = a[ (a['Tstat'] < -4) & (a['FDR'] < 0.1) & ((a['Ctrl']/a['Case']) > cutoff1)
    & ((a['Ctrl.ll']/a['Case']) > cutoff2) & ((a['Ctrl.h']/a['Case']) > cutoff2) & ((a['Ctrl.v']/a['Case']) > cutoff2)]


x = x.iloc[np.argsort(x['FDR']),]
x = x.reset_index(drop=True)

rec = pd.DataFrame(np.zeros((x.shape[0],6)),columns = ['chr1', 'x1', 'x2','chr2', 'y1', 'y2'])

for i in range(0,rec.shape[0]):
    rec.iloc[i,0:3] = ann.iloc[ np.where(ann['pos'] == x['bin1'][i])[0],0:3].to_numpy()
    rec.iloc[i,3:6] = ann.iloc[ np.where(ann['pos'] == x['bin2'][i])[0] , 0:3].to_numpy()

out =pd.concat([rec.reset_index(drop=True), x], axis=1)
out.to_csv('output_loop_candidate_129.txt', header=True, index=None, sep='\t')

