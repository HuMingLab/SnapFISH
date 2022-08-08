import numpy as np
import pandas as pd

def euc_dist(u,v):
    d = np.sqrt((u[0]-v[0])**2 + (u[1]-v[1])**2 + (u[2]-v[2])**2 )
    return round(d, 4)
    


a = pd.read_csv('063022_input_3D_coordinates_CAST.txt',sep='\t')
ann = pd.read_csv('063022_input_ann.txt',sep='\t')
NUM = ann.shape[0] # number of loci
BINSIZE = ann['end'][0] - ann['start'][0]

NUM_BINPAIR = NUM*(NUM-1)/2
rec = pd.DataFrame(np.zeros((int(NUM_BINPAIR), 2)),columns = ['bin1', 'bin2'])
k = 0

for i in range(1,NUM):
    for j in range((i+1),(NUM+1)):
        rec['bin1'][k] = i
        rec['bin2'][k] = j
        k += 1


NAME = sorted(list(set(a['cell_id'])))
colname = ['cell.'+ str(x) for x in NAME]

out=pd.DataFrame(np.zeros((int(NUM_BINPAIR), len(NAME))), columns = colname)

for ID in range(0,len(NAME)):
    x = a[a['cell_id'] == NAME[ID]].reset_index(drop=True)
    for i in range(0,(x.shape[0]-1)):
        for j in range(i+1,x.shape[0]):
            
            posi = x['pos'][i]
            posj = x['pos'][j]
            d = euc_dist(x.iloc[i, 2:5].to_numpy(), x.iloc[j, 2:5].to_numpy())
            out.loc[(rec['bin1'] == posi) & (rec['bin2'] == posj) ,ID] = d

out1=out.replace(0, np.nan)
final= pd.concat([rec.reset_index(drop=True), out1], axis=1)
final.to_csv('output_3D_dist_CAST.txt', header=True, index=None, sep='\t')

out_mean=final[final.columns.drop(['bin1','bin2'])].mean(axis=1)
out_mean=out_mean.to_frame('out.mean')
output=pd.concat([rec.reset_index(drop=True), out_mean], axis=1)

output.to_csv('output_3D_dist_CAST_avg.txt', header=True, index=None, sep='\t')
