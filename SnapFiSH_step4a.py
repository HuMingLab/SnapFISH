import numpy as np
import pandas as pd

GAP = 10e3 # 10Kb gap for  5Kb bin resolution

y0 = pd.read_csv('063022_output_loop_candidate_CAST.txt',sep='\t')

chr_name=sorted(list(set(y0['chr1'])))

final = pd.DataFrame()

for CHR in chr_name:
    y = y0[ y0['chr1']==CHR]

# step 1: find neighborhood
    for i in range(0,y.shape[0]):
        z = y[ (abs(y['x1'] - y['x1'][i])<= GAP) & (abs(y['y1'] - y['y1'][i]) <= GAP)]
        y.loc[i,'CountNei'] = z.shape[0]


# step 2: split singletons and peak clusters with >=3 bin pairs
    v = y[ y['CountNei'] >= 2] # peak cluster: sharp peak + broad peak

# step 3: keep singletons
    out = pd.DataFrame()

# step 4: for clusters, find summit
    if v.shape[0] > 0:
   # for peak cluster, assign label 1, 2, ..., N
        v['label'] = np.array(range(1,len(v)+1))

   # for all bin pairs within the neighborhood
   # assign the same cluster label, using the minimal label
        for i in range(0,len(v)):
      
            w = v[ (abs(v['x1'] - v['x1'][i])<= GAP) & (abs(v['y1'] - v['y1'][i]) <= GAP)]
            w_min = min(w['label'])
            w_label = sorted(list(set(w['label'])))
            for j in range(0,len(w_label)):
          
                v['label'][ v['label'] == w_label[j] ] = w_min
          
       

    # assign consecutive label number
        v_rec = pd.DataFrame(sorted(list(set( v['label'] ) )))
        v_rec =pd.concat([v_rec.reset_index(drop=True), pd.DataFrame(np.array(range(1,len(v_rec)+1)))], axis=1)
        for i in range(0,v.shape[0]):
         
            v.loc[i,'label'] = np.asarray(v_rec[ v_rec.iloc[:,1] == v['label'][i]])[0,1]
         

    # count cluster size, find cluster summit
        v['ClusterSize'] = 0
        v['summit'] = 0
        v['NegLog10FDR'] = 0
        for i in range(1,v_rec.shape[0]+1):
       
            vtmp = v[ v['label'] == i ]
            v['ClusterSize'][ v['label'] == i ] = vtmp.shape[0]
            v['NegLog10FDR'][ v['label'] == i ] = sum( -np.log10(vtmp['FDR']) )
       
        for i in range(0,v.shape[0]):
       
            vtmp = v[ v['label'] == v['label'][i]]
            
            if(v['Tstat'][i] == min(vtmp['Tstat'])):
          
               v.loc[i,'summit'] = 1
          
       

        out = pd.concat([out.reset_index(drop=True), v])
   
    final = pd.concat([final.reset_index(drop=True), out])

out = final[ final['summit']== 1 ]
out.to_csv('output_loop_summit_CAST.txt', header=True, index=None, sep='\t')
