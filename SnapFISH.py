import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.stats.multitest as multi
import math
import argparse

import sys
import os


    
    
def main():
    parser = create_parser()
    args = parser.parse_args()
    
    SnapFISH_dir = args.indirSF
    dir_in = args.indir
    ann_dir = args.anndir
    data_name = args.dataname
    
    out_dir = args.outdir
    
    save_pic = args.pic
    
    ann = pd.read_csv(ann_dir, header='infer', sep="\t", lineterminator='\n')
    
    output_avg_dist, final_3d_dist = step1a(dir_in,ann,out_dir,data_name)
    
    if (save_pic == 0):
        step1b(ann,output_avg_dist,data_name,out_dir)
    
    output_Ttest = step2(ann,final_3d_dist,data_name,out_dir)
    loop_candidate = step3a(output_Ttest,ann,data_name,out_dir)
    if (save_pic == 0):
        step3b(ann,loop_candidate,data_name,out_dir)
        
    loop_summit = step4a(loop_candidate,data_name,out_dir)
    
    if (save_pic == 0):
        step4b(ann,loop_summit ,out_dir,data_name)
        
        
    print("Exiting program", flush=True)
   
   
def euc_dist(u,v):
    d = np.sqrt((u[0]-v[0])**2 + (u[1]-v[1])**2 + (u[2]-v[2])**2 )
    return round(d, 4)

def step1a(dir_in,ann,out_dir,data_name):
    a = pd.read_csv(dir_in, header='infer', sep="\t")
    
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
    final_3d_dist= pd.concat([rec.reset_index(drop=True), out1], axis=1)
    #final.to_csv('output_3D_dist_129.txt', header=True, index=None, sep='\t')

    out_mean=final_3d_dist[final_3d_dist.columns.drop(['bin1','bin2'])].mean(axis=1)
    out_mean=out_mean.to_frame('out.mean')
    output_avg_dist=pd.concat([rec.reset_index(drop=True), out_mean], axis=1)

   # output.to_csv('output_3D_dist_129_avg.txt', header=True, index=None, sep='\t')
    newfilename1 = os.path.join(out_dir, "tempfile/output_3D_dist_%s.txt" % (data_name))
    newfilename2 = os.path.join(out_dir, "tempfile/output_3D_dist_avg_%s.txt" % (data_name))
    final_3d_dist.to_csv(newfilename1,header=True, index=None, sep='\t')
    output_avg_dist.to_csv(newfilename2, header=True, index=None, sep='\t')
    return output_avg_dist, final_3d_dist

def step1b(ann,a,data_name,out_dir):
    NUM = ann.shape[0] # number of loci

    x=np.zeros((NUM,NUM))
    a[['bin1','bin2']] -= 1

    for i in range(0,a.shape[0]):
        x[ int(a['bin1'][i]), int(a['bin2'][i]) ] = a['out.mean'][i]
        x[ int(a['bin2'][i]), int(a['bin1'][i]) ] = a['out.mean'][i]

    np.fill_diagonal(x, np.nan)

    plt.imshow(x, cmap='RdBu')
    plt.colorbar()
    newfilename1 = os.path.join(out_dir, "tempfile/heatmap_av_Euc_dist_%s.pdf" % (data_name))
    plt.savefig(newfilename1)

def step2(ann,a,data_name,out_dir):
    NUM = ann.shape[0] # number of loci
    BINSIZE = ann['end'][0] - ann['start'][0]

    NUMCELL =a.shape[1]-2

    CUT_UP = 50e3 # 50Kb upper bound for local background model
    CUT_LO = 20e3

    rec = pd.DataFrame(np.zeros((a.shape[0],8)),columns = ['Case', 'Ctrl', 'Tstat', 'Pvalue', 'Ctrl.ll', 'Ctrl.h', 'Ctrl.v','fdr'])


    v1=a[a.columns.drop(['bin1','bin2'])]

    for i in range(0,a.shape[0]):
        v = v1.iloc[i,:]
        v = v.dropna()

        x = a[ ((abs(a['bin1']-a['bin1'][i]) <= CUT_UP/BINSIZE) & (abs(a['bin2']-a['bin2'][i])<= CUT_UP/BINSIZE))
              & ~((abs(a['bin1']-a['bin1'][i]) <= CUT_LO/BINSIZE) & (abs(a['bin2']-a['bin2'][i])<= CUT_LO/BINSIZE))]
        
        # local background
        xll = x[ (x['bin1'] > a['bin1'][i]) & (x['bin2'] < a['bin2'][i])] # lower left
        xh  = x[ x['bin1'] == a['bin1'][i] ] # horizontal
        xv  = x[ x['bin2'] == a['bin2'][i] ] # vertical

        x_mean   = np.empty(NUMCELL)
        xll_mean = np.empty(NUMCELL)
        xh_mean  = np.empty(NUMCELL)
        xv_mean  = np.empty(NUMCELL)

        for j in range(0,NUMCELL): # calculate average distance in each cell for different local background models
        
            x_mean[j]   = x.iloc[:,j+2].mean(skipna = True)
            xll_mean[j] = xll.iloc[:,j+2].mean(skipna = True)
            xh_mean[j]  = xh.iloc[:,j+2].mean(skipna = True)
            xv_mean[j]  = xv.iloc[:,j+2].mean(skipna = True)
        
        x_mean   = x_mean[~np.isnan(x_mean)]
        xll_mean   = xll_mean[~np.isnan(xll_mean)]
        xh_mean   = xh_mean[~np.isnan(xh_mean)]
        xv_mean   = xv_mean[~np.isnan(xv_mean)]
        
        t_stat, p_val=stats.ttest_ind(v, x_mean, equal_var=False)

        rec['Case'][i] = np.mean(v)
        rec['Ctrl'][i] = np.mean(x_mean)
        rec['Tstat'][i] = t_stat
        rec['Pvalue'][i] = p_val
        rec['Ctrl.ll'][i] = np.mean(xll_mean)
        rec['Ctrl.h'][i] = np.mean(xh_mean)
        rec['Ctrl.v'][i] = np.mean(xv_mean)
        
        
    rej, pval_corr = multi.multipletests(rec['Pvalue'], method='fdr_bh')[:2]
    rec['fdr'] = pval_corr
    output=pd.concat([a[['bin1','bin2']], rec], axis=1)
    newfilename1 = os.path.join(out_dir, "tempfile/output_Ttest_%s.txt" % (data_name))
    output.to_csv(newfilename1,header=True, index=None, sep='\t')
    #output.to_csv('output_Ttest_CAST.txt', header=True, index=None, sep='\t')
    
    return output
    
def step3a(a,ann,data_name,out_dir):
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

    x = a[ (a['Tstat'] < -4) & (a['fdr'] < 0.1) & ((a['Ctrl']/a['Case']) > cutoff1)
        & ((a['Ctrl.ll']/a['Case']) > cutoff2) & ((a['Ctrl.h']/a['Case']) > cutoff2) & ((a['Ctrl.v']/a['Case']) > cutoff2)]


    x = x.iloc[np.argsort(x['fdr']),]
    x = x.reset_index(drop=True)

    rec = pd.DataFrame(np.zeros((x.shape[0],6)),columns = ['chr1', 'x1', 'x2','chr2', 'y1', 'y2'])

    for i in range(0,rec.shape[0]):
        rec.iloc[i,0:3] = ann.iloc[ np.where(ann['pos'] == x['bin1'][i])[0],0:3].to_numpy()
        rec.iloc[i,3:6] = ann.iloc[ np.where(ann['pos'] == x['bin2'][i])[0] , 0:3].to_numpy()

    out =pd.concat([rec.reset_index(drop=True), x], axis=1)
    newfilename1 = os.path.join(out_dir, "tempfile/output_loop_candidate_%s.txt" % (data_name))
    out.to_csv(newfilename1,header=True, index=None, sep='\t')
    #out.to_csv('output_loop_candidate_129.txt', header=True, index=None, sep='\t')
    
    return out

def step3b(ann,x,data_name,out_dir):
    
    NUM = ann.shape[0]

    x['NegLog10FDR'] = -np.log10(x['fdr'])

    u = np.zeros((NUM,NUM))
    for i in range(0,x.shape[0]):
        u[ int(x['bin1'][i]), int(x['bin2'][i]) ] = x['NegLog10FDR'][i]
        u[ int(x['bin2'][i]), int(x['bin1'][i]) ] = x['NegLog10FDR'][i]

    u[u == 0] = 'nan'


    plt.imshow(u, cmap='Reds')
    plt.colorbar()
    newfilename1 = os.path.join(out_dir, "tempfile/heatmap_loop_candidate_%s.pdf" % (data_name))
    plt.savefig(newfilename1)

def step4a(y0,data_name,out_dir):
    GAP = 10e3 # 10Kb gap for  5Kb bin resolution

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
                v['NegLog10FDR'][ v['label'] == i ] = sum( -np.log10(vtmp['fdr']) )
           
            for i in range(0,v.shape[0]):
           
                vtmp = v[ v['label'] == v['label'][i]]
                
                if(v['Tstat'][i] == min(vtmp['Tstat'])):
              
                   v.loc[i,'summit'] = 1
              
           

            out = pd.concat([out.reset_index(drop=True), v])
       
        final = pd.concat([final.reset_index(drop=True), out])

    out = final[ final['summit']== 1 ]
    newfilename1 = os.path.join(out_dir, "output_loop_summit_%s.txt" % (data_name))
    out.to_csv(newfilename1,header=True, index=None, sep='\t')
    #out.to_csv('output_loop_summit_CAST.txt', header=True, index=None, sep='\t')

    return out
    
def step4b(ann,x,out_dir,data_name):
    NUM = ann.shape[0]

    u = np.zeros((NUM,NUM))
    for i in range(0,x.shape[0]):
        u[ int(x['bin1'][i]), int(x['bin2'][i]) ] = 1
        u[ int(x['bin2'][i]), int(x['bin1'][i]) ] = 1
        

    plt.imshow(u, cmap='Reds')
    plt.colorbar()
    newfilename1 = os.path.join(out_dir, "tempfile/heatmap_loop_summits_%s.pdf" % (data_name))
    plt.savefig(newfilename1)

def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--indirSF', action = 'store', required = True, \
                        help = 'SnapFISH directory')
    parser.add_argument('-i', '--indir', action = 'store', required = True, \
                        help = 'input directory')
    parser.add_argument('-o', '--outdir', action = 'store', \
                        required = True, help = 'output directory')
    parser.add_argument('-a', '--anndir', action = 'store', \
                        required = True, help = 'file list directory')
    parser.add_argument('-p', '--pic', type = int, help = '0 to save figure, and 1 to not', \
                        required = True)
    parser.add_argument('-d', '--dataname', action = 'store', help = 'data name', \
                        required = True)
    return parser
    


    
if __name__ == "__main__":
    main()
