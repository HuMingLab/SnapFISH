import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
from natsort import natsorted
from scipy import stats
import statsmodels.stats.multitest as multi
import warnings
import glob
import concurrent.futures
from functools import partial

import argparse
import sys
import os

def main():
    parser = create_parser()
    args = parser.parse_args()
    
    SnapFISH_dir = args.indirSF
    
    ann_dir = args.anndir
    data_name = args.dataname
    
    out_dir = args.outdir
    dir_in_r1 = args.rep1
    dir_in_r2 = args.rep2
    max_worker = args.num_cpus
    
    ann = pd.read_csv(ann_dir, header='infer', sep="\t", lineterminator='\n')
    chrom=list(range(1, 20))
    with concurrent.futures.ProcessPoolExecutor(max_worker) as executor:
        for result in executor.map(partial(run_step123, dir_in_r1=dir_in_r1, dir_in_r2=dir_in_r2, \
                                              ann_dir = ann, out_dir=out_dir),chrom):
            print("..", flush=True)
        
    step4_result = SnapFISH25_step4(out_dir,ann)
    step5_result = SnapFISH25_step5(step4_result)
    SnapFISH25_step6(step5_result,out_dir)


    
def euc_dist(u,v):
    d = np.sqrt((u[0]-v[0])**2 + (u[1]-v[1])**2 + (u[2]-v[2])**2 )
    return round(d, 4)
  
def SnapFISH25_step1(dir_in_r1,dir_in_r2,chrom,out_dir):
    CHR = chrom
    for ii in range(1,3):
        if (ii == 1):
            dir_in = dir_in_r1
        else:
            dir_in = dir_in_r2
        a0 = pd.read_csv(dir_in, header='infer', sep="\t")
        a = a0[ a0['chromID']==CHR]

        CELL = natsorted(list(set(a['alleleNAME'])))
        NUM = len(CELL)
        rec =pd.DataFrame(np.zeros((1770, 2)), columns = ['bin1', 'bin2'])
        k=1
        
        for i in range(1,60):
            for j in range((i+1),61):
                rec['bin1'][k-1] = i
                rec['bin2'][k-1] = j
                k += 1
                
        out=pd.DataFrame(np.zeros((len(rec), len(CELL))), columns = CELL)

        for ID in range(0,len(CELL)):
            x = a[a['alleleNAME'] == CELL[ID]].reset_index(drop=True)
            for i in range(0,(x.shape[0]-1)):
                for j in range(i+1,x.shape[0]):
                    
                    hybi = x['regionID..hyb1.60.'][i]
                    hybj = x['regionID..hyb1.60.'][j]
                    d = euc_dist(x.iloc[i, 4:7].to_numpy(), x.iloc[j, 4:7].to_numpy())
                    out.iloc[(rec['bin1'] == hybi) & (rec['bin2'] == hybj) ,ID] = d


        out1=out.replace(0, np.nan)
        final= pd.concat([rec.reset_index(drop=True), out1], axis=1)
        newfilename1 = os.path.join(out_dir, "tempfile/Rep%s_PairwiseDist_chr%s.txt" % (ii,CHR))
        final.to_csv(newfilename1,header=True, index=None, sep='\t')
        
    
def SnapFISH25_step2(chrom, out_dir):
    CHR = chrom
    newfilename1 = os.path.join(out_dir, "tempfile/Rep%s_PairwiseDist_chr%s.txt" % (1,CHR))
    newfilename2 = os.path.join(out_dir, "tempfile/Rep%s_PairwiseDist_chr%s.txt" % (2,CHR))
    x = pd.read_csv(newfilename1,sep='\t')
    y = pd.read_csv(newfilename2,sep='\t')

    rec = x.iloc[:,0:2]
    rec['d'] = np.nan
    xx=x[x.columns.drop(['bin1','bin2'])]
    yy=y[y.columns.drop(['bin1','bin2'])]

    for i in range(0,len(x)):
        u = xx.loc[i,:]
        v = yy.loc[i,:]
        w = np.concatenate((u, v))
        w = w[~np.isnan(w)]
        rec['d'][i] = np.mean(w)
        
    newfilename3 = os.path.join(out_dir, "tempfile/AvgDist_chr%s.txt" % (CHR))
    rec.to_csv(newfilename3,header=True, index=None, sep='\t')
    return x,y

def SnapFISH25_step3(chrom, ann1, x, y,out_dir):
    CHR = chrom

    a0 = ann1
    ann = a0[a0['Chrom ID'] == CHR].reset_index(drop=True)
    ann['EndID'] = (ann['End'] - ann['Start'][0])/25000
      
    a1 = x
    a2 = y


    a1_col = list(a1.columns)[2:]
    a2_col = list(a2.columns)[2:]
    a1_col_rep = ["Rep1." + word for word in a1_col]
    a2_col_rep = ["Rep2." + word for word in a2_col]
    
    a = pd.DataFrame(columns=['bin1','bin2','bin1.EndID','bin2.EndID'])
    
    a[['bin1','bin2']]=a1.iloc[:,0:2]
    a['bin1.EndID']=0
    a['bin2.EndID']=0
    a[a1_col_rep]=a1.iloc[:,2:]
    a[a2_col_rep]=a2.iloc[:,2:]

    for i in range(0,len(a)):
        a['bin1.EndID'][i] = ann['EndID'][ int(a['bin1'][i])-1 ]
        a['bin2.EndID'][i] = ann['EndID'][ int(a['bin2'][i])-1 ]
        
        
    rec = pd.DataFrame(np.nan, index = np.arange(len(a)),
                       columns = ['Case', 'Ctrl', 'Tstat', 'Pvalue', 'Ctrl.ll', 'Ctrl.h', 'Ctrl.v'])

      
    for i in range(0,len(rec)):
        v = a.iloc[i, 4:]
        v = v[~np.isnan(v)]
        
        x = a[   ((abs(a['bin1.EndID']-a['bin1.EndID'][i])<=2) & (abs(a['bin2.EndID']-a['bin2.EndID'][i])<=2))
                   & ~((a['bin1.EndID'] == a['bin1.EndID'][i]) & (a['bin2.EndID'] == a['bin2.EndID'][i]) )]
        xll = x[ (x['bin1.EndID'] > a['bin1.EndID'][i]) & (x['bin2.EndID'] < a['bin2.EndID'][i])]
        xh  = x[ x['bin1.EndID'] == a['bin1.EndID'][i]]
        xv  = x[ x['bin2.EndID'] == a['bin2.EndID'][i]]
        
        
        if len(x)==24: # require 24 bin pairs in the local neighborhood, all without missing data
            x_mean = np.empty(x.shape[1]-4)
            #x_mean.fill(np.nan)
            x_ll   = np.empty(x.shape[1]-4)
            x_h    = np.empty(x.shape[1]-4)
            x_v    = np.empty(x.shape[1]-4)
            
            for j in range(0,len(x_mean)): # average per allele
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    
                    x_mean[j] = np.nanmean(x.iloc[:,j+4])
                    x_ll[j]   = np.nanmean(xll.iloc[:,j+4])
                    x_h[j]    = np.nanmean(xh.iloc[:,j+4])
                    x_v[j]    = np.nanmean(xv.iloc[:,j+4])

            x_mean = x_mean[ ~np.isnan(x_mean) ]
            x_ll   = x_ll[ ~np.isnan(x_ll) ]
            x_h    = x_h[ ~np.isnan(x_h) ]
            x_v    = x_v[ ~np.isnan(x_v) ]
          
            if((len(v)>=2) & (len(x_mean)>=2)):
          
                    
                t_stat, p_val=stats.ttest_ind(v, x_mean, equal_var=False)

                rec['Case'][i] = np.mean(v)
                rec['Ctrl'][i] = np.mean(x_mean)
                rec['Tstat'][i] = t_stat
                rec['Pvalue'][i] = p_val
                rec['Ctrl.ll'][i] = np.mean(x_ll)
                rec['Ctrl.h'][i] = np.mean(x_h)
                rec['Ctrl.v'][i] = np.mean(x_v)

    out = pd.concat([a[['bin1','bin2','bin1.EndID','bin2.EndID']], rec], axis=1)
    out = out[ ~np.isnan(out['Pvalue']) & ~np.isnan(out['Ctrl']) & ~np.isnan(out['Ctrl.ll']) & ~np.isnan(out['Ctrl.h'])
              & ~np.isnan(out['Ctrl.v'])]
    out = out[ (out['bin2.EndID']-out['bin1.EndID']>4) & (out['bin2.EndID']-out['bin1.EndID']<40)] # (100Kb, 1Mb)
     
    rej, pval_corr = multi.multipletests(out['Pvalue'], method='fdr_bh')[:2]

    out['FDR'] = pval_corr
    out['chr'] = CHR
    out['ratio']    = out['Ctrl']/out['Case']
    out['ratio.ll'] = out['Ctrl.ll']/out['Case']
    out['ratio.h']  = out['Ctrl.h']/out['Case']
    out['ratio.v']  = out['Ctrl.v']/out['Case']


    column_to_move = out.pop('chr')
    final = out
    final.insert(0, 'chr', column_to_move)
    newfilename1 = os.path.join(out_dir, "tempfile/Ttest_results_chr%s.txt" % (CHR))
    final.to_csv(newfilename1, header=True, index=None, sep='\t')

     
def SnapFISH25_step4(out_dir,ann):
    newfilename1 = os.path.join(out_dir, "tempfile/Ttest_results_chr*.txt")
    interesting_files = glob.glob(newfilename1)
    y = pd.concat((pd.read_csv(f, header = 0,sep='\t') for f in interesting_files))
    y = y.reset_index(drop=True)
    #y.to_csv('Ttest_results_allchr.txt', header=True, index=None, sep='\t')

    a = ann


    z = pd.DataFrame(np.zeros((len(y), 6)), columns = ['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2'])

    z['chr1'] = y['chr'].reset_index(drop=True)
    z['chr2'] = y['chr'].reset_index(drop=True)

    for i in range(0,len(z)):

        b = a[ a['Chrom ID'] == y['chr'][i] ]
        b=b.reset_index(drop=True)
        z['x1'][i] = b['Start'][ y['bin1'][i]-1 ]
        z['x2'][i] = z['x1'][i]+25e3
        z['y1'][i] = b['Start'][ y['bin2'][i]-1 ]
        z['y2'][i] = z['y1'][i]+25e3


    out = pd.concat([z, y.iloc[:,1:]], axis=1)
    
    return out
    
def SnapFISH25_step5(step4_result):
    yy = step4_result
    yy['d'] = yy['y1'] - yy['x1']


    cutoff1 = 1.1
    cutoff2 = 1.05
    u = yy[ (yy['Tstat'] < -4) & (yy['FDR']<0.1) & (yy['ratio'] > cutoff1)
            & (yy['ratio.ll']>cutoff2) & (yy['ratio.h']>cutoff2) & (yy['ratio.v']>cutoff2)]
    return u
    
def SnapFISH25_step6(step5_result,out_dir):
    y0 = step5_result
    #read.table('101522_loop_candidates.txt', head=T)
    chr_name=sorted(list(set(y0['chr1'])))

    final = pd.DataFrame()
    for CHR in chr_name:
        y = y0[ y0['chr1']==CHR]
        y=y.reset_index(drop=True)
        GAP = 25e3
    # step 1: find neighborhood
        for i in range(0,len(y)):
            z = y[ (abs(y['x1'] - y['x1'][i])<= GAP) & (abs(y['y1'] - y['y1'][i]) <= GAP)]
            y.loc[i,'CountNei'] = len(z)

    # step 2: split singletons and peak clusters with >=3 bin pairs
        v0 = y[ y['CountNei'] == 1]
        v = y[ y['CountNei'] >= 2] # peak cluster: sharp peak + broad peak

    # step 3: keep singletons
        out = pd.DataFrame()

        if len(v0)>=1: # add singletons

            v0['label'] = np.nan
            v0['ClusterSize'] = 1
            v0['summit'] = np.nan
            v0['NegLog10FDR'] = -np.log10(v0['FDR'])
            out = pd.concat([out.reset_index(drop=True), v0])



    # step 4: for clusters, find summit
        if len(v)>0:
      
       # for peak cluster, assign label 1, 2, ..., N
            v['label'] = np.array(range(1,len(v)+1))

       # for all bin pairs within the neighborhood
       # assign the same cluster label, using the minimal label
            for i in range(0,len(v)):
          
                w = v[ (abs(v['x1'] - v['x1'][i])<=GAP) & (abs(v['y1'] - v['y1'][i])<=GAP)]
                w_min = min(w['label'])
                w_label = sorted(list(set(w['label'])))
                
                for j in range(1,len(w_label)):
              
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


    ming = final[ (final['summit']== 1) & (final['ClusterSize']>=2)]
    jie  = final[ final['ClusterSize'] == 1]



    out = pd.concat([ming,jie])
    out=out.sort_values('chr1')
    
    newfilename1 = os.path.join(out_dir, "Final_loop_summits.txt")
    out.to_csv(newfilename1, header=True, index=None, sep='\t')
    
def run_step123(chrom,dir_in_r1,dir_in_r2,ann_dir,out_dir):
    SnapFISH25_step1(dir_in_r1,dir_in_r2,chrom,out_dir)
    x,y = SnapFISH25_step2(chrom, out_dir)
    SnapFISH25_step3(chrom, ann_dir, x, y, out_dir)



def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--indirSF', action = 'store', required = True, \
                        help = 'SnapFISH directory')
    parser.add_argument('--rep1', action = 'store', required = True, \
                        help = 'input rep 1 data directory for 25Kb')
    parser.add_argument('--rep2', action = 'store', required = True, \
                        help = 'input rep 2 data directory for 25Kb')
    parser.add_argument('-o', '--outdir', action = 'store', \
                        required = True, help = 'output directory')
    parser.add_argument('-a', '--anndir', action = 'store', \
                        required = True, help = 'file list directory')
    parser.add_argument('-d', '--dataname', action = 'store', help = 'data name', \
                        required = True)
    parser.add_argument('-n', '--num_cpus', type = int , action = 'store', required = True, \
                        help = 'number of CPUS used for parallel computing')

    return parser
    

if __name__ == "__main__":
    main()

