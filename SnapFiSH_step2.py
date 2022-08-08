import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.stats.multitest as multi


a = pd.read_csv('063022_output_3D_dist_CAST.txt',sep='\t')
ann = pd.read_csv('063022_input_ann.txt',sep='\t')
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
output.to_csv('output_Ttest_CAST.txt', header=True, index=None, sep='\t')
