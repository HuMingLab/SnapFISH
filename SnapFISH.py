import argparse
import sys
import os
from itertools import combinations
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.stats.multitest as multi
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

class SnapFISH:
    def __init__(self, coor_path, ann_path, path, suf, save_pic):
        self.path, self.suf = path, suf
        self.id_col, self.pos_col = "cell_id", "pos"
        self.data = pd.read_csv(coor_path, sep="\t")
        self.ann = pd.read_csv(ann_path, sep="\t")
        self.BINSIZE = self.ann['end'][0] - self.ann['start'][0]
        self.save_pic = save_pic

    def run_SnapFISH(self, paired=False):
        self.SnapFISH_step1()
        self.SnapFISH_step2(paired)
        self.SnapFISH_step3()
        self.SnapFISH_step4()

    def SnapFISH_step1(self):
        CELL_IDS = np.unique(self.data[self.id_col])
        NUM_POS = self.ann.shape[0]

        cell_ids = np.repeat(CELL_IDS, NUM_POS)
        pos_ids = np.tile(self.ann[self.pos_col], len(CELL_IDS))
        missed_rows = set(zip(cell_ids, pos_ids)).difference(
            set(zip(self.data[self.id_col], self.data[self.pos_col])))
        if len(missed_rows) != 0:
            missed_rows = np.array(list(missed_rows))
            missed_df = pd.DataFrame(missed_rows, columns=[self.id_col, self.pos_col])
            self.data = pd.concat([self.data, missed_df])
        self.data.sort_values([self.id_col, self.pos_col], inplace=True)
        self.data = self.data.reset_index(drop=True)

        def apply_dist(x): # calculate the distance matrix for each cell
            diff = np.array([*combinations(x[["x", "y", "z"]].values, 2)])
            return np.sqrt(np.sum(np.square(diff[:,0,:] - diff[:,1,:]), axis=1))
        dists = self.data.groupby(self.id_col).apply(apply_dist)
        self.final = pd.DataFrame(np.array(dists.values.tolist()).T, columns=dists.index)
        bins = combinations(self.ann[self.pos_col], 2)
        self.bin_cols = pd.DataFrame(list(bins), columns=["bin1", "bin2"])

        self.out_3D_dist = pd.concat([self.bin_cols, self.final], axis=1)
        out_path = os.path.join(self.path, f"output_3D_dist_{self.suf}.txt")
        self.out_3D_dist.round(4).to_csv(out_path, sep="\t", index=False)

        mean_dist = self.final.mean(axis=1, skipna=True).to_frame("out.mean")
        out_3D_dist_mean = pd.concat([self.bin_cols, mean_dist], axis=1)
        mean_path = os.path.join(self.path, f"output_3D_dist_{self.suf}_avg.txt")
        out_3D_dist_mean.round(4).to_csv(mean_path, sep="\t", index=False)
        
        if (self.save_pic == 0):
          NUM = self.ann.shape[0] # number of loci
          a = out_3D_dist_mean
          x=np.zeros((NUM,NUM))
          a[['bin1','bin2']] -= 1
      
          for i in range(0,a.shape[0]):
              x[ int(a['bin1'][i]), int(a['bin2'][i]) ] = a['out.mean'][i]
              x[ int(a['bin2'][i]), int(a['bin1'][i]) ] = a['out.mean'][i]
      
          np.fill_diagonal(x, np.nan)
      
          plt.imshow(x, cmap='RdBu')
          plt.colorbar()
          newfilename1 = os.path.join(out_dir, "heatmap_av_Euc_dist_%s.pdf" % (data_name))
          plt.savefig(newfilename1)
          plt.clf()

    def SnapFISH_step2(self, paired=False):
        CUT_UP = 50e3
        CUT_LO = 20e3 if self.BINSIZE == 5e3 else 25e3
        L, U = CUT_UP/self.BINSIZE, CUT_LO/self.BINSIZE

        t_test_ls = []
        for i in range(self.out_3D_dist.shape[0]):
            v = self.final.iloc[i,:].dropna()
            b1, b2 = self.bin_cols['bin1'][i], self.bin_cols['bin2'][i]
            bin1 = abs(self.bin_cols['bin1'] - b1)
            bin2 = abs(self.bin_cols['bin2'] - b2)
            filter_x = (bin1 <= L)&(bin2 <= L)&~((bin1 <= U)&(bin2 <= U))
            x_bins = self.bin_cols[filter_x]
            x_vals = self.final[filter_x]
            
            # local background
            xll = x_vals[(x_bins['bin1'] > b1)&(x_bins['bin2'] < b2)].mean().dropna() # lower left
            xh  = x_vals[x_bins['bin1'] == b1].mean().dropna() # horizontal
            xv  = x_vals[x_bins['bin2'] == b2].mean().dropna() # vertical

            x_mean = x_vals.mean().dropna()
            if not paired:
                t_stat, p_val = stats.ttest_ind(v, x_mean, equal_var=False)
            else:
                t_stat, p_val = stats.ttest_rel(v, x_mean)
            t_test_ls.append(list(map(np.nanmean, [v, x_mean, xll, xh, xv])) + [t_stat, p_val])

        t_cols = ['Case', 'Ctrl', 'Ctrl.ll', 'Ctrl.h', 'Ctrl.v', 'Tstat', 'Pvalue']
        t_test_df = pd.DataFrame(t_test_ls, columns=t_cols)
        pval_corr = multi.multipletests(t_test_df['Pvalue'], method='fdr_bh')[1]
        t_test_df["fdr"] = pval_corr
        t_cols = ['Case', 'Ctrl', 'Tstat', 'Pvalue', 'Ctrl.ll', 'Ctrl.h', 'Ctrl.v','fdr']

        self.out_Ttest = pd.concat([self.bin_cols, t_test_df[t_cols]], axis=1)
        Ttest_path = os.path.join(self.path, f"output_Ttest_{self.suf}.txt")
        self.out_Ttest.to_csv(Ttest_path, index=False, sep="\t")

    def SnapFISH_step3(self):
        a = self.out_Ttest
        # only consider bin pairs 100Kb ~ 1MB
        BOUND_LO, BOUND_UP = 1e5, 1e6
        # user defined cutoff
        cutoff1, cutoff2 = 1.1, 1.05

        a = a[ (a['bin2']-a['bin1']>=BOUND_LO/self.BINSIZE) & (a['bin2']-a['bin1']<=BOUND_UP/self.BINSIZE)
            & (~np.isnan(a['Ctrl'])) & (~np.isnan(a['Ctrl.ll'])) & (~np.isnan(a['Ctrl.h'])) & (~np.isnan(a['Ctrl.v']))]

        x = a[ (a['Tstat'] < -4) & (a['fdr'] < 0.1) & ((a['Ctrl']/a['Case']) > cutoff1)
            & ((a['Ctrl.ll']/a['Case']) > cutoff2) & ((a['Ctrl.h']/a['Case']) > cutoff2) & ((a['Ctrl.v']/a['Case']) > cutoff2)]

        x = x.iloc[np.argsort(x['fdr']),].reset_index(drop=True)

        bins = x[["bin1", "bin2"]].values.flatten()
        bin_locs = self.ann.set_index(self.pos_col).loc[bins, ["chr", "start", "end"]].values
        rec = pd.DataFrame(bin_locs.reshape((-1, 6)), 
                           columns = ['chr1','x1','x2','chr2','y1','y2'])

        self.out_candidate = pd.concat([rec, x], axis=1)
        candidate_path = os.path.join(self.path, f"output_loop_candidate_{self.suf}.txt")
        self.out_candidate.to_csv(candidate_path, index=False, sep="\t")

    def SnapFISH_step4(self):
        GAP = 10e3 if self.BINSIZE == 5e3 else 50e3
        def spread_label(row, df):
            neighbors = df[(abs(df['x1']-row['x1'])<=GAP) & (abs(df['y1']-row['y1'])<=GAP)]
            nan_row_idx = neighbors[pd.isna(neighbors["label"])].index
            if len(nan_row_idx) != 0:
                df.loc[nan_row_idx, "label"] = row["label"]
                for _, row in df.loc[nan_row_idx].iterrows():
                    spread_label(row, df)
        def apply_cluster(x):
            x[["NegLog10FDR","ClusterSize"]] = sum(-np.log10(x["fdr"])),len(x)
            return x.iloc[[np.argmin(x["Tstat"])]]

        result = []
        for _, df in self.out_candidate.groupby("chr1"):
            df["label"] = np.nan
            while df["label"].isna().any():
                na_rows = df.index[pd.isna(df["label"])]
                if len(na_rows) == len(df):
                    df.loc[na_rows[0], "label"] = 1
                else:
                    df.loc[na_rows[0], "label"] = np.max(df["label"]) + 1
                spread_label(df.loc[na_rows[0]], df)
            s_df = df.groupby("label", group_keys=False).apply(apply_cluster)
            result.append(s_df[s_df["ClusterSize"]>1])
        summit_cols = self.out_candidate.columns.tolist()+["NegLog10FDR", "ClusterSize"]
        self.out_summit = pd.concat(result) if len(result) != 0 \
            else pd.DataFrame(columns=summit_cols)

        summit_path = os.path.join(self.path, f"output_loop_summit_{self.suf}.txt")
        self.out_summit.to_csv(summit_path, index=False, sep="\t")
        
        if (self.save_pic == 0):
          NUM = self.ann.shape[0]
          x = self.out_summit
          u = np.zeros((NUM,NUM))
          for i in range(0,x.shape[0]):
              u[ int(x['bin1'][i]), int(x['bin2'][i]) ] = 1
              u[ int(x['bin2'][i]), int(x['bin1'][i]) ] = 1
              
      
          plt.imshow(u, cmap='Reds')
          plt.colorbar()
          newfilename1 = os.path.join(out_dir, "heatmap_loop_summits_%s.pdf" % (data_name))
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
    parser = create_parser()
    args = parser.parse_args()
    
    SnapFISH_dir = args.indirSF
    dir_in = args.indir
    ann_dir = args.anndir
    data_name = args.dataname
    
    out_dir = args.outdir
    
    save_pic = args.pic

    sf = SnapFISH(dir_in, ann_dir, out_dir, data_name, save_pic)
    sf.run_SnapFISH()
