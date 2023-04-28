import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import combinations
from scipy import stats
from scipy.spatial.distance import pdist
import statsmodels.stats.multitest as multi

class SingleChrProcessor:
    def __init__(self, data, ann, suf, sf):
        # add suffix so that files won't overwrite each other
        self.suf = sf.suf + "chr" + suf
        self.sf = sf
        # data and annotation file of a single chromosome
        self.data = data
        self.ann = ann
        self.BINSIZE = ann["end"][0] - ann["start"][0]


    def single_chr_wrapper(self):
        self.SnapFISH_step1()
        self.SnapFISH_step2(self.sf.paired)
        return self.out_test


    def insert_missing_rows(self):
        CELL_IDS = np.unique(self.data[self.sf.id_col])
        NUM_POS = self.ann.shape[0]

        cell_ids = np.repeat(CELL_IDS, NUM_POS)
        pos_ids = np.tile(self.ann[self.sf.pos_col], len(CELL_IDS))
        full_rids = set(zip(cell_ids, pos_ids))
        raw_rids = set(zip(self.data[self.sf.id_col], self.data[self.sf.pos_col]))
        # expected (cell id, locus id) - observed (cell id, locus id)
        missed_rows = full_rids - raw_rids
        
        if len(missed_rows) != 0:
            missed_rows = np.array(list(missed_rows))
            missed_df = pd.DataFrame(
                missed_rows, 
                columns=[self.sf.id_col, self.sf.pos_col]
            )
            missed_df[self.sf.pos_col] = missed_df[self.sf.pos_col].astype("int")
            self.data = pd.concat([self.data, missed_df])
        
        # insert all missing cell id and locus id
        self.data = self.data.sort_values([self.sf.id_col, self.sf.pos_col])
        self.data = self.data.reset_index(drop=True)

    
    def step1_heatmap(self, plot_ty):
        NUM = self.ann.shape[0]
        heat_ut = self.out_3D_dist_mean.pivot_table(
            plot_ty, "bin1", "bin2", dropna=False
        ).values
        x = np.diag(NUM*[0.0])
        x[np.triu_indices(NUM, 1)] = heat_ut[np.triu_indices_from(heat_ut)]
        x[np.tril_indices(NUM, -1)] = heat_ut.T[np.tril_indices_from(heat_ut)]
        np.fill_diagonal(x, np.nan)

        plt.imshow(x, cmap='RdBu')
        plt.colorbar()
        if plot_ty == "out.mean":
            heatmap_name = os.path.join(self.sf.path, f"heatmap_av_Euc_dist_{self.suf}.pdf")
        elif plot_ty == "contact.freq":
            heatmap_name = os.path.join(self.sf.path, f"heatmap_Contact_freq_{self.suf}.pdf")
        plt.savefig(heatmap_name)
        plt.clf()


    def calculate_contact_freqs(self):
        # use average spatial distance between two loci with 
        # 1D distance 25Kb as the cutoff for contact
        contact_dist = 25e3

        bin_dists = self.out_3D_dist["end2"] - self.out_3D_dist["end1"]
        filtered_bins = self.final[bin_dists.abs() == contact_dist]
        cutoff = np.nanmean(filtered_bins)
        # print the calculated cutoff
        print(f"{self.suf} cutoff: {round(cutoff, 4)}")

        contact_count = ((~self.final.isna())&(self.final < cutoff)).sum(axis = 1)
        notna_count = (~self.final.isna()).sum(axis = 1)
        return contact_count/notna_count


    def SnapFISH_step1(self):
        self.insert_missing_rows()
        # calculate the pairwise distances within each cell
        dists = self.data.groupby(self.sf.id_col).apply(
            lambda x: pdist(x[["x", "y", "z"]].values)
        )
        dist_vals = np.array(dists.values.tolist())
        self.final = pd.DataFrame(
            dist_vals.T, 
            columns=dists.index
        )

        # locus ID column, also 1D genomic locations since there
        # might be gaps between bins
        bins = np.array([*combinations(self.ann[["pos", "start"]].values, 2)])
        self.bin_cols = pd.DataFrame(
            bins.transpose(0, 2, 1).reshape((-1, 4)), 
            columns=["bin1", "bin2", "end1", "end2"]
        )
        self.out_3D_dist = pd.concat([self.bin_cols, self.final], axis=1)

        # calculate the mean distance and contact frequency for each bin pair
        mean_dist = self.final.mean(axis=1, skipna=True).to_frame("out.mean")
        self.out_3D_dist_mean = pd.concat([self.bin_cols, mean_dist], axis=1)
        self.out_3D_dist_mean["contact.freq"] = self.calculate_contact_freqs()

        mean_wo_na = self.out_3D_dist_mean.dropna()
        corr = stats.pearsonr(mean_wo_na["out.mean"], mean_wo_na["contact.freq"])
        print("Pearson's r (average.dist and contact.freq):", round(corr[0], 4))

        if self.sf.save_pic:
            self.step1_heatmap("out.mean") # pairwise distances
            self.step1_heatmap("contact.freq") # contact frequencies

        # return the 3D distances and mean distances if write is False
        if not self.sf.write:
            return self.out_3D_dist, self.out_3D_dist_mean
        out_path = os.path.join(self.sf.path, f"output_3D_dist_{self.suf}.txt")
        self.out_3D_dist.round(4).to_csv(out_path, sep="\t", index=False)
        mean_path = os.path.join(self.sf.path, f"output_3D_dist_{self.suf}_avg.txt")
        self.out_3D_dist_mean.round(4).to_csv(mean_path, sep="\t", index=False)

    
    def test_iterrow(self, CUT_UP, CUT_LO, NUM_NBRS, paired, i):
        row_i = self.final.iloc[i,:].dropna()

        # calculate the 1D distances between the ith bin pair and all other bins
        end1_i = self.bin_cols["end1"][i]
        bin1_dist = abs(self.bin_cols["end1"] - end1_i)
        end2_i = self.bin_cols["end2"][i]
        bin2_dist = abs(self.bin_cols["end2"] - end2_i)

        # filter by 1D distances
        filter_up = (bin1_dist <= CUT_UP)&(bin2_dist <= CUT_UP)
        filter_lo = (bin1_dist <= CUT_LO)&(bin2_dist <= CUT_LO)
        filter_x = filter_up&(~filter_lo)
        x_bins = self.bin_cols[filter_x]

        # local background
        x_mean = self.final[filter_x].mean().dropna()

        # skip if no enough background
        low_resol = self.BINSIZE >= 20e3 and len(x_bins) != NUM_NBRS
        high_resol = self.BINSIZE < 20e3 and len(x_bins) <= 2
        not_enough_nbr = low_resol or high_resol
        if len(x_mean) == 0 or len(row_i) == 0 or not_enough_nbr:
            return []

        if not paired:
            t_stat, p_val_t = stats.ttest_ind(row_i, x_mean, equal_var=False)
            w_stat, p_val_w = stats.ranksums(row_i, x_mean)
        else:
            t_stat, p_val_t = stats.ttest_rel(row_i, x_mean)
            w_stat, p_val_w = stats.wilcoxon(row_i, x_mean)

        return [
            *self.out_3D_dist_mean.iloc[i].tolist(),
            len(x_bins), np.nanmean(row_i), np.mean(x_mean), 
            t_stat, p_val_t, w_stat, p_val_w
        ]


    def SnapFISH_step2(self, paired=False):
        CUT_UP, CUT_LO = 50e3, 25e3
        # for 25Kb resolution, NUM_NBRS = 5^2-3^2=16
        # for  5Kb resolution, NUM_NBRS = 21^2-11^2 = 320
        NUM_NBRS = (CUT_UP/self.BINSIZE*2 + 1)**2 - (CUT_LO/self.BINSIZE*2 + 1)**2

        test_ls = [] # iterate through all bin pairs
        for i in range(self.out_3D_dist.shape[0]):
            test_row = self.test_iterrow(CUT_UP, CUT_LO, NUM_NBRS, paired, i)
            if len(test_row) != 0:
                test_ls.append(test_row)

        self.out_test = pd.DataFrame(test_ls, columns=[
            *self.out_3D_dist_mean.columns.tolist(),
            "num.circle", "case", "ctrl", "t.stat",
            "t.pval", "wilcox.stat", "wilcox.pval",
        ])
        # adjust for multiple comparisons
        self.out_test["t.fdr"] = multi.multipletests(
            self.out_test["t.pval"], method="fdr_bh"
        )[1]
        self.out_test["wilcox.fdr"] = multi.multipletests(
            self.out_test["wilcox.pval"], method="fdr_bh"
        )[1]
        self.out_test = self.out_test.astype(self.out_3D_dist_mean.dtypes)