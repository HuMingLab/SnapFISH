import os, argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import combinations
from scipy import stats
from scipy.spatial.distance import pdist
import statsmodels.stats.multitest as multi

class SnapFISH:
    def __init__(self, coor_path, ann_path, path, suf, 
                 w=True, save_pic=False, paired=False):
        self.path = path
        self.suf = suf
        self.write = w
        self.save_pic = save_pic
        self.id_col = "cell_id"
        self.pos_col = "pos"
        self.chr_col = "chr"
        self.paired = paired

        self.data = pd.read_csv(coor_path, sep="\t")
        self.ann = pd.read_csv(ann_path, sep="\t")
        self.BINSIZE = self.ann["end"][0] - self.ann["start"][0]

        if not os.path.exists(path):
            os.mkdir(path)


    def run_SnapFISH(self):
        self.separate_by_chrs()
        self.SnapFISH_step4()
        

    def separate_by_chrs(self):
        # group the 3D coordinates and annotation file by chromosome
        # calculate 3D distances, conduct tests, and identify loop
        # candidates for each chromosome separately then return the
        # results
        candidate_ls = []
        for chr_num, sub_data in self.data.groupby(self.chr_col):
            sub_ann = self.ann[self.ann[self.chr_col]==chr_num]
            scp = SingleChrProcessor(
                data=sub_data.reset_index(drop=True), 
                ann=sub_ann.reset_index(drop=True), 
                suf=str(chr_num), 
                sf=self
            )
            candidate_ls.append(scp.single_chr_wrapper())
        self.out_candidate = pd.concat(candidate_ls)

        if not self.write:
            return self.out_candidate
        # round floats to 4 digits
        out_candidate = self.out_candidate.copy()
        round_cols = ["out.mean", "contact.freq", "case", "ctrl", "t.stat", "wilcox.stat"]
        out_candidate[round_cols] = out_candidate[round_cols].round(4)
        candidate_path = os.path.join(self.path, f"output_loop_candidate_{self.suf}.txt")
        out_candidate.to_csv(candidate_path, index=False, sep="\t")


    def spread_label(self, row, df):
        # 20Kb gap for  5Kb bin resolution
        # 50Kb gap for 25Kb bin resolution
        GAP = 20e3 if self.BINSIZE == 5e3 else 50e3
        neighbors = df[(abs(df['x1']-row['x1'])<=GAP) & (abs(df['y1']-row['y1'])<=GAP)]
        nan_row_idx = neighbors[pd.isna(neighbors["label"])].index
        if len(nan_row_idx) != 0:
            df.loc[nan_row_idx, "label"] = row["label"]
            # recursively call it self until no un-labled neighbors
            for _, row in df.loc[nan_row_idx].iterrows():
                self.spread_label(row, df)


    def apply_cluster(self, x):
        x["ClusterSize"] = len(x)
        x["NegLog10FDR"] = sum(-np.log10(x["t.fdr"]))
        # keep only the one with the most significant test statistic 
        # in each cluster
        return x.iloc[[np.argmin(x["t.stat"])]]


    def SnapFISH_step4(self):
        result = []
        for _, df in self.out_candidate.groupby("chr1"):
            df["label"] = np.nan
            while df["label"].isna().any():
                na_rows = df.index[pd.isna(df["label"])]
                if len(na_rows) == len(df):
                    df.loc[na_rows[0], "label"] = 1 # first cluster
                else:
                    # new cluster, cluster id = max(all previous ID) + 1
                    df.loc[na_rows[0], "label"] = np.max(df["label"]) + 1
                # label all its neighbors
                self.spread_label(df.loc[na_rows[0]], df)
            # keep only the most significant bin pair
            s_df = df.groupby("label", group_keys=False).apply(self.apply_cluster)
            result.append(s_df)

        summit_cols = self.out_candidate.columns.tolist()+["ClusterSize", "NegLog10FDR"]
        if len(result) != 0:
            results = pd.concat(result)
            # filter by contact frequency
            # if singleton -> require contact frequency >= 3/6
            singletons = (results["contact.freq"] >= 3/6)&(results["ClusterSize"] == 1)
            # if cluster -> require contact frequency >= 2/6
            summits = (results["contact.freq"] >= 2/6)&(results["ClusterSize"] >= 2)
            self.out_summit = results[singletons|summits]
        else:
            self.out_summit = pd.DataFrame(columns=summit_cols)

        summit_path = os.path.join(self.path, f"output_loop_summit_{self.suf}.txt")
        if not self.write:
            return self.out_summit
        self.out_summit.to_csv(summit_path, index=False, sep="\t")

        if self.save_pic:
            self.step4_heatmap()


    def step4_heatmap(self):
        # plot the loop summit for each chromosome
        for chr, sub_summit in self.out_summit.groupby("chr1"):
            x = sub_summit[["bin1", "bin2"]]
            plt_pos = x.values - min(self.ann[self.pos_col])
            NUM = self.ann[self.ann["chr"]==chr].shape[0]
            u = np.diag(NUM*[0.0])
            for a, b in plt_pos:
                u[a, b] = 1
                u[b, a] = 1
            plt.imshow(u, cmap='Reds')
            plt.colorbar()
            chr_suf = self.suf + "chr" + str(chr)
            newfilename1 = os.path.join(self.path, "heatmap_loop_summits_%s.pdf" % (chr_suf))
            plt.savefig(newfilename1)
            plt.clf() # clear previous plot


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
        return self.SnapFISH_step3()


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

    
    def step1_heatmap(self):
        NUM = self.ann.shape[0]
        heat_ut = self.out_3D_dist_mean.pivot_table(
            "out.mean", "bin1", "bin2", dropna=False
        ).values
        x = np.diag(NUM*[0.0])
        x[np.triu_indices(NUM, 1)] = heat_ut[np.triu_indices_from(heat_ut)]
        x[np.tril_indices(NUM, -1)] = heat_ut.T[np.tril_indices_from(heat_ut)]
        np.fill_diagonal(x, np.nan)

        plt.imshow(x, cmap='RdBu')
        plt.colorbar()
        newfilename1 = os.path.join(self.sf.path, "heatmap_av_Euc_dist_%s.pdf" % (self.suf))
        plt.savefig(newfilename1)
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
            self.step1_heatmap()

        # return the 3D distances and mean distances if write is False
        if not self.sf.write:
            return self.out_3D_dist, self.out_3D_dist_mean
        out_path = os.path.join(self.sf.path, f"output_3D_dist_{self.suf}.txt")
        self.out_3D_dist.round(4).to_csv(out_path, sep="\t", index=False)
        mean_path = os.path.join(self.sf.path, f"output_3D_dist_{self.suf}_avg.txt")
        self.out_3D_dist_mean.round(4).to_csv(mean_path, sep="\t", index=False)

    
    def step2_iterrow(self, CUT_UP, CUT_LO, NUM_NBRS, paired, i):
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
            test_row = self.step2_iterrow(CUT_UP, CUT_LO, NUM_NBRS, paired, i)
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

        if not self.sf.write:
            return self.out_test
        out_test = self.out_test.copy()
        round_cols = ["out.mean", "contact.freq", "case", "ctrl", "t.stat", "wilcox.stat"]
        out_test[round_cols] = out_test[round_cols].round(4)
        Ttest_path = os.path.join(self.sf.path, f"output_test_{self.suf}.txt")
        out_test.to_csv(Ttest_path, index=False, sep="\t")


    def SnapFISH_step3(self):
        # only consider bin pairs 100Kb ~ 1MB
        BOUND_LO, BOUND_UP = 1e5, 1e6
        bin_dist = np.abs(self.out_test["end2"] - self.out_test["end1"])
        left_df = self.out_test[(bin_dist>=BOUND_LO)&(bin_dist<=BOUND_UP)]
        # require FDR < 0.1 and negetive test statistic
        candidates = left_df[(left_df["t.fdr"] < 0.1)&(left_df["t.stat"] < 0)]
        candidates = candidates.reset_index(drop = True)

        # add chromosome and 1D genomic location information
        bins = candidates[["bin1", "bin2"]].values.flatten()
        bin_locs = self.ann.set_index(self.sf.pos_col)
        bin_locs = bin_locs.loc[bins, ["chr", "start", "end"]]
        bin_locs = bin_locs.values.reshape((-1, 6))
        rec = pd.DataFrame(
            bin_locs, 
            columns = ["chr1", "x1", "x2", "chr2", "y1", "y2"]
        )
        return pd.concat([rec, candidates], axis=1)


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--indir', action = 'store', \
                        required = True, help = 'input directory')
    parser.add_argument('-o', '--outdir', action = 'store', \
                        required = True, help = 'output directory')
    parser.add_argument('-a', '--anndir', action = 'store', \
                        required = True, help = 'file list directory')
    parser.add_argument('-p', '--pic', action = 'store', type = int, \
                        required = True, help = '0 to save figure, and 1 to not')
    parser.add_argument('-d', '--dataname', action = 'store', \
                        required = True,  help = 'data name')
    return parser.parse_args()


if __name__ == "__main__":
    args = create_parser()
    sf = SnapFISH(
        coor_path = args.indir, 
        ann_path = args.anndir,
        path = args.outdir, 
        suf = args.dataname, 
        w = True,
        save_pic = bool(args.pic)
    )
    sf.run_SnapFISH()
