import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from src.SingleChrProcessor import SingleChrProcessor

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

        assert self.BINSIZE == 5e3 or self.BINSIZE == 25e3, "resolution not supported"

        if not os.path.exists(path):
            os.mkdir(path)


    def run_SnapFISH(self):
        self.SnapFISH_step1_2()
        self.SnapFISH_step3()
        self.SnapFISH_step4()

    
    def SnapFISH_step1(self):
        # separate by chromosome number and calculate pairwise 
        # distances for each chromosome separately
        for chr_num, sub_data in self.data.groupby(self.chr_col):
            sub_ann = self.ann[self.ann[self.chr_col]==chr_num]
            scp = SingleChrProcessor(
                data=sub_data.reset_index(drop=True), 
                ann=sub_ann.reset_index(drop=True), 
                suf=str(chr_num), 
                sf=self
            )
            scp.SnapFISH_step1()

    
    def SnapFISH_step1_2(self):
        # group the 3D coordinates and annotation file by chromosome
        # calculate 3D distances and conduct tests for each chromosome 
        # separately then return the results
        tests_ls = []
        for chr_num, sub_data in self.data.groupby(self.chr_col):
            sub_ann = self.ann[self.ann[self.chr_col]==chr_num]
            scp = SingleChrProcessor(
                data=sub_data.reset_index(drop=True), 
                ann=sub_ann.reset_index(drop=True), 
                suf=str(chr_num), 
                sf=self
            )
            test_df = scp.single_chr_wrapper()
            test_df.insert(0, "chr", chr_num)
            tests_ls.append(test_df.reset_index(drop=True))
        self.out_all_tests = pd.concat(tests_ls)

        if self.write:
            # round columns of step 2 output and store as a single file
            out_test = self.out_all_tests.copy()
            round_cols = ["out.mean", "contact.freq", "case", "ctrl", "t.stat", "wilcox.stat"]
            out_test[round_cols] = out_test[round_cols].round(4)
            test_path = os.path.join(self.path, f"output_test_{self.suf}.txt")
            out_test.to_csv(test_path, index=False, sep="\t")


    def SnapFISH_step3(self, test_method = "t", fdr_cut = 0.1):
        BOUND_LO, BOUND_UP = 1e5, 1e6
        bin_dist = np.abs(self.out_all_tests["end2"] - self.out_all_tests["end1"])
        left_df = self.out_all_tests[(bin_dist>=BOUND_LO)&(bin_dist<=BOUND_UP)]
        # require FDR < 0.1 and negetive test statistic
        if test_method == "t":
            candidates = left_df[(left_df["t.fdr"] < fdr_cut)&(left_df["t.stat"] < 0)]
        elif test_method == "w":
            candidates = left_df[(left_df["t.fdr"] < fdr_cut)&(left_df["t.stat"] < 0)]
        candidates = candidates.reset_index(drop = True)

        # generate unique ID for each bin
        id1 = list(zip(candidates["chr"], candidates["bin1"]))
        id2 = list(zip(candidates["chr"], candidates["bin2"]))
        id_ann = list(zip(self.ann["chr"], self.ann[self.pos_col]))
        bin_locs = self.ann.copy()
        bin_locs.index = id_ann

        # add chromosome and 1D genomic location information
        kept_cols = ["chr", "start", "end"]
        bin1_df = bin_locs.loc[id1, kept_cols]
        bin1_df = bin1_df.rename({
            "chr":"chr1", "start":"x1", "end":"x2"
        }, axis=1)
        bin2_df = bin_locs.loc[id2, kept_cols]
        bin2_df = bin2_df.rename({
            "chr":"chr2", "start":"y1", "end":"y2"
        }, axis=1)
        self.out_candidate = pd.concat([
            bin1_df.reset_index(drop = True),
            bin2_df.reset_index(drop = True),
            candidates.drop("chr", axis = 1)
        ], axis = 1)

        if self.write:
            # round floats to 4 digits
            out_candidate = self.out_candidate.copy()
            round_cols = ["out.mean", "contact.freq", "case", "ctrl", "t.stat", "wilcox.stat"]
            out_candidate[round_cols] = out_candidate[round_cols].round(4)
            candidate_path = os.path.join(self.path, f"output_loop_candidate_{self.suf}.txt")
            out_candidate.to_csv(candidate_path, index=False, sep="\t")

        if self.save_pic:
            self.step3_4_heatmap("candidate")


    def spread_label(self, row, df):
        # 20Kb gap for  5Kb bin resolution
        # 50Kb gap for 25Kb bin resolution
        GAP = 20e3 if self.BINSIZE < 20e3 else 50e3
        neighbors = df[(abs(df['x1']-row['x1'])<=GAP) & (abs(df['y1']-row['y1'])<=GAP)]
        nan_row_idx = neighbors[pd.isna(neighbors["label"])].index
        # recursively call it self until no un-labled neighbors
        if len(nan_row_idx) != 0:
            df.loc[nan_row_idx, "label"] = row["label"]
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
            self.step3_4_heatmap("summit")


    def step3_4_heatmap(self, plot_ty):
        # plot the loop candidates/summit for each chromosome
        if plot_ty == "candidate":
            data = self.out_candidate
        elif plot_ty == "summit":
            data = self.out_summit

        for chr, sub_summit in data.groupby("chr1"):
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
            if plot_ty == "candidate":
                out_path = os.path.join(self.path, f"heatmap_loop_candidates_{chr_suf}.pdf")
            elif plot_ty == "summit":
                out_path = os.path.join(self.path, f"heatmap_loop_summits_{chr_suf}.pdf")
            plt.savefig(out_path)
            plt.clf() # clear previous plot
