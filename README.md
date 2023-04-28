# SnapFISH: a computational pipeline to identify chromatin loops from multiplexed DNA FISH data

Find the preprint of SnapFISH [here](https://www.biorxiv.org/content/10.1101/2022.12.16.520793v1).

## Installation

SnapFISH is a published package on PyPI and can be downloaded by `pip install SnapFISH`. It is recommended to create a new conda environment before installing SnapFISH. This can be done by 

```bash
conda create --name SnapFISH_env python==3.6.8
conda activate SnapFISH_env
pip install SnapFISH
```

Alternatively, you can download the source code from this github site and run it using the shell script:

```bash
git clone https://github.com/lindsayhrlee/SnapFISH SnapFISH && cd SnapFISH
pip install -r requirements.txt
```

## Program Descriptions

The SnapFISH algorithm consists of four steps:

* Step 1a: take 3D coordinates of each locus in each cell as input, calculate Eucildean disance between any two loci pairs, and also calculate average Euclidean disance for each loci pair across all available cells to create distance matrix.          
* Step 1b: make heatmap of average Euclidean distance and contact frequencies                    
* Step 2: perform two sample T-test or Wilcoxon rank test for each loci pair based on local background model               
* Step 3a: identify loop candidates        
* Step 3b: make heatmap of loop candidates               
* Step 4a: group nearby loop candidates into clusters. remove singletons, only keep summits              
* Step 4b: make heatmap of loop summits

## Use SnapFISH As a Command Line Tool

You can use SnapFISH as a command line tool after activating the conda environment and installing SnapFISH by `pip install SnapFISH`. 

### Identification of Enhancer-Promoter Loops

To call enhancer-promoter loops, type the following in your terminal:

```bash
SnapFISH call-loops -i PATH/TO/3DCOOR.txt -o OUT/DIR -a PATH/TO/ANN.txt -p 1 -d SUFFIX
```

where

1. -i PATH/TO/3DCOOR.txt : the txt file containing the following columns (separated by tab)
  * Col #1 `chr`: chromosome number
  * Col #2 `cell_id`: unique ID for each cell
  * Col #3 `pos`: locus ID
  * Col #4~#6 `x`, `y`, `z`: 3D coordinates (X,Y,Z) for each bin, unit: nm.

2. -o OUT/DIR : the directory for output

3. -a PATH/TO/ANN.txt : the 1D genomic annotation file containing the genomic location for each bin
  * Col #1 `chr`: chromosome number
  * Col #2 `start`: the starting 1D genomic location for the corresponding bin specified by Col #4
  * Col #3 `end`: the ending 1D genomic location for the corresponding bin specified by Col #4
  * Col #4 `pos`: the unique ID for each bin

4. -p 1 : whether to save the heatmaps produced in step 1b, 3b, and 4b. 0 or 1, 1 to save.

5. -d SUFFIX : name of data for suffix of the output files

For windows, the file path needs to be separated by `\\` instead of `/`.

For example, running 

```bash
SnapFISH call-loops -i ext/129_3d_coor.txt -o tmp -a ext/input_ann.txt -p 1 -d 129
```

will print 

```
129chr3 cutoff: 349.4659
Pearson's r (average.dist and contact.freq): -0.9845
```

where the cutoff is the average Euclidean distance in nm between two bin pairs 25kb apart. The Pearson's r value is the correlation between the average distance matrix and the contact frequency matrix. The following files will also be produced

1. tmp/output_3D_dist_129chr3.txt: the pairwise distance between each bin in each cell. Missing values are skipped.

2. tmp/output_3D_dist_129chr3_avg.txt: the average pairwise distance between each bin across all cells. 

3. tmp/heatmap_av_Euc_dist_129chr3.pdf, tmp/heatmap_Contact_freq_129chr3.pdf: heatmaps of average pairwise distance and contact frequency.

5. tmp/output_test_129.txt: T-test and Wilcoxon rank test results from the local background model.

6. tmp/output_loop_candidate_129.txt, tmp/heatmap_loop_candidates_129chr3.pdf: loop candidates selected and the heatmap showing their positions.

7. tmp/output_loop_summit_129.txt, tmp/heatmap_loop_summits_129chr3.pdf: loop summits identified and the heatmap showing their positions.

### Pairwise Distances

In some cases, users might only be interested in getting the pairwise distances between each bin without actually calling loops. This can be done by

```bash
SnapFISH step1 -i PATH/TO/3DCOOR.txt -o OUT/DIR -a PATH/TO/ANN.txt -p 1 -d SUFFIX
```

The only outputs will be pairwise distances and the two associated heatmaps.

### Wilcoxon Rank Test

By default, SnapFISH identifies loop candidates by the T-test. However, if the underlying distribution deviates from the normal distribution significantly, users can also filter loops by the Wilcoxon rank test. After running the `call-loops` command, the user can adjust the filtering criteria by

```bash
SnapFISH test -i OUT/DIR/output_test_SUF.txt -o OUT/DIR -a PATH/TO/ANN.txt -p 1 -d SUFFIX -m w -h 0.1
```

where 

1. -i OUT/DIR/output_test_SUF.txt: the output from step2.

2. -m w: if `w`, will use the FDR from Wilcoxon rank tests; if `t`, will use the FDR from T-tests. Default is `t`.

3. -h 0.1: the filtering threshold of FDR. Only bin pairs with FDR less than the threshold will be kept. Default is 0.1.

For example,

```bash
SnapFISH -i tmp/output_test_129.txt -o tmp -a ext/input_ann.txt -p 1 -d 129 -m w -h 0.05
```

will only keep the bin pairs with FDR less than 0.05 from the Wilcoxon rank test.

## Run SnapFISH with Shell

SnapFISH can also be run by downloading the files in this github page and run the shellscript. The **required** input variables of the shell script are:

1. SnapFISH_dir : The directory of SnapFISH
2. input_path : The path of the file containing 3D coordinates.
3. output_dir : The directory for output
4. ann_file : The path of the file containing genomic location for each bin.
5. save_pic : 1 for outputting heatmaps, and 0 for not
6. data_name : Name of data for prefix

To run SnapFISH,
```
./run_SnapFISH.sh
```

If you are using Windows, please change `/` to `\\`.

## Requirements
SnapFISH was built using following Python packages.

1. Python 3.6.8
2. numpy 1.19.5
3. pandas 1.1.5
5. scipy 1.5.4
6. statsmodels 0.12.2
7. matplotlib 3.3.4

## Details

- Chromatin tracing data from Huang et al, NG, 2021 (PMID: 34002095)
- 5Kb bin resolution, 41 loci cover a 205Kb region in mESC *Sox2* locus. chr3:34,601,078–34,806,078 (ref: mm10).
- Cast allele contains the 7.5Kb 4CBS insertion.
- 129 allele does not contain the 7.5Kb 4CBS insertion. 
- After filtering, we keep 649 cells (649 129 alleles and 649 CAST alleles). On average, each allele contains 28 5Kb bins.
- SnapFISH outputs the list of summit. 

## Contact Us
For any questions regarding this software, contact Ming Hu (hum@ccf.org) or Lindsay Lee (leeh7@ccf.org).
