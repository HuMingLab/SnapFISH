# SnapFISH: single nucleus analysis pipeline for super-resolution imaging data

Find the preprint of SnapFISH [here](https://www.biorxiv.org/content/10.1101/2022.12.16.520793v1)

### Requirements
SnapFiSH was built using following Python packages.

1. Python 3.6.8
2. numpy 1.19.5
3. pandas 1.1.5
5. scipy 1.5.4
6. statsmodels 0.12.2
7. matplotlib.pyplot
8. math


### Details

- Chromatin tracing data from Huang et al, NG, 2021 (PMID: 34002095)
- 5Kb bin resolution, 41 loci cover a 205Kb region in mESC *Sox2* locus. chr3:34,601,078–34,806,078 (ref: mm10).
- Cast allele contains the 7.5Kb 4CBS insertion.
- 129 allele does not contain the 7.5Kb 4CBS insertion. 
- After filtering, we keep 649 cells (649 129 alleles and 649 CAST alleles). On average, each allele contains 28 5Kb bins.
- SnapFISH outputs the list of summit. 

### File Descriptions

* SnapFish.py       
  * Step 1a: take 3D coordinates of each locus in each cell as input, calculate Eucildean disance between any two loci pairs, and also calculate average Euclidean disance for each loci pair across all available cells to create distance matrix.          
  * Step 1b: make heatmap of average Euclidean distance                       
  * Step 2: perform two sample T-test for each loci pair based on local background model               
  * Step 3a: identify loop candidates        
  * Step 3b: make heatmap of loop candidates               
  * Step 4a: group nearby loop candidates into clusters. remove singletons, only keep summits              
  * Step 4b: make heatmap of loop summits
  
The input files are located in the ext folder.

* input_ann.txt 
  * Input data, genomic location for each bin. 
  * The columns consist of chr, start, end, and pos
  * Col #4 "pos" is the unique ID for each bin.

* 129_3d_coor.txt
  * Input data for the 129 allele. 
  * Col #1: chromosome number
  * Col #2: unique ID for each cell
  * Col #3: locus ID
  * Col #4~#6: 3D coordinates (X,Y,Z) for each bin, unit: nm.



### Running SnapFiSH

The shellscript is included in this github page. The **required** input variables are:

1. SnapFISH_dir : The directory of SnapFISH
2. input_dir : The directory of input file
3. output_dir : The directory for output
4. ann_file : The directory of file containing genomic location for each bin.
5. save_pic : 0 for outputting heatmaps, and 1 for not
6. data_name : Name of data for prefix

To run SnapFISH,
```
./run_SnapFISH.sh
```


### Contact Us
For any questions regarding this software, contact Ming Hu (hum@ccf.org) or Lindsay Lee (leeh7@ccf.org).
