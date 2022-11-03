# SnapFISH: single nucleus analysis pipeline for super-resolution imaging data



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
- 5Kb bin resolution, 41 loci cover a 205Kb region in mESC Sox2 locus. chr3:34,601,078–34,806,078 (ref: mm10).
- Cast allele contains the 7.5Kb 4CBS insertion.
- 129 allele does not contain the 7.5Kb 4CBS insertion. 
- after filtering, we keep 649 cells (649 129 alleles and 649 CAST alleles). On average, each allele contains 28 5Kb bins.

### File Details
The input files are located in the ext folder.

* 063022_input_ann.txt 
  * input data, genomic location for each 5Kb bin. Col #4 "pos" is the unique ID for each 5Kb bin.

* 063022_input_3D_coordinates_129.txt
  * input data for the 129 allele. 
  * Col #1: unique ID for each cell
  * Col #2: unique ID for each 5Kb bin
  * Col #3~#5: 3D coordinates (X,Y,Z) for each 5Kb bin, unit: nm.

* 063022_input_3D_coordinates_CAST.txt    
  * input data for the CAST allele. 

* SnapFish_step1a.py       
  * take 3D coordinates of each locus in each cell as input, calculate Eucildean disance between any two loci pairs, and also calculate average Euclidean disance for each loci pair across all available cells to create distance matrix.

* SnapFiSH_step1b.py            
  * make heatmap of average Euclidean distance

* SnapFiSH_step2.py                         
  * perform two sample T-test for each loci pair based on local background model

* SnapFiSH_step3a.py                
  * identify loop candidates

* SnapFiSH_step3b.py            
  * make heatmap of loop candidates

* SnapFiSH_step4a.py                    
  * group nearby loop candidates into clusters. remove singletons, only keep summits

* SnapFiSH_step4b.py                
  * make heatmap of loop summits

### Running SnapFiSH

The shellscript is included in this github page. The required inputs variables are:

1. SnapFISH_dir : The directory of SnapFISH
2. input_dir : The directory of input file
3. output_dir : The directory for output
4. ann_file : The directory of file containing genomic location for each 5Kb bin. Look under "063022_input_ann.txt" of File Details section for more details 
5. save_pic : 0 for outputting heatmaps, and 1 for not
6. data_name : Name of data

### Contact Us
For any questions regarding this software, contact Ming Hu (hum@ccf.org) or Lindsay Lee (leeh7@ccf.org).
