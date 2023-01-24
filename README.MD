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

For the 25Kb version, one package is requried in an addition to above

1. futures 3.0.5


### Details

- Chromatin tracing data from Huang et al, NG, 2021 (PMID: 34002095)
- 5Kb bin resolution, 41 loci cover a 205Kb region in mESC Sox2 locus. chr3:34,601,078–34,806,078 (ref: mm10).
- Cast allele contains the 7.5Kb 4CBS insertion.
- 129 allele does not contain the 7.5Kb 4CBS insertion. 
- After filtering, we keep 649 cells (649 129 alleles and 649 CAST alleles). On average, each allele contains 28 5Kb bins.
- SnapFISH outputs the list of summit. If you choose to output heatmaps, the PDF file will be located in /dir/to/SnapFISHoutput/output/tempfile. 

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

* SnapFish.py       
  * Step 1a: take 3D coordinates of each locus in each cell as input, calculate Eucildean disance between any two loci pairs, and also calculate average Euclidean disance for each loci pair across all available cells to create distance matrix.          
  * Step 1b: make heatmap of average Euclidean distance                       
  * Step 2: perform two sample T-test for each loci pair based on local background model               
  * Step 3a: identify loop candidates        
  * Step 3b: make heatmap of loop candidates               
  * Step 4a: group nearby loop candidates into clusters. remove singletons, only keep summits              
  * Step 4b: make heatmap of loop summits

### Running SnapFiSH

The shellscript is included in this github page. The required input variables are:

1. SnapFISH_dir : The directory of SnapFISH
2. input_dir : The directory of input file
3. output_dir : The directory for output
4. ann_file : The directory of file containing genomic location for each 5Kb bin. Look under "063022_input_ann.txt" of File Details section for more details 
5. save_pic : 0 for outputting heatmaps, and 1 for not
6. data_name : Name of data
7. bin_size : The resolution of bin size

For 25Kb version, the required input variables are:

1. SnapFISH_dir : The directory of SnapFISH for 25Kb
2. rep1_dir : The directory of input replicate 1 file
3. rep2_dir : The directory of input replicate 2 file
3. output_dir : The directory for output
4. ann_file : The directory of file containing genomic location for each 25Kb bin.
5. num_CPUs : The number of CPUs one would like to use. One can check how many CPUs are available by "lscpu". If num_CPUs = 1, the program will run as a single processor. When using a HPC with job scheduler, make sure to ask for 1 node.


### Contact Us
For any questions regarding this software, contact Ming Hu (hum@ccf.org) or Lindsay Lee (leeh7@ccf.org).
