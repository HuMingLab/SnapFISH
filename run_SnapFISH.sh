
#!/bin/bash

############################################################################
###                            User Variables                            ###
############################################################################


SnapFISH_dir="/home/SnapFiSH/python/"
input_dir="/home/SnapFiSH/python/063022_input_3D_coordinates_129.txt"
output_dir="/home/SnapFiSH/python/output"
ann_file="/home/SnapFiSH/python/063022_input_ann.txt"
save_pic=1
data_name="129"
bin_size=5000

########## For 25Kb version ########

rep1_dir="/home/SnapFiSH/101522_mESC_DNAseqFISH+_25Kb_autosomes_268718loci_Rep1.txt"                                                               
rep2_dir="/home/SnapFiSH//101522_mESC_DNAseqFISH+_25Kb_autosomes_394347loci_Rep2.txt" 
num_CPUs=10


############################################################################
if [ ! -d "$output_dir" ]; then
  mkdir $output_dir
fi

if [ ! -d "${output_dir}/tempfile" ]; then
  mkdir ${output_dir}/tempfile
fi


if [[ $bin_size -eq 25000 ]]; then
python ${SnapFISH_dir}/SnapFISH25/SnapFISH25Kb.py -s $SnapFISH_dir  -o $output_dir -a $ann_file -d $data_name --rep1 $rep1_dir --rep2 $rep2_dir -n $num_CPUs                     
fi

if [[ $bin_size -ne 25000 ]]; then
python ${SnapFISH_dir}/SnapFISH.py -s $SnapFISH_dir -i $input_dir -o $output_dir -a $ann_file -p $save_pic -d $data_name                                                
fi


