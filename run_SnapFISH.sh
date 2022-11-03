
#!/bin/bash

############################################################################
###                            User Variables                            ###
############################################################################


SnapFISH_dir="/home/leeh7/SnapFiSH/python/"
input_dir="/home/leeh7/SnapFiSH/python/063022_input_3D_coordinates_129.txt"
output_dir="/home/leeh7/SnapFiSH/python/output"
ann_file="/home/leeh7/SnapFiSH/python/063022_input_ann.txt"
save_pic=0
data_name="129"


############################################################################
if [ ! -d "$output_dir" ]; then
  mkdir $output_dir
fi

if [ ! -d "${output_dir}/tempfile" ]; then
  mkdir ${output_dir}/tempfile
fi

python ${SnapFISH_dir}/SnapFISH.py -s $SnapFISH_dir -i $input_dir -o $output_dir -a $ann_file -p $save_pic -d $data_name

