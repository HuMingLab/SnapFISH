
#!/bin/bash

############################################################################
###                            User Variables                            ###
############################################################################

SnapFISH_dir="/home/SnapFISH/directory/"
input_path="/home/input/directory/129_3d_coor.txt"
output_dir="/home/output/directory"
ann_file="/home/input/directory/input_ann.txt"
save_pic=0     # 1 for outputting heatmaps, and 0 for not
data_name="129"

############################################################################
if [ ! -d "$output_dir" ]; then
  mkdir $output_dir
fi

python3 $SnapFISH_dir/SnapFISH.py -i $input_path -o $output_dir -a $ann_file -p $save_pic -d $data_name
