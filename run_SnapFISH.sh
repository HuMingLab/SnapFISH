
#!/bin/bash

############################################################################
###                            User Variables                            ###
############################################################################


input_dir="/home/input/directory/129_3d_coor.txt"
output_dir="/home/output/directory"
ann_file="/home/input/directory/input_ann.txt"
save_pic=1      # 0 for outputting heatmaps, and 1 for not
data_name="129"

############################################################################
if [ ! -d "$output_dir" ]; then
  mkdir $output_dir
fi

python3 SnapFISH.py -i $input_dir -o $output_dir -a $ann_file -p $save_pic -d $data_name
