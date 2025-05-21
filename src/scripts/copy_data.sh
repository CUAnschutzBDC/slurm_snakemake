#!/bin/bash 


# Use this script to make soft links if the lane info isn't included
# Change the directories below to fit your directory scruture
original_dir=/pl/active/Anschutz_BDC/data/sussel/raw_data/20250416_LH00407_0129_B22YWC2LT3
new_dir=/pl/active/Anschutz_BDC/analysis/wells/analysis/sussel/maddy/acetate_RNA_seq/raw_data

cd $original_dir
for f in *_R1_001.fastq.gz *_R2_001.fastq.gz; do
  target=${original_dir}/$(readlink "$f")
  newname="${new_dir}/${f/_R/_L001_R}"
  ln -s "$target" "$newname"
done