#!bin/bash
set -e
set -u
set -o pipefail

SAMPLE_LIST=~/sex_ratio_EE/data/sample_list.txt

# extract names of fastq files:
for file in /.../EE_SR/*.fastq.gz; do
    name=$(basename "$file" | sed -E 's/TJ-2710-([^_]+).*/\1/')
    echo "$name" 
done | uniq > ${SAMPLE_LIST}
