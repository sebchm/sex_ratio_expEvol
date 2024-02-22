#!/bin/bash
set -e
set -u
set -o pipefail

SAMPLE_LIST=~/sex_ratio_EE/data/sample_list.txt
FASTQ_DIR=/media/nas/pool_GWAS
for IND in $(cat ${SAMPLE_LIST})
    do
    source ~/sex_ratio_EE/bin/2.trim_map_markDup_reads/rr_trim_map.sh ${IND} *${IND}*R1*.fastq.gz *${IND}*R2*.fastq.gz
    done
