#!/bin/bash
set -e
set -u
set -o pipefail

SAMPLE_LIST=(...)/sample_list.txt
FASTQ_DIR=(...)
for IND in $(cat ${SAMPLE_LIST})
    do
    source rr_trim_map.sh ${IND} VA-3199-${IND}_*R1_001.fastq.gz VA-3199-${IND}_*R1_001.fastq.gz
    done
