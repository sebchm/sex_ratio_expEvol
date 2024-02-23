#!/bin/bash
set -e
set -u
set -o pipefail

SAMPLE_LIST=(...)/sample_list.txt
FASTQ_DIR=(...)
for IND in $(cat ${SAMPLE_LIST})
    do
    source (...)/rr_trim_map.sh ${IND} *${IND}*R1*.fastq.gz *${IND}*R2*.fastq.gz
    done
