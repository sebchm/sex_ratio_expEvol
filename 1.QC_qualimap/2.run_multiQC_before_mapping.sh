#!bin/bash

REPORT_DIR=~/sex_ratio_EE/results/1.multiQC
SAMPLE_LIST=~/sex_ratio_EE/data/sample_list.txt

# activate conda if it's not activated

# prepare report for each sample
for SAMPLE in $(cat ${SAMPLE_LIST})
        do
fastqc /media/nas/EE_SR/VA-3199-${SAMPLE}* -o ${REPORT_DIR}
        done

# prepare report for all samples
multiqc ${REPORT_DIR}/ -o ${REPORT_DIR}
