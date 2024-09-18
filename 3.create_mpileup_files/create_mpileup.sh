#! /bin/bash
set -e
set -u
set -o pipefail

mkdir ~/sex_ratio_EE/results/3.create_mpileup_files
    # 4.1 create bam list for mpileup
    BAM_LIST=~/sex_ratio_EE/results/3.create_mpileup_files/bam_list_for_mpileup.txt

    # 4.2 extract names of fastq files:
    for SAMPLE in $(cat ${SAMPLE_LIST}); do
        echo "~/sex_ratio_EE/results/2.trim_map_markDup_reads/${SAMPLE}/${SAMPLE}.dupmarked.AmbigRm.bam" 
    done > ${BAM_LIST}

    # 4.3 create mpileup file for all samples
    REF=~/sex_ratio_EE/data/Rhrob_anchored.fa # reference genome

    samtools mpileup -b ${BAM_LIST} -B -Q 0 -f ${REF} -o ~/sex_ratio_EE/results/3.create_mpileup_files/EE_SR.mpileup # gzip later

