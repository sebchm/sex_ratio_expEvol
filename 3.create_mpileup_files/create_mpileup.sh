#! /bin/bash
set -e
set -u
set -o pipefail

# create bam list for mpileup
SAMPLE_LIST=(...)/sex_ratio_EE/data/sample_list.txt
BAM_LIST=(...)/3.create_mpileup_files/bam_list_for_mpileup.txt

# extract names of fastq files:
for SAMPLE in $(cat ${SAMPLE_LIST}); do
    echo "(...)/2.trim_map_markDup_reads/${SAMPLE}/${SAMPLE}.dupmarked.AmbigRm.bam" 
done > ${BAM_LIST}

# create mpileup file for all samples
REF=(...)/sex_ratio_EE/data/Rhrob_anchored.fa

samtools mpileup -b ${BAM_LIST} -B -Q 0 -f ${REF} -o (...)/EE_SR.mpileup
gzip (...)/EE_SR.mpileup
