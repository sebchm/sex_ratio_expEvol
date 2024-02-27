#!/bin/bash
set -e
set -u
set -o pipefail

# run qualimap for sorted bams
QMAP=~/software/qualimap_v2.2.1/qualimap
        # lists with sample name and bam location:
QUALIMAP_SORTED_BAMS=(...)/2_A.qualimap/bam_list_rawMappedReads

# remove bam lists if exist:
if [ -f "$QUALIMAP_SORTED_BAMS" ]; then
    rm "$QUALIMAP_SORTED_BAMS"
fi

# prepare files with sample names and corresponding bam locations:
for i in $(cat ~/sex_ratio_EE/data/sample_list.txt)
        do
        echo -e "${i}\t(...)/2.trim_map_markDup_reads/${i}/${i}_resorted.bam" >> ${QUALIMAP_SORTED_BAMS}
        done

# run qualimap
${QMAP} multi-bamqc -r -d ${QUALIMAP_SORTED_BAMS} -outformat HTML --java-mem-size=20G -outdir (...)/2_A.qualimap/sorted_bams
