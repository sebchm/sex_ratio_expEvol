#!/bin/bash
set -e
set -u
set -o pipefail

QMAP=~/software/qualimap_v2.2.1/qualimap
# lists with sample name and bam location:
QUALIMAP_INPUT_RAW=(...)/bam_list_rawMappedReads
QUALIMAP_INPUT_AMBRM_DUPMARKED=(...)/bam_list_AfterAmbReadsRm_DupMarked

# remove input lists if exist:
if [ -f "$QUALIMAP_INPUT_RAW" ]; then
    rm "$QUALIMAP_INPUT_RAW"
fi

if [ -f "$QUALIMAP_INPUT_AMBRM_DUPMARKED" ]; then
    rm "$QUALIMAP_INPUT_AMBRM_DUPMARKED"
fi

# prepare files with sample names and corresponding bam locations:
for i in $(cat ~/sex_ratio_EE/data/sample_list.txt)
        do
        #samtools sort  (...)/${i}/${i}.out_sorted.bam > (...)/2.trim_map_markDup_reads/${i}/${i}.out_resorted.bam
        echo -e "${i}\t/(...)/2.trim_map_markDup_reads/${i}/${i}.out_sorted.bam" >> ${QUALIMAP_INPUT_RAW}
        echo -e "${i}\t/(...)/2.trim_map_markDup_reads/${i}/${i}.dupmarked.AmbigRm.bam" >> ${QUALIMAP_INPUT_AMBRM_DUPMARKED}
        done

# prepare qualimap reports:
${QMAP} multi-bamqc -r -d ${QUALIMAP_INPUT_RAW} -outformat HTML --java-mem-size=10G
${QMAP} multi-bamqc -r -d ${QUALIMAP_INPUT_AMBRM_DUPMARKED} -outformat HTML --java-mem-size=10G
