#!/bin/bash
set -e
set -u
set -o pipefail
# run qualimap for bams after removing AmbReads and marking duplicates
QMAP=~/software/qualimap_v2.2.1/qualimap
QUALIMAP_AMBRM_DUPMARKED_BAMS=~/sex_ratio_EE/results/2_A.qualimap/bam_list_AfterAmbReadsRm_DupMarked

# remove bam lists if exist:
if [ -f "$QUALIMAP_AMBRM_DUPMARKED_BAMS" ]; then
    rm "$QUALIMAP_AMBRM_DUPMARKED_BAMS"
fi

# prepare files with sample names and corresponding bam locations:
for i in $(cat ~/sex_ratio_EE/data/sample_list.txt)
        do
        echo -e "${i}\t/media/raid/home/schmielewski/sex_ratio_EE/results/2.trim_map_markDup_reads/${i}/${i}.dupmarked.AmbigRm.bam" >> ${QUALIMAP_AMBRM_DUPMARKED_BAMS}
        done

# run qualimap
${QMAP} multi-bamqc -r -d ${QUALIMAP_AMBRM_DUPMARKED_BAMS} -outformat HTML --java-mem-size=20G -outdir ~/sex_ratio_EE/results/2_A.qualimap/ambRm_dumparked_bams
