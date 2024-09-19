#!/bin/bash
set -e
set -u
set -o pipefail

##########################################
### 9) GENOME-WIDE DIVERSITY ESTIMATES ###
##########################################

mkdir ~/sex_ratio_EE/bin/9.genome_wide_diversity

    # 9.1 calculate pi, theta and Tajima's D at introns and exons
POPOOLATION_PATH=/media/raid/home/schmielewski/software/popoolation_1.2.2
RESULTS_DIR=/media/raid/home/schmielewski/sex_ratio_EE/results/9.genome_wide_diversity/exon_intron/mpileup_files_SS_the_same_snps

cd /media/raid/home/schmielewski/sex_ratio_EE/results/8.subsampleCoverage_mpileup/mpileup_files_SS_the_same_snps

for MEASURE in pi theta d
do
        for GEN in 1 28
        do
            for SELECTION in F M
            do
                for LINE in A B C D
                do
                    for SEQ_TYPE in intron exon
                    do

# GTF file- has been filtered to contain only genes significantly expresses in adults (from Plesnar-Bielak et al. 2024)
GTF_FILE=/media/raid/home/schmielewski/sex_ratio_EE/bin/9.genome_wide_diversity/GTF_files/Rhrob_anchored_${SEQ_TYPE}_filtered.gtf

# min-count has to be set to 2 when d is calculated
    if [ "$MEASURE" == "d" ]; then
        MIN_COUNT=2
    else
        MIN_COUNT=3
    fi

# Syn-nonsyn-at-position
perl ${POPOOLATION_PATH}/Variance-at-position.pl \
       --fastq-type sanger \
       --measure ${MEASURE} \
       --gtf ${GTF_FILE} \
       --pileup filteredSameSNPsAcrossLines_EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_${SELECTION}-${LINE}_F${GEN}.autosomes.mpileup \
       --output ${RESULTS_DIR}/filteredSameSNPsAcrossLines_${SELECTION}-${LINE}_F${GEN}.autosomes.SS.${MEASURE}_${SEQ_TYPE} \
       --pool-size 400 \
       --min-count ${MIN_COUNT} \
       --min-coverage 53 \
       --max-coverage 215 \
       --min-qual 20 \
       --snp-output ${RESULTS_DIR}/filteredSameSNPsAcrossLines_${SELECTION}-${LINE}_F${GEN}.SynVSNon.autosomes.${MEASURE}_${SEQ_TYPE}.SS.SNP
                done
            done
        done
    done
done


    # 9.2 calculate pi, theta and Tajima's D at synonymous and non-synonymous sites
POPOOLATION_PATH=/media/raid/home/schmielewski/software/popoolation_1.2.2
GTF_FILE=/media/raid/home/schmielewski/sex_ratio_EE/bin/9.genome_wide_diversity/GTF_files/Rhrob_anchored_CDS_filtered.gtf
RESULTS_DIR=/media/raid/home/schmielewski/sex_ratio_EE/results/9.genome_wide_diversity/syn_nonsym/mpileup_files_SS_the_same_snps

MPILEUP_PATH=/media/raid/home/schmielewski/sex_ratio_EE/results/8.subsampleCoverage_mpileup/mpileup_files_SS_the_same_snps

for MEASURE in pi theta d
do
    for GEN in 1 28
    do
        for SELECTION in F M
        do
            for LINE in A B C D
            do

# Syn-nonsyn-at-position
perl ${POPOOLATION_PATH}/syn-nonsyn/Syn-nonsyn-at-position.pl \
       --fastq-type sanger \
       --measure ${MEASURE} \
       --gtf ${GTF_FILE} \
       --pileup ${MPILEUP_PATH}/filteredSameSNPsAcrossLines_EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_${SELECTION}-${LINE}_F${GEN}.autosomes.mpileup \
       --codon-table ${POPOOLATION_PATH}/syn-nonsyn/codon-table.txt \
       --nonsyn-length-table ${POPOOLATION_PATH}/syn-nonsyn/nsl_p1.txt \
       --output ${RESULTS_DIR}/${MEASURE}/filteredSameSNPsAcrossLines_${SELECTION}-${LINE}_F${GEN}.SynVSNon.autosomes.SS.${MEASURE} \
       --pool-size 400 \
       --min-count 3 \
       --min-coverage 53 \
       --max-coverage 215 \
       --min-qual 20 \
       --snp-output ${RESULTS_DIR}/${MEASURE}/filteredSameSNPsAcrossLines_${SELECTION}-${LINE}_F${GEN}.SynVSNon.autosomes.${MEASURE}.SS.SNP \
       --region-output ${RESULTS_DIR}/${MEASURE}/filteredSameSNPsAcrossLines_${SELECTION}-${LINE}_F${GEN}.SynVSNon.autosomes.${MEASURE}.SS.Region
            done
        done
    done
done

    # 9.3 calculate pi, theta and Tajima's D in 10kb windows
RESULTS_DIR=/media/raid/home/schmielewski/sex_ratio_EE/results/9.genome_wide_diversity/windowed_analysis/mpileup_files_SS_the_same_snps

cd /media/raid/home/schmielewski/sex_ratio_EE/results/8.subsampleCoverage_mpileup/mpileup_files_SS_the_same_snps

for MEASURE in pi d theta
do
    for GEN in 1 28
    do
        for SELECTION in F M
        do
            for LINE in A B C D
            do

# min-count has to be set to 2 when d is calculated
                    if [ "$MEASURE" == "d" ]; then
                        MIN_COUNT=2
                    else
                        MIN_COUNT=3
                    fi

# Syn-nonsyn-at-position
 perl ~/software/popoolation_1.2.2/Variance-sliding.pl \
       --window-size    10000 \
       --step-size 10000 \
       --fastq-type sanger \
       --measure ${MEASURE} \
       --input ${MPILEUP}/filteredSameSNPsAcrossLines_EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_${SELECTION}-${LINE}_F${GEN}.autosomes.mpileup \
       --output ${RESULTS_DIR}/filteredSameSNPsAcrossLines_${SELECTION}-${LINE}_F${GEN}.10kbWindows.autosomes.SS.${MEASURE} \
       --pool-size 400 \
       --min-count ${MIN_COUNT} \
       --min-coverage 53 \
       --max-coverage 215 \
       --min-covered-fraction 0.2 \
       --min-qual 20 \
       --snp-output ${RESULTS_DIR}/filteredSameSNPsAcrossLines_${SELECTION}-${LINE}_F${GEN}.10kbWindows.autosomes.${MEASURE}.SS.SNP
            done
        done
    done
done
