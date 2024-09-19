#!/bin/bash
set -e
set -u
set -o pipefail
###############################################
### 10) EFFECTIVE POPULATION SIZE ESTIMATION ###
###############################################

    # 10.1 prepare sync files: each sync should contain only one line (F1 and evolved samples)

        # 10.1.1: autosomes: convert SS mpileups to sync

    # mpileup files have one generation only, and I need sync file with 2 generations (CHR, POS, REF, F1, F28 in columns)
    # so, firstly, I will convert all mpileup files into sync files and paste sync files.
PATH_MPILEUP=/media/raid/home/schmielewski/sex_ratio_EE/results/8.subsampleCoverage_mpileup/mpileup_files_SS_the_same_snps/autosomes
PATH_SYNC=/media/raid/home/schmielewski/sex_ratio_EE/results/8.subsampleCoverage_mpileup/mpileup_files_SS_the_same_snps/autosomes/sync_files
mkdir -p ${PATH_MPILEUP}
mkdir -p ${PATH_SYNC}

for i in $(cat ~/sex_ratio_EE/data/line_names.txt)
    do
PILEUP_INPUT=${PATH_MPILEUP}/filteredSameSNPsAcrossLines_EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_${i}.autosomes.mpileup
SYNC_OUTPUT=${PATH_SYNC}/filteredSameSNPsAcrossLines_EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_${i}.autosomes.sync

java -ea -Xmx60g -jar ~/software/popoolation2_1201/mpileup2sync.jar \
        --input ${PILEUP_INPUT} \
        --output ${SYNC_OUTPUT} \
        --fastq-type sanger \
        --min-qual 20 \
        --threads 40
    done

cd ${PATH_SYNC}
# append F28 to the sync file from F1:

for i in $(cat ~/sex_ratio_EE/data/line_names_noGen.txt)
    do
SYNC_F1=${PATH_SYNC}/filteredSameSNPsAcrossLines_EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_${i}_F1.autosomes.sync
SYNC_F28=${PATH_SYNC}/filteredSameSNPsAcrossLines_EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_${i}_F28.autosomes.sync

awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[NR]=$4; next} {print $0, a[FNR]}' ${SYNC_F28} ${SYNC_F1} > ${PATH_SYNC}/filteredSameSNPsAcrossLines_EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_${i}_F1_F28.autosomes.sync
    done

        # 10.1.2: sex chromosome: convert SS mpileups to sync

PATH_MPILEUP=/media/raid/home/schmielewski/sex_ratio_EE/results/8.subsampleCoverage_mpileup/mpileup_files_SS_the_same_snps/sex_chr
PATH_SYNC=/media/raid/home/schmielewski/sex_ratio_EE/results/8.subsampleCoverage_mpileup/mpileup_files_SS_the_same_snps/sex_chr/sync_files
mkdir -p ${PATH_MPILEUP}
mkdir -p ${PATH_SYNC}

for i in $(cat ~/sex_ratio_EE/data/line_names.txt)
    do

PILEUP_INPUT=${PATH_MPILEUP}/EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_${i}.mpileup.gz_commonSNPsOnly_sex_chr.mpileup
SYNC_OUTPUT=${PATH_SYNC}/EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_${i}_commonSNPsOnly_sex_chr.sync

java -ea -Xmx60g -jar ~/software/popoolation2_1201/mpileup2sync.jar \
        --input ${PILEUP_INPUT} \
        --output ${SYNC_OUTPUT} \
        --fastq-type sanger \
        --min-qual 20 \
        --threads 45
    done

cd ${PATH_SYNC}
# append F28 to the sync file from F1:

for i in $(cat ~/sex_ratio_EE/data/line_names_noGen.txt)
    do
SYNC_F1=${PATH_SYNC}/EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_${i}_F1_commonSNPsOnly_sex_chr.sync
SYNC_F28=${PATH_SYNC}/EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_${i}_F28_commonSNPsOnly_sex_chr.sync

awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[NR]=$4; next} {print $0, a[FNR]}' ${SYNC_F28} ${SYNC_F1} > ${PATH_SYNC}/EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_${i}_F1_F28_commonSNPsOnly_sex_chr.sync
    done

    # 10.2 use estimateNe from poolSeq package to Estimate Ne
Rscript estimate_Ne.R
