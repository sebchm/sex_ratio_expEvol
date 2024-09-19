
#!/bin/bash
set -e
set -u
set -o pipefail


#####################
### 11) ACER TEST ###
#####################
    # 11.1 separate male-biased from female-biased lines
PATH_BASE=~/sex_ratio_EE/results/6.filter_by_coverage/sync_after_filtering
mkdir -p ${PATH_BASE}/female-biased_lines/split_files
mkdir -p ${PATH_BASE}/male-biased_lines/split_files

        # 11.1.1 AUTOSOMES: separate male-biased from female-biased lines
       # get only female-biased lines from the sync file
cut -f 1,2,3,4,5,6,7,12,13,14,15 ${PATH_BASE}/autosomes/autosomes_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_allChr_filteredCovMAF.sync > ${PATH_BASE}/female-biased_lines/autosomes_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_filteredCovMAF_femaleBiasedLines.sync

       # get only male-biased lines from the sync file
cut -f 1,2,3,8,9,10,11,16,17,18,19 ${PATH_BASE}/autosomes/autosomes_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_allChr_filteredCovMAF.sync > ${PATH_BASE}/male-biased_lines/autosomes_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_filteredCovMAF_maleBiasedLines.sync

        # 11.1.2 SEX CHROMOSOME: separate male-biased from female-biased lines
        # get only female-biased lines from the sync file
cut -f 1,2,3,4,5,6,7,12,13,14,15 ${PATH_BASE}/sexChr/sexChr_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_allChr_filteredCovMAF.sync > ${PATH_BASE}/female-biased_lines/sexChr_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_filteredCovMAF_femaleBiasedLines.sync

       # get only male-biased lines from the sync file
cut -f 1,2,3,8,9,10,11,16,17,18,19 ${PATH_BASE}/sexChr/sexChr_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_allChr_filteredCovMAF.sync > ${PATH_BASE}/male-biased_lines/sexChr_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_filteredCovMAF_maleBiasedLines.sync

    # 11.2 split the files into chunks containing only 10k lines
       # 11.2.1 AUTOSOMES: split the input file to pieces containing 10k lines
cd ${PATH_BASE}/female-biased_lines/split_files
split -l 10000 ../autosomes_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_filteredCovMAF_femaleBiasedLines.sync autosomes_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_filteredCovMAF_femaleBiasedLines_segment_

       cd ${PATH_BASE}/male-biased_lines/split_files
split -l 10000 ../autosomes_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_filteredCovMAF_maleBiasedLines.sync autosomes_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_filteredCovMAF_maleBiasedLines_segment_

        # 11.2.2 SEX CHROMOSOME: split the input file to pieces containing 10k lines
cd ${PATH_BASE}/female-biased_lines/split_files
split -l 10000 ../sexChr_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_filteredCovMAF_femaleBiasedLines.sync sexChr_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_filteredCovMAF_femaleBiasedLines_segment_

cd ${PATH_BASE}/male-biased_lines/split_files
split -l 10000 ../sexChr_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_filteredCovMAF_maleBiasedLines.sync sexChr_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_filteredCovMAF_maleBiasedLines_segment_

       # 11.3 ACER test: female autosomes
R_script=~/sex_ratio_EE/bin/11.acer_test/afterSyncFiltering/acer_F-biased_lines.R

for file in ~/sex_ratio_EE/results/6.filter_by_coverage/sync_after_filtering/female-biased_lines/split_files/autosomes/*
 do
        Rscript --vanilla $R_script $file
 done > ~/sex_ratio_EE/results/11.acer_test/ACER_output_femaleBiasedLines.txt

        # 11.3 ACER test: male autosomes
R_script=~/sex_ratio_EE/bin/11.acer_test/afterSyncFiltering/acer_M-biased_lines.R

for file in ~/sex_ratio_EE/results/6.filter_by_coverage/sync_after_filtering/male-biased_lines/split_files/autosomes/*
 do
        Rscript --vanilla $R_script $file
 done > ~/sex_ratio_EE/results/11.acer_test/ACER_output_maleBiasedLines.txt

        # 11.4 ACER test: female sex chromosome
 R_script=~/sex_ratio_EE/bin/11.acer_test/afterSyncFiltering/acer_F-biased_lines_sexChr.R

for file in ~/sex_ratio_EE/results/6.filter_by_coverage/sync_after_filtering/female-biased_lines/split_files/sexChr/*
 do
        Rscript --vanilla $R_script $file
 done > ~/sex_ratio_EE/results/11.acer_test/ACER_output_femaleBiasedLines_sexChr.txt

        # 11.5 ACER test: male sex chromosome
 R_script=~/sex_ratio_EE/bin/11.acer_test/afterSyncFiltering/acer_M-biased_lines_sexChr.R

for file in ~/sex_ratio_EE/results/6.filter_by_coverage/sync_after_filtering/male-biased_lines/split_files/sexChr/*
 do
        Rscript --vanilla $R_script $file
 done > ~/sex_ratio_EE/results/11.acer_test/ACER_output_maleBiasedLines_sexChr.txt
