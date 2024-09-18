#!/bin/bash
set -e
set -u
set -o pipefail

### 7) GLM ON ALLELE COUNTS ###

    # the aim is to find loci which diverged between 4 male- and 4 female-biased lines consistently among replicates

# 7.1 split the sync file into smaller chunks- the input is too big to load it into R
  mkdir -p ~/sex_ratio_EE/results/7.glm_afc/split_files_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_allChr_filteredCov_tabSep
  cd ~/sex_ratio_EE/results/7.glm_afc/split_files_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_allChr_filteredCov_tabSep
  # split into chunks containing 10k lines each
  split -l 10000 ~/sex_ratio_EE/results/5.create_sync_file/EE_SR_IndelsRm_RepeatsRm_sexMergedCov_allChr_filteredCov.sync EE_SR_IndelsRm_RepeatsRm_sexMergedCov_allChr_filteredCov_segment_

# 7.2 remove monomorphic positions and keep only positions where coverage in all lines is between 50-202X (or 37-150X for the sex chromosome)
R_script=~/sex_ratio_EE/bin/7.glm_afc/filter_SNPs_keep_polymorphic.R
for file in ~/sex_ratio_EE/results/7.glm_afc/split_files_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_allChr_filteredCov_tabSep/*
 do
        Rscript --vanilla $R_script $file
 done > ~/sex_ratio_EE/results/7.glm_afc/polymorphicSnps_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_allChr_filteredCov.txt

# 7.3 prepare input for GLM
 # remove remaining headers from the SNP file (contig contig_pos SNPID line  A  T  C  G coverage MajAlleleCount_line selection)
  cd ~/sex_ratio_EE/results/7.glm_afc
  grep -v "coverage" polymorphicSnps_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_allChr_filteredCov.txt > polymorphicSnps_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_allChr_filteredCov.txt_noheader
  mv polymorphicSnps_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_allChr_filteredCov.txt_noheader polymorphicSnps_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_allChr_filteredCov.txt

# split again the input file for GLM to pieces containign 10k lines:
  mkdir ~/sex_ratio_EE/results/7.glm_afc/split_files_GLM_input
  cd ~/sex_ratio_EE/results/7.glm_afc/split_files_GLM_input
  split -l 10000 ../polymorphicSnps_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_allChr_filteredCov.txt polymorphicSnps_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_allChr_filteredCov_segment_

    # 7.4 run GLM on allele counts
R_script_GLM=~/sex_ratio_EE/bin/7.glm_afc/GLM.R
for file in ~/sex_ratio_EE/results/7.glm_afc/split_files_GLM_input/*
 do
        Rscript --vanilla R_script_GLM $file
 done > ~/sex_ratio_EE/results/7.glm_afc/GLM_output_quasibinomial.txt

    # 7.5 calculate q-values and draw Manhattan plots
Rsript ~/sex_ratio_EE/bin/7.glm_afc/4.calculate_qvalues_ManhattanPlots.R
