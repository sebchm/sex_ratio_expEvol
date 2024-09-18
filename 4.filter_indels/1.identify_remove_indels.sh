#!bin/bash
set -e
set -u
set -o pipefail

    # 4.1 identify indels
perl ~/software/popoolation_1.2.2/basic-pipeline/identify-genomic-indel-regions.pl \
      --indel-window 5 \
      --min-count 5 \
      --input ~/sex_ratio_EE/results/3.create_mpileup_files/EE_SR.mpileup \
      --output ~/sex_ratio_EE/results/4.filter_indels/EE_SR_indels.gtf

    # 4.2 remove indels
perl ~/software/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl \
      --input ~/sex_ratio_EE/results/3.create_mpileup_files/EE_SR.mpileup \
      --gtf ~/sex_ratio_EE/results/4.filter_indels/EE_SR_indels.gtf \
      --output ~/sex_ratio_EE/results/4.filter_indels/EE_SR_IndelsRm.mpileup
      
gzip ~/sex_ratio_EE/results/3.create_mpileup_files/EE_SR.mpileup

    # remove repeat sequences (satellites, transposons etc.), based on repeat annotation from Chmielewski et al. 2024
perl ~/software/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl \
      --input ~/sex_ratio_EE/results/4.filter_indels/EE_SR_IndelsRm.mpileup \
      --gtf ~/sex_ratio_EE/results/4.filter_indels/AMU_Rhrob_2024_RepeatModeller_LG.gff \
      --output ~/sex_ratio_EE/results/4.filter_indels/EE_SR_IndelsRm_RepeatsRm.mpileup
