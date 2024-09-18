#!/bin/bash
set -e
set -u
set -o pipefail

# filter the sync file based on the coverage. Low coverage regions are more likely to contain repeats, while high coverage indicate copy number variation
mkdir ~/sex_ratio_EE/results/6.filter_by_coverage
  
    # 6.1 merge coverage of male and female samples
awk -f ~/sex_ratio_EE/bin/5.create_sync_file/merge_samples.awk \
  ~/sex_ratio_EE/results/5.create_sync_file/EE_SR_IndelsRm_RepeatsRm.sync > \ # sync file
  ~/sex_ratio_EE/results/5.create_sync_file/EE_SR_IndelsRm_RepeatsRm_sexMergedCov.sync

    # 6.2 split autosomes and sex chromosome
        # sex chromosome is expected to have lower coverage, as males are hemizygous for X (bulb mites have X0 system). Sex chromosome is LG4 (=chr6)
            # create sync with sex chromosome positions
grep "^LG4" ~/sex_ratio_EE/results/5.create_sync_file/EE_SR_IndelsRm_RepeatsRm_sexMergedCov.sync > ~/sex_ratio_EE/results/5.create_sync_file/sexChr/EE_SR_IndelsRm_RepeatsRm_sexMergedCov_sexChromosome.sync
            # create sync with autosomes 
grep -v "^LG4" ~/sex_ratio_EE/results/5.create_sync_file/EE_SR_IndelsRm_RepeatsRm_sexMergedCov.sync > ~/sex_ratio_EE/results/5.create_sync_file/autosomes/EE_SR_IndelsRm_RepeatsRm_sexMergedCov_autosomes.sync

    # 6.3 subsample the sync file after merging F and M coverage 
    # whole coverage filtering is based on the concept of the 'target coverage' --> see manuscript. Here, I generate data to to visualise the distribution of the coverage in order to obtain the target coverage
        # 6.3.1 before merging female and male coverage
          #  6.3.1.1 sex chromosome (select random 0.6% lines)
grep "^LG4" ~/sex_ratio_EE/results/5.create_sync_file/EE_SR_IndelsRm_RepeatsRm.sync | awk 'BEGIN {srand()} !/^$/ { if (rand() <= .0006) print $0}' > ~/sex_ratio_EE/results/5.create_sync_file/sexChr/subsampled_files/EE_SR_IndelsRm_RepeatsRm_sexChromosome_subsampled.sync
          #  6.3.1.2 autosomes
grep -v "^LG4" ~/sex_ratio_EE/results/5.create_sync_file/EE_SR_IndelsRm_RepeatsRm | awk 'BEGIN {srand()} !/^$/ { if (rand() <= .0006) print $0}' > ~/sex_ratio_EE/results/5.create_sync_file/autosomes/subsampled_files/EE_SR_IndelsRm_RepeatsRm_autosomes_subsampled.sync
        
        # 6.3.2 after merging female and male coverage
          # 6.3.2.1 sex chromosome
awk 'BEGIN {srand()} !/^$/ { if (rand() <= .0006) print $0}' ~/sex_ratio_EE/results/5.create_sync_file/sexChr/EE_SR_IndelsRm_RepeatsRm_sexMergedCov_sexChromosome.sync > ~/sex_ratio_EE/results/5.create_sync_file/sexChr/subsampled_files/EE_SR_IndelsRm_RepeatsRm_sexMergedCov_sexChromosome_subsampled.sync
          # 6.3.2.2 autosomes
awk 'BEGIN {srand()} !/^$/ { if (rand() <= .0006) print $0}' ~/sex_ratio_EE/results/5.create_sync_file/autosomes/EE_SR_IndelsRm_RepeatsRm_sexMergedCov_autosomes.sync > ~/sex_ratio_EE/results/5.create_sync_file/autosomes/subsampled_files/EE_SR_IndelsRm_RepeatsRm_sexMergedCov_autosomes_subsampled.sync

    # 6.4 run the Rmd script to draw distribution of coverage and determine the taget coverage. Coverage filtering threshold are based on this target (=peak) coverage
    # ~/sex_ratio_EE/results/5.create_sync_file/draw_distribution_coverage.Rmd

    # 6.5 filter the sync files based on coverage
python ~/sex_ratio_EE/bin/6.filter_by_coverage/FilterPositionOnCoverage.py /media/raid/home/schmielewski/sex_ratio_EE/results/5.create_sync_file/EE_SR_IndelsRm_sexMergedCov_Autosomes.sync 51 202 > /media/raid/home/schmielewski/sex_ratio_EE/results/6.filter_by_coverage/EE_SR_IndelsRm_sexMergedCov_Autosomes_FilteredCov.sync
python ~/sex_ratio_EE/bin/6.filter_by_coverage/FilterPositionOnCoverage.py /media/raid/home/schmielewski/sex_ratio_EE/results/5.create_sync_file/EE_SR_IndelsRm_RepeatsRm_sexMergedCov_sexChromosome.sync 38 150 > /media/raid/home/schmielewski/sex_ratio_EE/results/6.filter_by_coverage/EE_SR_IndelsRm_sexMergedCov_sexChromosome_FilteredCov.sync

    # 6.6 draw coverage distribution after filtering on coverage
        # autosomes
awk 'BEGIN {srand()} !/^$/ { if (rand() <= .0006) print $0}' /media/raid/home/schmielewski/sex_ratio_EE/results/5.create_sync_file/autosomes/EE_SR_IndelsRm_RepeatsRm_sexMergedCov_autosomes_filteredCov.sync > /media/raid/home/schmielewski/sex_ratio_EE/results/5.create_sync_file/autosomes/EE_SR_IndelsRm_RepeatsRm_sexMergedCov_autosomes_filteredCov_subsampled.sync
        # sex chromosome
awk 'BEGIN {srand()} !/^$/ { if (rand() <= .0001) print $0}' /media/raid/home/schmielewski/sex_ratio_EE/results/5.create_sync_file/sexChr/EE_SR_IndelsRm_RepeatsRm_sexMergedCov_sexChromosome_filteredCov.sync > /media/raid/home/schmielewski/sex_ratio_EE/results/5.create_sync_file/sexChr/EE_SR_IndelsRm_RepeatsRm_sexMergedCov_sexChromosome_filteredCov_subsampled.sync

    # ~/sex_ratio_EE/results/5.create_sync_file/draw_distribution_coverage.Rmd

    # 6.7 merge autosomes and sex chromosome
cat /media/raid/home/schmielewski/sex_ratio_EE/results/5.create_sync_file/autosomes/EE_SR_IndelsRm_RepeatsRm_sexMergedCov_autosomes_filteredCov_subsampled.sync /media/raid/home/schmielewski/sex_ratio_EE/results/5.create_sync_file/sexChr/EE_SR_IndelsRm_RepeatsRm_sexMergedCov_sexChromosome_filteredCov.sync >  ~/sex_ratio_EE/results/5.create_sync_file/EE_SR_IndelsRm_RepeatsRm_sexMergedCov_allChr_filteredCov.sync # gzip later (!)
