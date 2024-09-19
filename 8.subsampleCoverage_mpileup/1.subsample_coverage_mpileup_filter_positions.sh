##################################################
### 8) SUBSAMPLE MPILEUP TO A UNIFORM COVERAGE ###
##################################################
    # mpileup files with subsampled coverage will be used to estimate Ne and genome-wide diversity measures (pi, theta, D)
PATH_MPILEUP=~/sex_ratio_EE/results/4.filter_indels/mpileup_separateLines_sexMerged_repeatsRm_indelsRm
mkdir ${PATH_MPILEUP}/temp_folder
mkdir ~/sex_ratio_EE/results/8.subsampleCoverage_mpileup

    # 8.1 subsample to uniform coverage
for line in $(cat ~/sex_ratio_EE/data/line_names.txt)
        do
                # unzip the mpileup file since PoPoolation doesn't accept zcat...
gunzip -c ${PATH_MPILEUP}/EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_${line}.mpileup.gz > ${PATH_MPILEUP}/temp_folder/EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_${line}.mpileup

        # run PoPoolation
perl ~/software/popoolation_1.2.2/basic-pipeline/subsample-pileup.pl \
        --min-qual 20 \
        --method withoutreplace \
        --max-coverage 215 \
        --fastq-type sanger \
        --target-coverage 53 \
        --input ${PATH_MPILEUP}/temp_folder/EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_${line}.mpileup \
        --output ~/sex_ratio_EE/results/8.subsampleCoverage_mpileup/EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_${line}.mpileup

        # gzip the output
gzip ~/sex_ratio_EE/results/8.subsampleCoverage_mpileup/EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_${line}.mpileup

        # remove unzipped mpileup
rm ${PATH_MPILEUP}/temp_folder/EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_${line}.mpileup
        done

    # 8.2 keep only positions common in all samples
        # positions with coverage < target-coverage are removed. Due to variance in coverage between samples, the number of positions in the output is different for all samples. To ensure comparability between samples, only SNPs present in all mpileup will be retained
mkdir mpileup_files_SS_the_same_snps
        # 8.2.1 get chr and pos from each subsampled mpileup
> ~/sex_ratio_EE/results/8.subsampleCoverage_mpileup/chr_pos_line_SS_mpileups.txt

for i in $(cat ~/sex_ratio_EE/data/line_names.txt)
do
    cut -f1,2 ~/sex_ratio_EE/results/8.subsampleCoverage_mpileup/EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_${i}.autosomes.mpileup | sed "s/$/\t${i}/" >> ~/sex_ratio_EE/results/8.subsampleCoverage_mpileup/chr_pos_line_SS_mpileups.txt
done

        # 8.2.2 make SNPID column
awk '{print $1, $2, $1"_"$2}' chr_pos_line_SS_mpileups.txt > chr_pos_SNPID_line_SS_mpileups.txt

        # 8.2.3 count number of occurences of each SNP - I expect 16 if a SNP is present in all 16 samples
cut -f3 chr_pos_SNPID_line_SS_mpileups.txt | sort | uniq -c > counts_chr_pos_SNPID_line_SS_mpileups.txt

        # 8.2.4 count occurences of each snp
awk '$1 == 16 {print $2 "\t" $3}' counts_chr_pos_SNPID_line_SS_mpileups.txt > common_SNPs_SS_mpileups.txt

        # 8.2.5 keep only SNPs with SNPs present in all 16 samples
for i in ~/sex_ratio_EE/results/8.subsampleCoverage_mpileup/*mpileup; do
    awk 'NR==FNR { snps[$1 FS $2]; next } ($1 FS $2) in snps' common_SNPs_SS_mpileups.txt ${i} > mpileup_files_SS_the_same_snps/"${file%.mpileup}_commonSNPs.mpileup"
done
