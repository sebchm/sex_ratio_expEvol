# based on the F28 data, remove SNPs which are monomorphic (MAF < 0.05), and keep only positions with coverage between 53 and 215. 
# colons in the sync file have to be replaced with a field separator, e.g. sed "s/:/\\t/g" prior to loading into R and split into smaller files 
# the output is in the longer format

# author: Sebastian Chmielewski, Evolutionary Biology Group, AMU; sebchm_a_amu.edu.pl

# run it with a command:
# for file in <dir with split files>/*
#   do
# Rscript --vanilla <this script> $file
# done > polymorphicSnps_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_allChr_filteredCov.txt

maxcov <- 215 # upper coverage threshold for autosomes
mincov <- 53 # lower coverage threshold for autosomes
maxcov_sex <- 147 # upper coverage threshold for the sex chromosome
mincov_sex <- 36 # lower coverage threshold for the sex chromosome
minorAF <- 0.05 # minor allele frequency threshold

# suppressWarnings(library(dplyr)) # load libraries and suppress warnings
suppressWarnings(library(tidyr))

args = commandArgs(trailingOnly=TRUE)
infile=args[1]

AFCdata = read.table(infile, header=FALSE, sep='\t')
# test
#AFCdata = read.table("~/sex_ratio_EE/results/7.glm_afc/split_files_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_allChr_filteredCov_tabSep/polymorphicSnps_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_allChr_filteredCov_segment_aa", header=FALSE, sep='\t')

# assign column names
  #sample_names <- as.factor(unique(c("FB-A_F1", "FB-A_F1", "FB-B_F1", "FB-B_F1", "FB-C_F1", "FB-C_F1"  , "FB-D_F1", "FB-D_F1", "MB-A_F1", "MB-A_F1", "MB-B_F1", "MB-B_F1", "MB-C_F1", "MB-C_F1", "MB-D_F1", "MB-D_F1", "FB-A_F28", "FB-A_F28", "FB-B_F28", "FB-B_F28", "FB-C_F28", "FB-C_F28" , "FB-D_F28", "FB-D_F28", "MB-A_F28", "MB-A_F28", "MB-B_F28", "MB-B_F28", "MB-C_F28", "MB-C_F28", "MB-D_F28", "MB-D_F28")))
  #col_names_in_each_sample <- c("_A", "_T", "_C", "_G", "_N", "_d") 
  #col_names_all_samples <- paste(rep(sample_names, each = 6), col_names_in_each_sample, sep = "")
  #colnames(AFCdata) <- c("contig", "contig_pos", "ref_allele", col_names_all_samples)
colnames(AFCdata) <- c("contig", "contig_pos", "ref_allele", "FB-A_F1_A", "FB-A_F1_T", "FB-A_F1_C", "FB-A_F1_G", "FB-A_F1_N", "FB-A_F1_d", "FB-B_F1_A", "FB-B_F1_T", "FB-B_F1_C", "FB-B_F1_G", "FB-B_F1_N", "FB-B_F1_d", "FB-C_F1_A", "FB-C_F1_T", "FB-C_F1_C", "FB-C_F1_G", "FB-C_F1_N", "FB-C_F1_d", "FB-D_F1_A", "FB-D_F1_T", "FB-D_F1_C", "FB-D_F1_G", "FB-D_F1_N", "FB-D_F1_d", "MB-A_F1_A", "MB-A_F1_T", "MB-A_F1_C", "MB-A_F1_G", "MB-A_F1_N", "MB-A_F1_d", "MB-B_F1_A", "MB-B_F1_T", "MB-B_F1_C", "MB-B_F1_G", "MB-B_F1_N", "MB-B_F1_d", "MB-C_F1_A", "MB-C_F1_T", "MB-C_F1_C", "MB-C_F1_G", "MB-C_F1_N", "MB-C_F1_d", "MB-D_F1_A", "MB-D_F1_T", "MB-D_F1_C", "MB-D_F1_G", "MB-D_F1_N", "MB-D_F1_d", "FB-A_F28_A", "FB-A_F28_T", "FB-A_F28_C", "FB-A_F28_G", "FB-A_F28_N", "FB-A_F28_d", "FB-B_F28_A", "FB-B_F28_T", "FB-B_F28_C", "FB-B_F28_G", "FB-B_F28_N", "FB-B_F28_d", "FB-C_F28_A", "FB-C_F28_T", "FB-C_F28_C", "FB-C_F28_G", "FB-C_F28_N", "FB-C_F28_d", "FB-D_F28_A", "FB-D_F28_T", "FB-D_F28_C", "FB-D_F28_G", "FB-D_F28_N", "FB-D_F28_d", "MB-A_F28_A", "MB-A_F28_T", "MB-A_F28_C", "MB-A_F28_G", "MB-A_F28_N", "MB-A_F28_d", "MB-B_F28_A", "MB-B_F28_T", "MB-B_F28_C", "MB-B_F28_G", "MB-B_F28_N", "MB-B_F28_d", "MB-C_F28_A", "MB-C_F28_T", "MB-C_F28_C", "MB-C_F28_G", "MB-C_F28_N", "MB-C_F28_d", "MB-D_F28_A", "MB-D_F28_T", "MB-D_F28_C", "MB-D_F28_G", "MB-D_F28_N", "MB-D_F28_d")
  
AFCdata_AF_AC <- AFCdata[,-c(4:51,56,57,62,63,68,69,74,75,80,81,86,87,92,93,98,99)] # remove columns with F1, counts of deletions and Ns

AFCdata_AF_AC$sum_A = rowSums(AFCdata_AF_AC[,c(4,8,12,16,20,24,28,32)]) # get total A coverage
AFCdata_AF_AC$sum_T = rowSums(AFCdata_AF_AC[,c(5,9,13,17,21,25,29,33)])
AFCdata_AF_AC$sum_C = rowSums(AFCdata_AF_AC[,c(6,10,14,18,22,26,30,34)])
AFCdata_AF_AC$sum_G = rowSums(AFCdata_AF_AC[,c(7,11,15,19,23,27,31,35)])
AFCdata_AF_AC$F28.Total.cov = AFCdata_AF_AC$sum_A + AFCdata_AF_AC$sum_T + AFCdata_AF_AC$sum_C + AFCdata_AF_AC$sum_G # get total coverage
AFCdata_AF_AC$majAlleleCount = pmax(AFCdata_AF_AC$sum_A, AFCdata_AF_AC$sum_T, AFCdata_AF_AC$sum_C, AFCdata_AF_AC$sum_G) # get major allele count
AFCdata_AF_AC$majAllele = apply(AFCdata_AF_AC[, 36:39], 1, function(x) gsub("sum_","",names(x)[which.max(x)])) # get major allele
AFCdata_AF_AC$F28MAF = AFCdata_AF_AC$majAlleleCount / AFCdata_AF_AC$F28.Total.cov # calculate major allele frequency
AFCdata_AF_AC$F28Minorcount = AFCdata_AF_AC$F28.Total.cov - AFCdata_AF_AC$majAlleleCount # get minor allele count
AFCdata_AF_AC$F28MinAF = AFCdata_AF_AC$F28Minorcount/AFCdata_AF_AC$F28.Total.cov # calculate minor allele frequency
AFCdata_AF_AC$FB_A_cov_F28 = rowSums(AFCdata_AF_AC[,c(4,5,6,7)]) # calculate sum of coverage of female-biased A
AFCdata_AF_AC$FB_B_cov_F28 = rowSums(AFCdata_AF_AC[,c(8,9,10,11)])
AFCdata_AF_AC$FB_C_cov_F28 = rowSums(AFCdata_AF_AC[,c(12,13,14,15)])
AFCdata_AF_AC$FB_D_cov_F28 = rowSums(AFCdata_AF_AC[,c(16,17,18,19)])
AFCdata_AF_AC$MB_A_cov_F28 = rowSums(AFCdata_AF_AC[,c(20,21,22,23)])
AFCdata_AF_AC$MB_B_cov_F28 = rowSums(AFCdata_AF_AC[,c(24,25,26,27)])
AFCdata_AF_AC$MB_C_cov_F28 = rowSums(AFCdata_AF_AC[,c(28,29,30,31)])
AFCdata_AF_AC$MB_D_cov_F28 = rowSums(AFCdata_AF_AC[,c(32,33,34,35)])

# coverage filtering
AFCdata_AF_AC <- subset(AFCdata_AF_AC, contig != "LG4" & # keep only SNPs with coverage in all lines higher than a mincov
                      FB_A_cov_F28 >= mincov & FB_B_cov_F28 >= mincov & 
                      FB_C_cov_F28 >= mincov & FB_D_cov_F28 >= mincov & 
                      MB_A_cov_F28 >= mincov & MB_B_cov_F28 >= mincov & 
                      MB_C_cov_F28 >= mincov & MB_D_cov_F28 >= mincov &
                      FB_A_cov_F28 <= maxcov & FB_B_cov_F28 <= maxcov & # ... and lower than maxcov
                      FB_C_cov_F28 <= maxcov & FB_D_cov_F28 <= maxcov & 
                      MB_A_cov_F28 <= maxcov & MB_B_cov_F28 <= maxcov & 
                      MB_C_cov_F28 <= maxcov & MB_D_cov_F28 <= maxcov |
                      contig == "LG4" & 
                      FB_A_cov_F28 >= mincov_sex & FB_B_cov_F28 >= mincov_sex &  # for sex chromosome, use diferent coverage tresholds
                      FB_C_cov_F28 >= mincov_sex & FB_D_cov_F28 >= mincov_sex & 
                      MB_A_cov_F28 >= mincov_sex & MB_B_cov_F28 >= mincov_sex & 
                      MB_C_cov_F28 >= mincov_sex & MB_D_cov_F28 >= mincov_sex &
                      FB_A_cov_F28 <= maxcov_sex & FB_B_cov_F28 <= maxcov_sex &
                      FB_C_cov_F28 <= maxcov_sex & FB_D_cov_F28 <= maxcov_sex & 
                      MB_A_cov_F28 <= maxcov_sex & MB_B_cov_F28 <= maxcov_sex & 
                      MB_C_cov_F28 <= maxcov_sex & MB_D_cov_F28 <= maxcov_sex)
AFCdata_AF_AC <- subset(AFCdata_AF_AC, F28MinAF >= minorAF) # keep only SNPs with MinAF >= 0.05

AFCdata_AF_AC_longer <- AFCdata_AF_AC[,c(1:35, 42)] %>%
  tidyr::pivot_longer(cols = contains("_F28_"),
               names_sep = "_F28_",
               names_to = c("line", "allele"))
AFCdata_AF_AC_longer$SNPID = paste(AFCdata_AF_AC_longer$contig, AFCdata_AF_AC_longer$contig_pos, sep = "_")
AFCdata_AF_AC_longer <- AFCdata_AF_AC_longer %>%
  tidyr::pivot_wider(id_cols = c(contig, contig_pos, SNPID, line, majAllele),
              names_from = "allele",
              values_from = "value")

AFCdata_AF_AC_longer$coverage = AFCdata_AF_AC_longer$A + AFCdata_AF_AC_longer$T + AFCdata_AF_AC_longer$C + AFCdata_AF_AC_longer$G
AFCdata_AF_AC_longer$selection <- ifelse(grepl("^F", AFCdata_AF_AC_longer$line), "female_biased", "male_biased")
AFCdata_AF_AC_longer$MajAlleleCount_line <- apply(AFCdata_AF_AC_longer, 1, function(x) {
  x[as.character(x["majAllele"])]
})

options(max.print=1000000, width = 10000 )

if (nrow(AFCdata_AF_AC_longer) > 0) {
  print(as.data.frame(AFCdata_AF_AC_longer, col.names = NA), row.names=FALSE)
}

#### OLD VERSION BELOW. It's written using dplyr which makes it suuuper slow. 
# # modify the data
# AFCdata_AF_AC <- AFCdata %>%
#    select(!contains("F1")) %>% # remove first generation
#    select(!ends_with(c("_d", "_N"))) %>% # remove counts of deletions and Ns
#    mutate(sum_A = rowSums(select(., ends_with("F28_A"))), # sum all As
#           sum_T = rowSums(select(., ends_with("F28_T"))),
#           sum_C = rowSums(select(., ends_with("F28_C"))),
#           sum_G = rowSums(select(., ends_with("F28_G")))) %>%
#    mutate(F28.Total.cov = rowSums(select(., starts_with("sum_"))), # calculate total coverage#
#           majAlleleCount = pmax(sum_A, sum_T, sum_C, sum_G)) %>%
#    rowwise() %>%
#    mutate(majAllele = gsub("sum_", "", names(.)[36:39][which.max(c_across(c(sum_A, sum_T, sum_C, sum_G)))])) %>% # get major allele
#    ungroup() %>%
#    mutate(F28MAF = majAlleleCount/F28.Total.cov, # calculate major allele frequency
#             F28Minorcount = F28.Total.cov - majAlleleCount, # get minor allele count
#             F28MinAF = F28Minorcount/F28.Total.cov, # calculate minor allele frequency
#             FB_A_cov_F28 = rowSums(select(., starts_with("FB-A"))), # coverage sum of female-biased-A
#             FB_B_cov_F28 = rowSums(select(., starts_with("FB-B"))),
#             FB_C_cov_F28 = rowSums(select(., starts_with("FB-C"))),
#             FB_D_cov_F28 = rowSums(select(., starts_with("FB-D"))),
#             MB_A_cov_F28 = rowSums(select(., starts_with("MB-A"))), # coverage sum of male-biased-A
#             MB_B_cov_F28 = rowSums(select(., starts_with("MB-B"))),
#             MB_C_cov_F28 = rowSums(select(., starts_with("MB-C"))),
#             MB_D_cov_F28 = rowSums(select(., starts_with("MB-D")))) %>%
#    filter((contig == "LG4" & if_all(ends_with("_cov_F28"), ~ . >= mincov_sex & . <= maxcov_sex)) | #sex chromosome filtering- # keep only SNPs with coverage in all lines between mincov and maxcov
#             (contig != "LG4" & if_all(ends_with("_cov_F28"), ~ . >= mincov & . <= maxcov))) %>% # autosomal coverage filtering
#    filter(F28MinAF >= minorAF) # keep only SNPs with MinAF >= 0.05
# 
# 
#  ### change format to longer and keep only columns used for the GLM
# 
#  AFCdata_AF_AC_longer <- AFCdata_AF_AC %>%
#    select(c(1:35)) %>%
#    mutate(SNPID = paste(contig, contig_pos, sep = "_")) %>%
#    pivot_longer(cols = contains("_F28_"),
#                 names_sep = "_F28_",
#                 names_to = c("line", "allele")) %>%
#    pivot_wider(id_cols = c(contig, contig_pos, SNPID, line),
#                names_from = "allele",
#                values_from = "value") %>%
#    mutate(coverage = A + T + C + G,
#           MajAlleleCount_line = pmax(A, T, C, G),
#           selection = if_else(stringr::str_detect(line, "^F"), "female_biased", "male_biased"))

# options(max.print=1000000, width = 10000 )
# 
# if (nrow(AFCdata_AF_AC_longer) > 0) {
#   print(as.data.frame(AFCdata_AF_AC_longer, col.names = NA), row.names=FALSE)
# }

