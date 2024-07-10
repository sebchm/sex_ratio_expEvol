
### run adapted CMH-test using ACER package
  # female-biased lines (diff in the Ne compared to male-biased lines)
# author: Sebastian Chmielewski, Evolutionary Biology Group, AMU; sebchm_a_amu.edu.pl

  # install poolSeq package
    # install.packages("~/software/poolseq_R_package/v0.3.5.tar.gz", repos=NULL, type="source")
    # install.packages("~/software/ACER_R_package/ACER_1.0.3.tar.gz", repos=NULL, type="source")

suppressWarnings(library(ACER))
suppressWarnings(library(poolSeq))

args = commandArgs(trailingOnly=TRUE)
infile=args[1]

mySync <- read.sync(infile, gen=c(rep(1,4), rep(28,4)), repl=rep(c(1,2,3,4),2),polarization = "minor")

# validate the sync 
  # sync_t <- read.table(file="~/sex_ratio_EE/results/6.filter_by_coverage/sync_after_filtering/female-biased_lines/split_files/autosomes_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_filteredCovMAF_femaleBiasedLines_segment_aa")
  # colnames(sync_t) <- c("CHR", "POS", "REF", unique(c("FB-A_F1", "FB-A_F1", "FB-B_F1", "FB-B_F1", "FB-C_F1", "FB-C_F1", "FB-D_F1", "FB-D_F1", "FB-A_F28", "FB-B_F28", "FB-B_F28", "FB-C_F28", "FB-C_F28" , "FB-D_F28", "FB-D_F28")))
  # colnames(sync_t) <- c("CHR", "POS", "REF", unique(c("FB-A_F1", "FB-A_F1", "FB-B_F1", "FB-B_F1", "FB-C_F1", "FB-C_F1"  , "FB-D_F1", "FB-D_F1", "MB-A_F1", "MB-A_F1", "MB-B_F1", "MB-B_F1", "MB-C_F1", "MB-C_F1", "MB-D_F1", "MB-D_F1", "FB-A_F28", "FB-A_F28", "FB-B_F28", "FB-B_F28", "FB-C_F28", "FB-C_F28" , "FB-D_F28", "FB-D_F28", "MB-A_F28", "MB-A_F28", "MB-B_F28", "MB-B_F28", "MB-C_F28", "MB-C_F28", "MB-D_F28", "MB-D_F28")))

# allele frequency matrix
af <- af(mySync, repl = rep(c(1,2,3,4),2), gen = c(rep(1,4), rep(28,4)))

# coverage matrix
cov <- coverage(mySync, repl = rep(c(1,2,3,4),2), gen = c(rep(1,4), rep(28,4)))

# load effective population size
Ne <- read.table("~/sex_ratio_EE/results/10.estimate_Ne/Ne_est_all_lines.txt", h = T)

pval_ACER <- adapted.cmh.test(freq = af, 
                           coverage = cov, 
                           gen = c(1,28),
                           repl = c(1,2,3,4), 
                           Ne=c(Ne[1,2], Ne[2,2], Ne[3,2], Ne[4,2]),
                           poolSize = rep(200, 8),
                           IntGen = FALSE, 
                           order = 0, 
                           mincov = 53, 
                           RetVal = 0, 
                           correct = TRUE)

SNP_pvalue <- data.frame(SNP = row.names(af), 
                         pvalue = pval_ACER)

options(max.print=1000000, width = 10000)


if (nrow(SNP_pvalue) > 0) {
  print(SNP_pvalue, row.names=FALSE)
}
