### calculate q-values from GLM and create manhattan plots
  # browseVignettes("qvalue")

library(qvalue)
library(dplyr)
library(ggplot2)

setwd("~/sex_ratio_EE/results/7.glm_afc/")
GLM_output <- read.table("GLM_output_quasibinomial.txt", h = F, col.names = c("SNPID.ID", "CHR", "SNPpos", "pvalues"))

#GLM_output$SNPexactPos <- paste(GLM_output$V4, GLM_output$V5)
#GLM_output <- GLM_output[,-c(4,5)]

hist(GLM_output$pvalues)
qobj <- qvalue(p = GLM_output$pvalues)

summary(qobj)
pi0 <- qobj$pi0
localFDR <- qobj$lfdr

plot(qobj)
hist(qobj)
GLM_output$qvalues <- qobj$qvalues



# remove contigs starting with 'sq' (these were not incorporated to the chromosome-scale assembly) and rename chromosomes
GLM_output <- GLM_output  %>%
  filter(str_detect(CHR,"^LG")) %>%
  mutate(CHR = recode(CHR,
                     "LG7" = "chr1",
                     "LG5" = "chr2",
                     "LG3" = "chr3",
                     "LG1" = "chr4",
                     "LG8" = "chr5",
                     "LG4" = "chr6",
                     "LG2" = "chr7",
                     "LG6" = "chr8")) %>%
  mutate(CHR = factor(CHR, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8"))) %>% 
  arrange(CHR)

ggplot(GLM_output, aes(SNPpos, -log10(pvalues))) +
  geom_point() +
  facet_wrap(~CHR, scales = "free_x") +
  geom_hline(yintercept = -log10(6.209532e-08))

GLM_output %>%
  filter(pvalues < 0.05) %>%
  ggplot(aes(SNPpos, -(qvalues))) +
  geom_point() +
  facet_wrap(~CHR, scales = "free_x") +
  geom_hline(yintercept = -0.05, linetype = "dashed")

# analyse 7 SNPs with the lowest q-value. Only seven, because many SNPs have 8th lowest q-value 
SNPs_lowest_Qval <- GLM_output %>% slice_min(qvalues, n = 7)
  # save the most differentiated SNPs:
  # write.table(SNPs_lowest_Qval[1:5,1], "~/sex_ratio_EE/results/7.glm_afc/SNPs_lowest_qvalue.txt", col.names = F, row.names = F, quote = F)

  # get coverage data for differentiated SNPs
    # cd ~/sex_ratio_EE/results/7.glm_af
    # grep -w -f SNPs_lowest_qvalue.txt polymorphicSnps_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_allChr_filteredCov.txt > SNPs_lowest_q-value_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_allChr_filteredCov.txt
diff_SNPs <- read.table("~/sex_ratio_EE/results/7.glm_afc/SNPs_lowest_q-value_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_allChr_filteredCov.txt", col.names = c("contig", "contig_pos", "SNPID.ID", "line", "ref","A", "T", "C", "G", "coverage", "selection", "MajAlleleCount"))
diff_SNPs <- diff_SNPs %>%
  left_join(SNPs_lowest_Qval[,c(1,5)], by = "SNPID.ID") %>% # add q-value
  pivot_longer(cols = c(A,T,C,G), 
               names_to = "allele", 
               values_to = "n_reads") %>%
  mutate(qval_SNPID.ID = paste(round(qvalues, 2), SNPID.ID, sep = "\n"))

ggplot(diff_SNPs, aes(line, n_reads, fill = allele)) +
  geom_col() +
  facet_grid(qval_SNPID.ID ~ selection, scales = "free_x") +
  theme(legend.position = "top")
