### calculate q values and create manhattan plots
  # browseVignettes("qvalue")
# author: Sebastian Chmielewski, Evolutionary Biology Group, AMU; sebchm_a_amu.edu.pl

library(qvalue)
library(dplyr)
library(ggplot2)

setwd("~/sex_ratio_EE/results/7.glm_afc/")
GLM_output <- read.table("GLM_output_quasibinomial.txt", h = F, col.names = c("SNPID.ID", "CHR", "SNPpos", "pvalues"))

hist(GLM_output$pvalues)
qobj <- qvalue(p = GLM_output$pvalues)

summary(qobj)
pi0 <- qobj$pi0
localFDR <- qobj$lfdr

plot(qobj)
hist(qobj)
GLM_output$qvalues <- qobj$qvalues

# plot q-values VS p-values:
GLM_output %>% 
  filter(pvalues <= 0.05) %>%
  ggplot(aes(pvalues, qvalues)) +
  geom_point()

# remove sqs and rename chromosomes
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

ggplot(GLM_output, aes(SNPpos, -p_adjusted_BH)) +
  geom_point() +
  facet_wrap(~CHR, scales = "free_x") +
  geom_hline(yintercept = -0.05)
