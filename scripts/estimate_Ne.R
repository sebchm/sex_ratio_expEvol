
### estimate Ne based on subsampled coverage files
# author: Sebastian Chmielewski, Evolutionary Biology Group, AMU; sebchm_a_amu.edu.pl

library(dplyr)
library(ggplot2)
library(poolSeq)
# Estimate Ne

#### PART 1: autosomes ####

# female-biased A
FB_A_sync <- read.sync("~/sex_ratio_EE/results/8.subsampleCoverage_mpileup/mpileup_files_SS_the_same_snps/autosomes/sync_files/filteredSameSNPsAcrossLines_EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_F-A_F1_F28.autosomes.sync",gen=c(1,28),repl=c(1,1),polarization = "minor")
AF <- af(FB_A_sync,repl=c(1,1),gen=c(1,28))
Cov <- coverage(FB_A_sync,repl=c(1,1),gen=c(1,28))
FB_A_Ne <- estimateNe(p0=AF[,1],pt=AF[,2],cov0=Cov[,1],covt=Cov[,2],t=27, method="P.planII",poolSize = c(200,200),truncAF = 0.01)

# female-biased B
FB_B_sync <- read.sync("~/sex_ratio_EE/results/8.subsampleCoverage_mpileup/mpileup_files_SS_the_same_snps/autosomes/sync_files/filteredSameSNPsAcrossLines_EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_F-B_F1_F28.autosomes.sync",gen=c(1,28),repl=c(1,1),polarization = "minor")
AF <- af(FB_B_sync,repl=c(1,1),gen=c(1,28))
Cov <- coverage(FB_B_sync,repl=c(1,1),gen=c(1,28))
FB_B_Ne <- estimateNe(p0=AF[,1],pt=AF[,2],cov0=Cov[,1],covt=Cov[,2],t=27, method="P.planII",poolSize = c(200,200),truncAF = 0.01)

# female-biased C
FB_C_sync <- read.sync("~/sex_ratio_EE/results/8.subsampleCoverage_mpileup/mpileup_files_SS_the_same_snps/autosomes/sync_files/filteredSameSNPsAcrossLines_EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_F-C_F1_F28.autosomes.sync",gen=c(1,28),repl=c(1,1),polarization = "minor")
AF <- af(FB_C_sync,repl=c(1,1),gen=c(1,28))
Cov <- coverage(FB_C_sync,repl=c(1,1),gen=c(1,28))
FB_C_Ne <- estimateNe(p0=AF[,1],pt=AF[,2],cov0=Cov[,1],covt=Cov[,2],t=27, method="P.planII",poolSize = c(200,200),truncAF = 0.01)

# female-biased D
FB_D_sync <- read.sync("~/sex_ratio_EE/results/8.subsampleCoverage_mpileup/mpileup_files_SS_the_same_snps/autosomes/sync_files/filteredSameSNPsAcrossLines_EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_F-D_F1_F28.autosomes.sync",gen=c(1,28),repl=c(1,1),polarization = "minor")
AF <- af(FB_D_sync,repl=c(1,1),gen=c(1,28))
Cov <- coverage(FB_D_sync,repl=c(1,1),gen=c(1,28))
FB_D_Ne <- estimateNe(p0=AF[,1],pt=AF[,2],cov0=Cov[,1],covt=Cov[,2],t=27, method="P.planII",poolSize = c(200,200),truncAF = 0.01)

# male-biased A
MB_A_sync <- read.sync("~/sex_ratio_EE/results/8.subsampleCoverage_mpileup/mpileup_files_SS_the_same_snps/autosomes/sync_files/filteredSameSNPsAcrossLines_EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_M-A_F1_F28.autosomes.sync",gen=c(1,28),repl=c(1,1),polarization = "minor")
AF <- af(MB_A_sync,repl=c(1,1),gen=c(1,28))
Cov <- coverage(MB_A_sync,repl=c(1,1),gen=c(1,28))
MB_A_Ne <- estimateNe(p0=AF[,1],pt=AF[,2],cov0=Cov[,1],covt=Cov[,2],t=27, method="P.planII",poolSize = c(200,200),truncAF = 0.01)

# male-biased B
MB_B_sync <- read.sync("~/sex_ratio_EE/results/8.subsampleCoverage_mpileup/mpileup_files_SS_the_same_snps/autosomes/sync_files/filteredSameSNPsAcrossLines_EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_M-B_F1_F28.autosomes.sync",gen=c(1,28),repl=c(1,1),polarization = "minor")
AF <- af(MB_B_sync,repl=c(1,1),gen=c(1,28))
Cov <- coverage(MB_B_sync,repl=c(1,1),gen=c(1,28))
MB_B_Ne <- estimateNe(p0=AF[,1],pt=AF[,2],cov0=Cov[,1],covt=Cov[,2],t=27, method="P.planII",poolSize = c(200,200),truncAF = 0.01)

# male-biased C
MB_C_sync <- read.sync("~/sex_ratio_EE/results/8.subsampleCoverage_mpileup/mpileup_files_SS_the_same_snps/autosomes/sync_files/filteredSameSNPsAcrossLines_EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_M-C_F1_F28.autosomes.sync",gen=c(1,28),repl=c(1,1),polarization = "minor")
AF <- af(MB_C_sync,repl=c(1,1),gen=c(1,28))
Cov <- coverage(MB_C_sync,repl=c(1,1),gen=c(1,28))
MB_C_Ne <- estimateNe(p0=AF[,1],pt=AF[,2],cov0=Cov[,1],covt=Cov[,2],t=27, method="P.planII",poolSize = c(200,200),truncAF = 0.01)

# male-biased D
MB_D_sync <- read.sync("~/sex_ratio_EE/results/8.subsampleCoverage_mpileup/mpileup_files_SS_the_same_snps/autosomes/sync_files/filteredSameSNPsAcrossLines_EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_M-D_F1_F28.autosomes.sync",gen=c(1,28),repl=c(1,1),polarization = "minor")
AF <- af(MB_D_sync,repl=c(1,1),gen=c(1,28))
Cov <- coverage(MB_D_sync,repl=c(1,1),gen=c(1,28))
MB_D_Ne <- estimateNe(p0=AF[,1],pt=AF[,2],cov0=Cov[,1],covt=Cov[,2],t=27, method="P.planII",poolSize = c(200,200),truncAF = 0.01)

Ne_all_lines <- data.frame(line = c("FB_A_Ne", "FB_B_Ne", "FB_C_Ne", "FB_D_Ne", "MB_A_Ne", "MB_B_Ne", "MB_C_Ne", "MB_D_Ne"), 
                           Ne = c(FB_A_Ne, FB_B_Ne, FB_C_Ne, FB_D_Ne, MB_A_Ne, MB_B_Ne, MB_C_Ne, MB_D_Ne))

Ne_all_lines <- Ne_all_lines %>%
  mutate(selection = ifelse(stringr::str_detect(line, "^FB"), "female_biased", "male_biased"),
         type = "autosomes")

# write.table(Ne_all_lines, "~/sex_ratio_EE/results/10.estimate_Ne/Ne_est_all_lines_subSampledCoverage.txt", sep = "\t", row.names = F, quote = F)
Ne_all_lines <- read.table("~/sex_ratio_EE/results/10.estimate_Ne/Ne_est_all_lines_subSampledCoverage.txt", h = T)

#### PART 2: SEX CHROMOSOME ####
# poolSize has been reduced from 200 to 150, as each pool comprises 100 females and 100 males, but males carry only half of the female's sex chromosomes
# female-biased A
FB_A_sync <- read.sync("~/sex_ratio_EE/results/8.subsampleCoverage_mpileup/mpileup_files_SS_the_same_snps/sex_chr/sync_files/EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_F-A_F1_F28_commonSNPsOnly_sex_chr.sync",gen=c(1,28),repl=c(1,1),polarization = "minor")
AF <- af(FB_A_sync,repl=c(1,1),gen=c(1,28))
Cov <- coverage(FB_A_sync,repl=c(1,1),gen=c(1,28))
FB_A_Ne_sexChr <- estimateNe(p0=AF[,1],pt=AF[,2],cov0=Cov[,1],covt=Cov[,2],t=27, method="P.planII",poolSize = c(150,150),truncAF = 0.01)

# female-biased B
FB_B_sync <- read.sync("~/sex_ratio_EE/results/8.subsampleCoverage_mpileup/mpileup_files_SS_the_same_snps/sex_chr/sync_files/EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_F-B_F1_F28_commonSNPsOnly_sex_chr.sync",gen=c(1,28),repl=c(1,1),polarization = "minor")
AF <- af(FB_B_sync,repl=c(1,1),gen=c(1,28))
Cov <- coverage(FB_B_sync,repl=c(1,1),gen=c(1,28))
FB_B_Ne_sexChr <- estimateNe(p0=AF[,1],pt=AF[,2],cov0=Cov[,1],covt=Cov[,2],t=27, method="P.planII",poolSize = c(150,150),truncAF = 0.01)

# female-biased C
FB_C_sync <- read.sync("~/sex_ratio_EE/results/8.subsampleCoverage_mpileup/mpileup_files_SS_the_same_snps/sex_chr/sync_files/EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_F-C_F1_F28_commonSNPsOnly_sex_chr.sync",gen=c(1,28),repl=c(1,1),polarization = "minor")
AF <- af(FB_C_sync,repl=c(1,1),gen=c(1,28))
Cov <- coverage(FB_C_sync,repl=c(1,1),gen=c(1,28))
FB_C_Ne_sexChr <- estimateNe(p0=AF[,1],pt=AF[,2],cov0=Cov[,1],covt=Cov[,2],t=27, method="P.planII",poolSize = c(150,150),truncAF = 0.01)

# female-biased D
FB_D_sync <- read.sync("~/sex_ratio_EE/results/8.subsampleCoverage_mpileup/mpileup_files_SS_the_same_snps/sex_chr/sync_files/EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_F-D_F1_F28_commonSNPsOnly_sex_chr.sync",gen=c(1,28),repl=c(1,1),polarization = "minor")
AF <- af(FB_D_sync,repl=c(1,1),gen=c(1,28))
Cov <- coverage(FB_D_sync,repl=c(1,1),gen=c(1,28))
FB_D_Ne_sexChr <- estimateNe(p0=AF[,1],pt=AF[,2],cov0=Cov[,1],covt=Cov[,2],t=27, method="P.planII",poolSize = c(150,150),truncAF = 0.01)

# male-biased A
MB_A_sync <- read.sync("~/sex_ratio_EE/results/8.subsampleCoverage_mpileup/mpileup_files_SS_the_same_snps/sex_chr/sync_files/EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_M-A_F1_F28_commonSNPsOnly_sex_chr.sync",gen=c(1,28),repl=c(1,1),polarization = "minor")
AF <- af(MB_A_sync,repl=c(1,1),gen=c(1,28))
Cov <- coverage(MB_A_sync,repl=c(1,1),gen=c(1,28))
MB_A_Ne_sexChr <- estimateNe(p0=AF[,1],pt=AF[,2],cov0=Cov[,1],covt=Cov[,2],t=27, method="P.planII",poolSize = c(150,150),truncAF = 0.01)

# male-biased B
MB_B_sync <- read.sync("~/sex_ratio_EE/results/8.subsampleCoverage_mpileup/mpileup_files_SS_the_same_snps/sex_chr/sync_files/EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_M-B_F1_F28_commonSNPsOnly_sex_chr.sync",gen=c(1,28),repl=c(1,1),polarization = "minor")
AF <- af(MB_B_sync,repl=c(1,1),gen=c(1,28))
Cov <- coverage(MB_B_sync,repl=c(1,1),gen=c(1,28))
MB_B_Ne_sexChr <- estimateNe(p0=AF[,1],pt=AF[,2],cov0=Cov[,1],covt=Cov[,2],t=27, method="P.planII",poolSize = c(150,150),truncAF = 0.01)

# male-biased C
MB_C_sync <- read.sync("~/sex_ratio_EE/results/8.subsampleCoverage_mpileup/mpileup_files_SS_the_same_snps/sex_chr/sync_files/EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_M-C_F1_F28_commonSNPsOnly_sex_chr.sync",gen=c(1,28),repl=c(1,1),polarization = "minor")
AF <- af(MB_C_sync,repl=c(1,1),gen=c(1,28))
Cov <- coverage(MB_C_sync,repl=c(1,1),gen=c(1,28))
MB_C_Ne_sexChr <- estimateNe(p0=AF[,1],pt=AF[,2],cov0=Cov[,1],covt=Cov[,2],t=27, method="P.planII",poolSize = c(150,150),truncAF = 0.01)

# male-biased D
MB_D_sync <- read.sync("~/sex_ratio_EE/results/8.subsampleCoverage_mpileup/mpileup_files_SS_the_same_snps/sex_chr/sync_files/EE_SR_IndelsRm_RepeatsRm_filteredCov_sexMerged_subSampledCov_M-D_F1_F28_commonSNPsOnly_sex_chr.sync",gen=c(1,28),repl=c(1,1),polarization = "minor")
AF <- af(MB_D_sync,repl=c(1,1),gen=c(1,28))
Cov <- coverage(MB_D_sync,repl=c(1,1),gen=c(1,28))
MB_D_Ne_sexChr <- estimateNe(p0=AF[,1],pt=AF[,2],cov0=Cov[,1],covt=Cov[,2],t=27, method="P.planII",poolSize = c(150,150),truncAF = 0.01)

Ne_sexChr_all_lines <- data.frame(line = c("FB_A_Ne_sexChr", "FB_B_Ne_sexChr", "FB_C_Ne_sexChr", "FB_D_Ne_sexChr", "MB_A_Ne_sexChr", "MB_B_Ne_sexChr", "MB_C_Ne_sexChr", "MB_D_Ne_sexChr"),
                                  Ne = c(FB_A_Ne_sexChr, FB_B_Ne_sexChr, FB_C_Ne_sexChr, FB_D_Ne_sexChr, MB_A_Ne_sexChr, MB_B_Ne_sexChr, MB_C_Ne_sexChr, MB_D_Ne_sexChr))

Ne_sexChr_all_lines <- Ne_sexChr_all_lines %>%
  mutate(selection = ifelse(stringr::str_detect(line, "^FB"), "female_biased", "male_biased"),
         type = "sex_chromosome")

#write.table(Ne_sexChr_all_lines, "~/sex_ratio_EE/results/10.estimate_Ne/subsampled_coverage/Ne_est_sexChr_all_lines.txt", sep = "\t", row.names = F, quote = F)
Ne_sexChr_all_lines <- read.table("~/sex_ratio_EE/results/10.estimate_Ne/subsampled_coverage/Ne_est_sexChr_all_lines.txt", h = T)

Ne <- rbind(Ne_all_lines, Ne_sexChr_all_lines)
Ne$line <- sub("_Ne$", "", Ne$line)
Ne$line <- sub("_Ne_sexChr$", "", Ne$line)


#### plotting

ne_mean <- Ne_all_lines %>%
  group_by(selection) %>%
  summarise(mean_ne_gen = mean(Ne),
            sd_ne = sd(Ne),
            n = n(),
            se_ne = sd_ne / sqrt(n)) %>%
  mutate(selection = ifelse(selection == "female_biased", "female-biased", "male-biased"))

Ne_autosomes <- Ne %>%
  filter(type == "autosomes") %>%
  mutate(selection = ifelse(selection == "female_biased", "female-biased", "male-biased"), 
         line = gsub("B_", "-", line))


ggplot(ne_mean, aes(selection, mean_ne_gen, color = selection)) +
  geom_pointrange(aes(ymin = mean_ne_gen - sd_ne, ymax = mean_ne_gen + sd_ne),
                  position = position_dodge(width = 0.5), linewidth = 1.4) +
  ggrepel::geom_label_repel(data = Ne_autosomes, aes(selection, Ne, label = line, color = selection),
                   position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5, seed = 123 ), size = 3, fill = NA) +
  geom_point(data = Ne_autosomes, aes(selection, Ne, color = selection),
             alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5, seed = 123 ), size = 3) +
  theme_bw() +
  theme(legend.position = "none", 
        legend.title = element_blank()) +
  scale_color_manual(values = c("#67AC55", "#F97A69")) +
  labs(y = "Ne") +
  geom_errorbar(aes(ymin = mean_ne_gen - sd_ne, ymax = mean_ne_gen + sd_ne),
                width = 0.05, position = position_dodge(width = 0.1)) 
