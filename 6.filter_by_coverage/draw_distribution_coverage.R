
# coverage_check
# author: "S. Chmielewski"
# date: "2024-02-28"

library(ggplot2)
library(dplyr)
library(viridis)
library(plotly)
library(tidyr)
library(stringr)

# This scripts plots the coverage of 98k random positions, based on mpileup with indels removed. 
# Autosomal coverage

#  autosomes}
# cov is a sync file after removind indels
cov <- read.table("~/sex_ratio_EE/results/5.create_sync_file/autosomes/EE_SR_IndelsRm_RepeatsRm_autosomes_subsampled.sync")

# remove N and del (last 2 digits from each cell)
remove_last_digits <- function(x) {
  gsub(".{0,4}$", "", x)
}
cov[4:ncol(cov)] <- lapply(cov[4:ncol(cov)], remove_last_digits)

# add col names
sample_names <- read.table("~/sex_ratio_EE/data/sample_info.txt", sep = "\t", h = T)[,5]
colnames(cov) <- c("LG", "pos", "ref_allele", sample_names)

# sum values within each cell (sum the coverage of all nucleotides)
sum_split_values <- function(x) {
  sapply(strsplit(x, ":"), function(y) sum(as.numeric(y)))
}
  # apply the function
cov <- cov %>%
  mutate(across(contains("males"), sum_split_values))

# change the format to longer
cov_longer_autosomes <- cov %>%
  arrange(LG, pos) %>% # sort the mpileup according to position
  pivot_longer(cols = tidyselect::contains("male"), # change format to longer
    names_to = c("line", "sex", "generation"),
    names_sep = "_",
    values_to = "coverage")%>%
  mutate(line_generation = paste0(line, "_", generation, "_", sex)) 



### Calculate mean coverage per line

# calculate mean coverage per line:
cov_longer_autosomes %>%
  group_by(line, generation, sex) %>%
  summarise(mean_coverage = mean(coverage)) 


### Plot raw coverage of each line

# plot coverage
cov_plot <- ggplot(cov_longer_autosomes, aes(coverage)) +
  geom_density(aes(color = line_generation)) +
  geom_density(linewidth = 1.5) +
  xlim(0,200) +
  #ylim(0,0.055) +
  scale_color_viridis(discrete = T) +
  theme(legend.position = "none")+
  #facet_wrap(~generation) +
  ggtitle("Before summing coverage across sexes; autosomes only")

plot(cov_plot)

# convert to a plotly plot
interactive_cov_plot <- ggplotly(cov_plot)
interactive_cov_plot

# plot coverage (only female samples)
cov_longer_autosomes %>%
  filter(sex == "females") %>%
  ggplot(aes(coverage)) +
  geom_density(aes(color = line_generation)) +
  geom_density(linewidth = 1.5) +
  xlim(0,200) +
  #ylim(0,0.055) +
  scale_color_viridis(discrete = T) +
  theme(legend.position = "none")+
  #facet_wrap(~generation) +
  ggtitle("Before summing coverage across sexes; female autosomes only") 

## sum autosomal coverage of males and females
The target coverage for autosomes is 100.98 (vertical line)
Therefore, the mean coverage across all lines has to be between 50X and 201X.

#  after summing male and female coverage}
# cov_check is a sync file after merging male and female coverage
cov_merged <- read.table("~/sex_ratio_EE/results/5.create_sync_file/autosomes/EE_SR_IndelsRm_RepeatsRm_sexMergedCov_autosomes_subsampled.sync")

# remove N and del (last 2 digits from each cell)
cov_merged[4:ncol(cov_merged)] <- lapply(cov_merged[4:ncol(cov_merged)], remove_last_digits)

# add col names
sample_names_check <- unique(c("FB-A_F1", "FB-A_F1", "FB-B_F1", "FB-B_F1", "FB-C_F1", "FB-C_F1"  , "FB-D_F1", "FB-D_F1", "MB-A_F1", "MB-A_F1", "MB-B_F1", "MB-B_F1", "MB-C_F1", "MB-C_F1", "MB-D_F1", "MB-D_F1", "FB-A_F28", "FB-A_F28", "FB-B_F28", "FB-B_F28", "FB-C_F28", "FB-C_F28" , "FB-D_F28", "FB-D_F28", "MB-A_F28", "MB-A_F28", "MB-B_F28", "MB-B_F28", "MB-C_F28", "MB-C_F28", "MB-D_F28", "MB-D_F28"))
colnames(cov_merged) <- c("LG", "pos", "ref_allele", sample_names_check)

# sum values within each cell (sum the coverage of all nucleotides)
cov_merged <- cov_merged %>%
  mutate(across(contains("_F"), sum_split_values))

# change the format to longer
cov_merged_longer_autosomes <- cov_merged %>%
  arrange(LG, pos) %>% # sort the mpileup according to position
  pivot_longer(cols = tidyselect::contains("_F"), # change format to longer
    names_to = c("line", "generation"),
    names_sep = "_",
    values_to = "coverage")%>%
  mutate(line_generation = paste0(line, "_", generation)) 

cov_merged_longer_autosomes %>%
  group_by(line, generation) %>%
  summarise(mean_cov = mean(coverage))


# plot coverage
plot_cov_merged_longer_autosomes <- ggplot(cov_merged_longer_autosomes, aes(coverage)) +
  geom_density(aes(color = line_generation)) +
  geom_density(linewidth = 1.5) +
  xlim(0,400) +
  #ylim(0,0.055) +
  scale_color_viridis(discrete = T) +
  theme(legend.position = "none")+
  #facet_wrap(~generation) +
  ggtitle("After summing coverage across sexes; autosomes only")+
  geom_vline(xintercept = 100.98)

plot(plot_cov_merged_longer_autosomes)

# convert to a plotly plot
interactive_cov_plot_merged_autosomes <- ggplotly(plot_cov_merged_longer_autosomes)
interactive_cov_plot_merged_autosomes


plot_cov_merged_longer_autosomes2 <- ggplot(cov_merged_longer_autosomes, aes(coverage)) +
  geom_density(aes(color = line_generation)) +
  geom_density(linewidth = 1.5) +
  theme_bw() +
  xlim(0,571) +
    geom_vline(xintercept = 101) +
    geom_vline(xintercept = 0.5*101, linetype = "dashed")+ 
    geom_vline(xintercept = 2*101, linetype = "dashed") +
  theme(legend.position = "top", 
        legend.title = element_blank()) +
  guides(color = guide_legend(ncol = 7)) +
    annotate("text", x = 0, y = .0125, label = "A", size = 5) 



  # Sex chromosome
    # Read the mpileup file with sex chromosome
    # load coverage from the sex chromosome. Note that hash is not used as a comment character!
cov_sex <- read.table("~/sex_ratio_EE/results/5.create_sync_file/sexChr/EE_SR_IndelsRm_RepeatsRm_sexChromosome_subsampled.sync")

# remove N and del (last 2 digits from each cell)
remove_last_digits <- function(x) {
  gsub(".{0,4}$", "", x)
}
cov_sex[4:ncol(cov_sex)] <- lapply(cov_sex[4:ncol(cov_sex)], remove_last_digits)

# add col names
colnames(cov_sex) <- c("LG", "pos", "ref_allele", sample_names)

# sum values within each cell (sum the coverage of all nucleotides)
sum_split_values <- function(x) {
  sapply(strsplit(x, ":"), function(y) sum(as.numeric(y)))
}
  # apply the function
cov_sex <- cov_sex %>%
  mutate(across(contains("males"), sum_split_values))

# subset data and change format to longer
cov_longer_sexChr <- cov_sex %>%
  arrange(LG, pos) %>% # sort the mpileup according to position
  pivot_longer(cols = tidyselect::contains("male"), # change format to longer
    names_to = c("line", "sex", "generation"),
    names_sep = "_",
    values_to = "coverage")%>%
  mutate(line_generation = paste0(line, "_", generation, "_", sex)) %>%
  filter(LG == "LG4") # remove sex chromosome

# calculate mean coverage per line:
cov_longer_sexChr %>%
  group_by(line, generation, sex) %>%
  summarise(mean_coverage = mean(coverage)) 


## Female coverage of the sex chromosome
# Target coverage = 50.41X


# plot female coverage of sex chromosome
cov_plot_sexChr <- cov_longer_sexChr %>%
  filter(sex == "females") %>%
  ggplot(aes(coverage)) +
  geom_density(aes(color = line_generation)) +
  geom_density(linewidth = 1.5) +
  xlim(0,200) +
  ylim(0,0.055) +
  scale_color_viridis(discrete = T) +
  theme(legend.position = "none")+
  ggtitle("Sex chromosome, females") +
  geom_vline(xintercept = 52.4)

plot(cov_plot_sexChr)

# convert to a plotly plot
interactive_cov_plot_sexChr <- ggplotly(cov_plot_sexChr)
interactive_cov_plot_sexChr


Male coverage of the sex chromosome
Tagret coverage = 16.04X


# plot male coverage of sex chromosome
cov_plot_sexChr_males <- cov_longer_sexChr %>%
  filter(sex == "males") %>%
  ggplot(aes(coverage)) +
  geom_density(aes(color = line_generation)) +
  geom_density(linewidth = 1.5) +
  xlim(0,200) +
  ylim(0,0.055) +
  scale_color_viridis(discrete = T) +
  theme(legend.position = "none")+
  ggtitle("Sex chromosome, males") +
  geom_vline(xintercept = 19.6)

plot(cov_plot_sexChr_males)

# convert to a plotly plot
interactive_cov_plot_sexChr_males <- ggplotly(cov_plot_sexChr_males)
interactive_cov_plot_sexChr_males


## check the coverage distribution after male and female coverage on the sex chromosome was summed

# The target coverage for the sex chromosome is 77.5X (vertical line)
# Therefore, the mean coverage across all lines has to be between 38X and 155X.

#  sex chromosome after summing male and female coverage}
# cov_check is a sync file after merging male and female coverage
cov_merged_sex <- read.table("~/sex_ratio_EE/results/5.create_sync_file/sexChr/EE_SR_IndelsRm_RepeatsRm_sexMergedCov_sexChromosome_subsampled.sync")

# remove N and del (last 2 digits from each cell)
cov_merged_sex[4:ncol(cov_merged_sex)] <- lapply(cov_merged_sex[4:ncol(cov_merged_sex)], remove_last_digits)

# add col names
sample_names_check <- unique(c("FB-A_F1", "FB-A_F1", "FB-B_F1", "FB-B_F1", "FB-C_F1", "FB-C_F1"  , "FB-D_F1", "FB-D_F1", "MB-A_F1", "MB-A_F1", "MB-B_F1", "MB-B_F1", "MB-C_F1", "MB-C_F1", "MB-D_F1", "MB-D_F1", "FB-A_F28", "FB-A_F28", "FB-B_F28", "FB-B_F28", "FB-C_F28", "FB-C_F28" , "FB-D_F28", "FB-D_F28", "MB-A_F28", "MB-A_F28", "MB-B_F28", "MB-B_F28", "MB-C_F28", "MB-C_F28", "MB-D_F28", "MB-D_F28"))
colnames(cov_merged_sex) <- c("LG", "pos", "ref_allele", sample_names_check)

# sum values within each cell (sum the coverage of all nucleotides)
cov_merged_sex <- cov_merged_sex %>%
  mutate(across(contains("_F"), sum_split_values))

# change the format to longer
cov_merged_longer_sex <- cov_merged_sex %>%
  arrange(LG, pos) %>% # sort the mpileup according to position
  pivot_longer(cols = tidyselect::contains("_F"), # change format to longer
    names_to = c("line", "generation"),
    names_sep = "_",
    values_to = "coverage")%>%
  mutate(line_generation = paste0(line, "_", generation)) 

# plot coverage
plot_cov_merged_longer_sex <- ggplot(cov_merged_longer_sex, aes(coverage)) +
  geom_density(aes(color = line_generation)) +
  geom_density(linewidth = 1.5) +
  xlim(0,400) +
  #ylim(0,0.055) +
  scale_color_viridis(discrete = T) +
  theme(legend.position = "none")+
  #facet_wrap(~generation) +
  ggtitle("After summing coverage across sexes; sex chromosome only")+
  geom_vline(xintercept = 77.5)

plot(plot_cov_merged_longer_sex)

# convert to a plotly plot
interactive_cov_plot_merged_sex <- ggplotly(plot_cov_merged_longer_sex)
interactive_cov_plot_merged_sex


cov_merged_longer_sex2 <- ggplot(cov_merged_longer_sex, aes(coverage)) +
  geom_density(aes(color = line_generation)) +
  geom_density(linewidth = 1.5) +
  theme_bw() +
    geom_vline(xintercept = 75) +
    geom_vline(xintercept = 0.5*75, linetype = "dashed")+ 
    geom_vline(xintercept = 2*75, linetype = "dashed") +
  theme(legend.position = "none", 
        legend.title = element_blank()) +
    annotate("text", x = -10, y = .0155, label = "B", hjust = 0.5, vjust = 0.5, size = 5) +
    xlim(0,571)

#  combine plots

coverage_combined_plot <- plot_grid( 
  plot_cov_merged_longer_autosomes2, 
  cov_merged_longer_sex2, 
  ncol = 1, 
  rel_heights = c(3, 2))

#print(coverage_combined_plot)

ggsave("~/sex_ratio_EE/results/5.create_sync_file/coverage_distribution_autosomes_sexChr.png", plot = coverage_combined_plot, width = 8, height = 6, dpi = 300, units = "in")


# coverage distribution after filtering
## autosomes

#  autosomes cov after filtering 2}
# cov_check is a sync file after merging male and female coverage
cov_merged_covFilt <- read.table("~/sex_ratio_EE/results/6.filter_by_coverage/withoutAdditionalFiltering/autosomes/all_lines//EE_SR_IndelsRm_RepeatsRm_sexMergedCov_autosomes_filteredCov_subsampled.sync")

# remove N and del (last 2 digits from each cell)
cov_merged_covFilt[4:ncol(cov_merged_covFilt)] <- lapply(cov_merged_covFilt[4:ncol(cov_merged_covFilt)], remove_last_digits)

# add col names
sample_names_check <- unique(c("FB-A_F1", "FB-A_F1", "FB-B_F1", "FB-B_F1", "FB-C_F1", "FB-C_F1"  , "FB-D_F1", "FB-D_F1", "MB-A_F1", "MB-A_F1", "MB-B_F1", "MB-B_F1", "MB-C_F1", "MB-C_F1", "MB-D_F1", "MB-D_F1", "FB-A_F28", "FB-A_F28", "FB-B_F28", "FB-B_F28", "FB-C_F28", "FB-C_F28" , "FB-D_F28", "FB-D_F28", "MB-A_F28", "MB-A_F28", "MB-B_F28", "MB-B_F28", "MB-C_F28", "MB-C_F28", "MB-D_F28", "MB-D_F28"))
colnames(cov_merged_covFilt) <- c("LG", "pos", "ref_allele", sample_names_check)

# sum values within each cell (sum the coverage of all nucleotides)
cov_merged_covFilt <- cov_merged_covFilt %>%
  mutate(across(contains("_F"), sum_split_values))

# change the format to longer
cov_merged_covFilt_longer <- cov_merged_covFilt %>%
  arrange(LG, pos) %>% # sort the mpileup according to position
  pivot_longer(cols = tidyselect::contains("_F"), # change format to longer
    names_to = c("line", "generation"),
    names_sep = "_",
    values_to = "coverage")%>%
  mutate(line_generation = paste0(line, "_", generation)) 

# plot coverage
plot_cov_merged_covFilt_longer <- ggplot(cov_merged_covFilt_longer, aes(coverage)) +
  geom_density(aes(color = line_generation)) +
  geom_density(linewidth = 1.5) +
  xlim(0,400) +
  #ylim(0,0.055) +
  scale_color_viridis(discrete = T) +
  theme(legend.position = "none")+
  #facet_wrap(~generation) +
  ggtitle("After summing coverage across sexes; sex chromosome only", subtitle = "After filtering on coverage")+
  geom_vline(xintercept = 107.24) +
  geom_vline(xintercept = 107.24/2, linetype = "dashed") +
  geom_vline(xintercept = 107.24 * 2, linetype = "dashed")

plot(plot_cov_merged_covFilt_longer)

# convert to a plotly plot
interactive_plot_cov_merged_covFilt_longer <- ggplotly(plot_cov_merged_covFilt_longer)
interactive_plot_cov_merged_covFilt_longer

cov_merged_covFilt_longer %>% 
  group_by(generation, line) %>%
  summarise(mean_coverage_autosomes = mean(coverage))

pos <- position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5, seed = 2)

plot_cov <- cov_merged_covFilt_longer %>% 
  group_by(generation, line) %>%
  summarise(mean_coverage_autosomes = mean(coverage)) %>%
  mutate(treatment = ifelse(str_detect(line, "^F"), "female-biased", "male_biased")) %>%
  ggplot(aes(generation, mean_coverage_autosomes, color = treatment, label = line)) +
    geom_point(position = pos) +
    ggrepel::geom_label_repel(position = pos) +
  theme_minimal() +
    theme(legend.position = "none") 

ggsave("~/sex_ratio_EE/results/5.create_sync_file/mean_coverage_per_line.png", plot = plot_cov, width = 10, height = 7, dpi = 300)



## sex chromosome



#  autosomes cov after filtering}
# cov_check is a sync file after merging male and female coverage

cov_merged_covFilt_sex <- read.table("~/sex_ratio_EE/results/6.filter_by_coverage/withoutAdditionalFiltering/sexChr/all_lines/EE_SR_IndelsRm_RepeatsRm_sexMergedCov_sexChromosome_filteredCov_subsampled.sync")

# remove N and del (last 2 digits from each cell)
cov_merged_covFilt_sex[4:ncol(cov_merged_covFilt_sex)] <- lapply(cov_merged_covFilt_sex[4:ncol(cov_merged_covFilt_sex)], remove_last_digits)

# add col names
sample_names_check <- unique(c("FB-A_F1", "FB-A_F1", "FB-B_F1", "FB-B_F1", "FB-C_F1", "FB-C_F1"  , "FB-D_F1", "FB-D_F1", "MB-A_F1", "MB-A_F1", "MB-B_F1", "MB-B_F1", "MB-C_F1", "MB-C_F1", "MB-D_F1", "MB-D_F1", "FB-A_F28", "FB-A_F28", "FB-B_F28", "FB-B_F28", "FB-C_F28", "FB-C_F28" , "FB-D_F28", "FB-D_F28", "MB-A_F28", "MB-A_F28", "MB-B_F28", "MB-B_F28", "MB-C_F28", "MB-C_F28", "MB-D_F28", "MB-D_F28"))
colnames(cov_merged_covFilt_sex) <- c("LG", "pos", "ref_allele", sample_names_check)

# sum values within each cell (sum the coverage of all nucleotides)
cov_merged_covFilt_sex <- cov_merged_covFilt_sex %>%
  mutate(across(contains("_F"), sum_split_values))

# change the format to longer
cov_merged_covFilt_sex_longer <- cov_merged_covFilt_sex %>%
  arrange(LG, pos) %>% # sort the mpileup according to position
  pivot_longer(cols = tidyselect::contains("_F"), # change format to longer
    names_to = c("line", "generation"),
    names_sep = "_",
    values_to = "coverage")%>%
  mutate(line_generation = paste0(line, "_", generation)) 

# plot coverage
plot_cov_merged_covFilt_sex_longer <- ggplot(cov_merged_covFilt_sex_longer, aes(coverage)) +
  geom_density(aes(color = line_generation)) +
  geom_density(linewidth = 1.5) +
  xlim(0,400) +
  #ylim(0,0.055) +
  scale_color_viridis(discrete = T) +
  theme(legend.position = "none")+
  #facet_wrap(~generation) +
  ggtitle("After summing coverage across sexes; sex chromosome only", subtitle = "After filtering on coverage")+
  geom_vline(xintercept = 73.58)  +
  geom_vline(xintercept = 73.58/2, linetype = "dashed") +
  geom_vline(xintercept = 73.58 * 2, linetype = "dashed")

plot(plot_cov_merged_covFilt_sex_longer)

# convert to a plotly plot
interactive_plot_cov_merged_covFilt_sex_longer <- ggplotly(plot_cov_merged_covFilt_sex_longer)
interactive_plot_cov_merged_covFilt_sex_longer

cov_merged_covFilt_sex_longer %>%
  group_by(generation, line) %>%
  summarise(mean_coverage = mean(coverage)) %>%
  mutate(treatment = ifelse(str_detect(line, "^F"), "female-biased", "male-biased")) %>%
ggplot(aes(generation, mean_coverage, color = treatment)) +
  geom_jitter(width = .2, seed = 123) +
  ggrepel::geom_label_repel(aes(label =line), seed = 123)

cov_merged_covFilt_sex_longer %>%
  group_by(generation, line) %>%
  summarise(mean_coverage = mean(coverage)) %>%
  mutate(treatment = ifelse(str_detect(line, "^F"), "female-biased", "male-biased")) %>%
ggplot() +
geom_point(aes(generation, mean_coverage, color = treatment),
            alpha = .5, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5), size = 2) +
geom_text_repel(aes(generation, mean_coverage, color = treatment, label = line),
                  position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5), size = 3, fill = NA) +
  theme_bw()



