## Run GLM on allele frequency counts to identify diverged SNPs
# author: Jonathan Parrett jonathan.parrett(_a_)amu.edu.pl, with some changes/corrections from Sebastian Chmielewski
  # due to the size of the file, the input should be divided into smaller chunks prior to GLM and remove the header
args <- commandArgs(trailingOnly=TRUE)
infile <- args[1] # from the command line
data <- read.delim(infile, sep="", header = F)

# test
  # data <- read.delim("~/sex_ratio_EE/results/7.glm_afc/polymorphicSnps_EE_SR_IndelsRm_RepeatsRm_sexMergedCov_allChr_filteredCov_head.txt", sep = "", h = F) # test

colnames(data) <- c("Contig", "Pos", "SNPID", "line", "MajAllele", "A", "T", "C", "G", "coverage", "selection", "MajAlleleCount")
data$selection <- as.factor(data$selection)
# add SNP id
SNPIDs <- unique(data$SNPID) 

# tests
data$check1<-ifelse(data$coverage-data$MajAlleleCount==0,0,1) # check if the Maj allele count equals the line coverage
data$check2<-ifelse(data$MajAlleleCount==0,0,1) # check if Maj allele count is 0

# initialise vectors
pvalues<-c()
SNPcontig <- c()
SNPpos<-c()
SNPID.ID<-c()
SNPexactPos<-c()

for (i in SNPIDs){
  temp_data <- data[data$SNPID==i,] # make df with only 1 SNP
  check1sum <- sum(temp_data$check1) 
  check2sum <- sum(temp_data$check2)
  if(check1sum==8 & check2sum==8){model <- glm(cbind(MajAlleleCount,(coverage-MajAlleleCount)) ~ selection, family=quasibinomial, data=temp_data)
  } # standard GLM
  if(check1sum<8 | check2sum<8){ # add one to the coverage and MajAlleleCount if not all lines fulfill the statements
    temp_data$MajAlleleCount <- temp_data$MajAlleleCount + 1
    temp_data$coverage <- temp_data$coverage + 2
    model <- glm(cbind(MajAlleleCount,(coverage-MajAlleleCount)) ~ selection, family=quasibinomial, data=temp_data)
  } 
  pvalues <- c(pvalues, coef(summary(model))[2, "Pr(>|t|)"])
  SNPcontig <- c(SNPcontig, as.character(temp_data$Contig[1]))
  SNPpos <- c(SNPpos, temp_data$Pos[1])
  SNPID.ID <- c(SNPID.ID, temp_data$SNPID[1])
  SNPexactPos <- c(SNPexactPos, as.character(temp_data$Exact.Pos[1]))
}
options(max.print=1000000, width = 10000)

pvaluesandSNPID2<-data.frame(SNPID.ID,SNPcontig,SNPpos,pvalues)

if (nrow(pvaluesandSNPID2) > 0) {
  print(pvaluesandSNPID2, row.names=FALSE)
}
