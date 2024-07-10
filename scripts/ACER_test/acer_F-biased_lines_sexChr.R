### run adapted CMH-test using ACER package
# author: Sebastian Chmielewski, Evolutionary Biology Group, AMU; sebchm_a_amu.edu.pl
# female-biased lines (diff in the Ne compared to male-biased lines)

suppressWarnings(library(ACER))
suppressWarnings(library(poolSeq))

args = commandArgs(trailingOnly=TRUE)
infile=args[1]

mySync <- read.sync(infile, gen=c(rep(1,4), rep(28,4)), repl=rep(c(1,2,3,4),2),polarization = "minor")

# allele frequency matrix
af <- af(mySync, repl = rep(c(1,2,3,4),2), gen = c(rep(1,4), rep(28,4)))

# coverage matrix
cov <- coverage(mySync, repl = rep(c(1,2,3,4),2), gen = c(rep(1,4), rep(28,4)))

# load effective population size
Ne <- read.table("~/sex_ratio_EE/results/10.estimate_Ne/Ne_est_sexChr_all_lines.txt", h = T)

pval_ACER <- adapted.cmh.test(freq = af, 
                           coverage = cov, 
                           gen = c(1,28),
                           repl = c(1,2,3,4), 
                           Ne=c(Ne[1,2], Ne[2,2], Ne[3,2], Ne[4,2]),
                           poolSize = rep(150, 8),
                           IntGen = FALSE, 
                           order = 0, 
                           mincov = 38, 
                           RetVal = 0, 
                           correct = TRUE)

SNP_pvalue <- data.frame(SNP = row.names(af), 
                         pvalue = pval_ACER)

options(max.print=1000000, width = 10000)


if (nrow(SNP_pvalue) > 0) {
  print(SNP_pvalue, row.names=FALSE)
}
