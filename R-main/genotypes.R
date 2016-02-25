library(stringr)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)

setwd("~/Genetics_Lab_Data/rproj/process_vcf")

#### Filtering the vcf and outputting genotypes ####

# Filter the vcf file
system('/usr/local/bin/vcflib/vcffilter -f "TYPE = snp & DP > 300 & AB > 0.1 & AB < 0.9" inst/extdata/sebastes_snps_only.recode.vcf > intermediates/sebastes_snps_filtered.vcf')
#system('/usr/local/bin/vcflib/vcffilter -f "NS > 32" smarmo.vcf > smarmo_py_filtered_32.vcf')

# Export plink ped and map
#system('/usr/local/bin/vcftools/vcftools --vcf output/smar_filtered_16.vcf --out output/smar_all_plink_16 --plink')


# Dump the genotypes out in 0, 1, 2 format
#system('/usr/local/bin/vcftools/vcftools --vcf output/smar_filtered_60.vcf --out output/smar --012')
#system('/usr/local/bin/vcftools/vcftools --vcf output/smar_filtered_ns32_ac8.vcf --out output/smar --positions output/white_1snp_ns32_ac8.txt --012')
system('/usr/local/bin/vcftools --vcf intermediates/sebastes_spp_filtered.vcf --out intermediates/sebastes --012')

# get the indivs
gsIndivs <- scan("intermediates/sebastes.012.indv", what = "character")
N <- length(gsIndivs)

# get the markers
tmp <- scan("intermediates/sebastes.012.pos", what = "character")
gsMarkers <- paste("X", tmp[c(T,F)], tmp[c(F,T)], sep = "_")
M <- length(gsMarkers)



## GENEPOP ##
# now make a hashie thing to change some allele names and make it all two-columny
alles <- c("0101", "0102", "0202", "0000")
names(alles) <- c("0", "1", "2", "-1")

# then read them in and make a matrix of two-column genotypes
genos <- alles[scan("output/smar.012", what = "character")] %>% 
  matrix(., nrow = N, byrow = TRUE)  

rownames(genos) <- paste(gsIndivs, ",", sep=" ")
genos <- genos[,-1]  # note that we tear off the column that gives the index of each individual

# now we need to sort the rows so that populations are together in case they arent already grouped
# togehter.  Just a simple order on the rownames.
genos <- genos[order(rownames(genos)),]

# once we have done that we should modify the row names to contain the genepop comma 
# Note that this is specific to the naming format with 4 characters (i.e. huda01)
pops <- str_sub(rownames(genos), 1, 4)
firsties <- c(1, 1 + which(diff(as.numeric(factor(pops))) != 0))
rownames(genos)[firsties] <- paste("POP", "\n", rownames(genos)[firsties], sep = "")

# then we just spooge out the preamble and the genotypes.
cat("Title line:\"smar_pymap_ns32_ac8_1snp.genepop\"", "\n", sep = "", file = "output/smar_pymap_ns32_ac8_1snp.genepop")
cat(gsMarkers, sep = "\n", file = "output/smar_pymap_ns32_ac8_1snp.genepop", append = TRUE)
write.table(genos, col.names = FALSE, quote = FALSE, sep = " ", file = "output/smar_pymap_ns32_ac8_1snp.genepop", append = TRUE)



## Two-column format ##
# now make a hashie thing to change some allele names and make it all two-columny
alles <- c("1\t1", "1\t2", "2\t2", "0\t0")
names(alles) <- c("0", "1", "2", "-1")

# then read them in and make a matrix of two-column genotypes
genos <- alles[scan("intermediates/sebastes.012", what = "character")] %>% 
  matrix(., nrow = N, byrow = TRUE)  

rownames(genos) <- gsIndivs
genos <- genos[,-1]  # note that we tear off the column that gives the index of each individual

# We also need a header for our two-column format with the locus name duplicated (i.e. loc1 loc1)

# We will pilfer a function for interleaving
interleave <- function(v1,v2)
{
  ord1 <- 2*(1:length(v1))-1
  ord2 <- 2*(1:length(v2))
  c(v1,v2)[order(c(ord1,ord2))]
}

twocolMarkers <- interleave(gsMarkers,gsMarkers) 


# now we need to sort the rows so that populations are together in case they arent already grouped
# togehter.  Just a simple order on the rownames.
genos <- genos[order(rownames(genos)),]

# then we just spooge out the locus list and the genotypes.
cat("Inds", twocolMarkers, "\n", sep = "\t", file = "output/sebastes_refiltered_TK.txt")
write.table(genos, col.names = FALSE, quote = FALSE, sep = "\t", file = "output/sebastes_refiltered_TK.txt", append = TRUE)



## We may also like a dataset with only 1 SNP per locus
cp <- vcf %>% select(CHROM, POS)
one <- cp %>% group_by(CHROM) %>% sample_n(., 1)
write.table(one, file="output/onesnp_ns32_ac8.txt", quote=F, row.names=F, col.names=F, sep="\t")
