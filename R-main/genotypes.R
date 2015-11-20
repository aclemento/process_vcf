#### Filtering the vcf and outputting genotypes ####

# Filter the vcf file
system('/usr/local/bin/vcflib/vcffilter -f "TYPE = snp & MQM > 30 & MQMR > 30 & NS > 15 & AC > 1 & DP > 499 & DP < 10000" inst/extdata/smar_freebayes.vcf > output/smar_filtered.vcf')

# Export plink ped and map
system('/usr/local/bin/vcftools/vcftools --vcf output/smar_filtered.vcf --out output/smar_plink --plink')



## GENEPOP ##
# Dump the genotypes out in 0, 1, 2 format
system('/usr/local/bin/vcftools/vcftools --vcf output/smar_filtered.vcf --out output/smar --012')

# get the indivs
gsIndivs <- scan("output/smar.012.indv", what = "character")
N <- length(gsIndivs)

# get the markers
tmp <- scan("output/smar.012.pos", what = "character")
gsMarkers <- paste("X", tmp[c(T,F)], tmp[c(F,T)], sep = "_")
M <- length(gsMarkers)

# now make a hashie thing to change some allele names and make it all two-columny
alles <- c("0101", "0102", "0202", "0000")
names(alles) <- c("0", "1", "2", "-1")

# then read them in and make a matrix of two-column genotypes
genos <- alles[scan("output/smar.012", what = "character")] %>% 
  matrix(., nrow = N, byrow = TRUE)  

rownames(genos) <- gsIndivs
genos <- genos[,-1]  # note that we tear off the column that gives the index of each individual

# now we need to sort the rows so that populations are together in case they arent already grouped
# togehter.  Just a simple order on the rownames.
genos <- genos[order(rownames(genos)),]

# once we have done that we should modify the row names to contain the genepop comma 
# Note that this is specific to the naming format with 4 characters (i.e. huda01)
pops <- str_sub(rownames(genos), 1, 4)
firsties <- c(1, 1 + which(diff(as.numeric(factor(pops))) != 0))
rownames(genos)[firsties] <- paste("POP ", "\n", rownames(genos)[firsties], sep = "", " , ")

# then we just spooge out the preamble and the genotypes.
cat("Title line:\"smar_all_GP2.txt\"", "\n", sep = "", file = "output/smar_all_GP2.txt")
cat(gsMarkers, sep = "\n", file = "output/smar_all_GP2.txt", append = TRUE)
write.table(genos, col.names = FALSE, quote = FALSE, sep = " ", file = "output/smar_all_GP2.txt", append = TRUE)



## Two-column format ##
# now make a hashie thing to change some allele names and make it all two-columny
alles <- c("1 1", "1 2", "2 2", "0 0")
names(alles) <- c("0", "1", "2", "-1")

# then read them in and make a matrix of two-column genotypes
genos <- alles[scan("output/smar.012", what = "character")] %>% 
  matrix(., nrow = N, byrow = TRUE)  

rownames(genos) <- gsIndivs
genos <- genos[,-1]  # note that we tear off the column that gives the index of each individual

# We also need a header for our two-column format with the locus name duplicated (i.e. loc loc.1)
# We will write a function for interleaving
interleave <- function(v1,v2)
{
  ord1 <- 2*(1:length(v1))-1
  ord2 <- 2*(1:length(v2))
  c(v1,v2)[order(c(ord1,ord2))]
}
# And create another version of gsMarkers

interleave(rep(1,5),rep(3,8))



# now we need to sort the rows so that populations are together in case they arent already grouped
# togehter.  Just a simple order on the rownames.
genos <- genos[order(rownames(genos)),]

# once we have done that we should modify the row names to contain POP and comma delimiting individuals 
# Note that this is specific to the naming format with 4 characters (i.e. huda01)
pops <- str_sub(rownames(genos), 1, 4)
firsties <- c(1, 1 + which(diff(as.numeric(factor(pops))) != 0))
rownames(genos)[firsties] <- paste("POP ", "\n", rownames(genos)[firsties], sep = "", " , ")

# then we just spooge out the preamble and the genotypes.
cat("Title line:\"smar_all_GP2.txt\"", "\n", sep = "", file = "output/smar_all_GP2.txt")
cat(gsMarkers, sep = "\n", file = "output/smar_all_GP2.txt", append = TRUE)
write.table(genos, col.names = FALSE, quote = FALSE, sep = " ", file = "output/smar_all_GP2.txt", append = TRUE)



## We would also like a dataset with only 1 SNP per locus
