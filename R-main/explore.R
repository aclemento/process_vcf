library(stringr)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(GGally)

setwd("~/Genetics_Lab_Data/rproj/process_vcf")


#### read in vcf file ####


# Download vcf from Google Drive and put it in inst/extdata/

# Decompress vcf file for reading into R
system('/bin/zcat inst/extdata/smar_freebayes.vcf.gz > inst/extdata/smar_freebayes.vcf')

# find header line of vcf file then read it in
#x <- readLines("inst/extdata/satro6_filter_QUAL20_AF02_TYPEsnp.vcf", n = 1000)
x <- readLines("inst/extdata/satro_144_noMNP_noComplex.vcf", n = 1000)
header_line <- min(which(str_detect(x, "#CHROM")))
vcf <- read_tsv("inst/extdata/satro_144_noMNP_noComplex.vcf", skip = header_line - 1) %>%
#vcf <- read_tsv("output/smar_filtered_ns32_ac8.vcf", skip = header_line - 1) %>%
  tbl_df
names(vcf)[1] <- "CHROM"


#### process the INFO fields ####


# pick CHROM POS QUAL and INFO
cpi <- vcf %>%
 select(CHROM, POS, QUAL, INFO)

# parse the INFO string in the first row and hope it is the same for all of them.  Return
# a string of names of the INFO fields
info_names <- str_split(cpi$INFO[1], "[;=]")[[1]] %>% 
  matrix(ncol = 2, byrow = TRUE) %>%
  "["(,1)

tmp <- cpi %>%
  separate(INFO, sep = ";", into = info_names) %>% 
  gather(key = "info_field", value = "info_value", -CHROM, -POS) %>%
  mutate(info_value = str_replace(info_value, "^.*=", "")) 
# This is the fully melted INFO fields and a good place to summarize TYPE
filter(tmp, info_field=="TYPE") %>% group_by(info_value) %>% summarise(count = n()) %>% View()
filter(tmp, info_field=="TYPE") %>% group_by(info_value) %>% View()

qstrs <- tmp %>%
  filter(info_field %in% c("TYPE", "CIGAR")) %>%
  unite(CHROMPOS, CHROM, POS, sep =":") %>%
  spread(info_field, info_value) %>%
  separate(CHROMPOS, c("CHROM", "POS"), sep = ":", convert = TRUE)

nums <- tmp %>%
  filter(!(info_field %in% c("TYPE", "CIGAR", "technology.ILLUMINA"))) %>% 
  mutate(info_numeric = as.numeric(info_value))

togeth <- left_join(nums, qstrs)

clean <- togeth %>%
  filter(TYPE == "snp")

# plot all fields
g <- ggplot(clean, aes(x = info_numeric)) +
  geom_histogram(colour = "blue") +
  facet_wrap(~ info_field, ncol = 7, scales = "free") 
g

ggsave(g, filename = "histo_matrix.pdf", width = 20, height = 18)

# plot a single field
ggplot(filter(clean, info_field == "DP"), aes(x = info_numeric)) +
  geom_histogram(binwidth = 100, colour = "blue") + xlim(0,25000)


#### Make a plotmatrix of interesting info values ####


interesting <- c("QUAL", "DP", "AC", "AF", "AB", "NS", "MQM", "MQMR")

inty_wide <- clean %>%
  filter(info_field %in% interesting) %>%
  unite(CHROMPOS, CHROM, POS, sep =":", remove = FALSE) %>%
  select(CHROMPOS, CHROM, POS, info_field, info_numeric) %>%
  spread(info_field, info_numeric) %>%
  mutate(status = ifelse(QUAL > 20 & AF > 0.02 & AF < 0.98 & AB > 0.1 & AB < 0.9 & MQM > 20, "keep", "toss"))

pdf(file = "output/satro_144_vcf_matrix.pdf", width = 24, height = 20)
ggpairs(inty_wide, 
        2:8,
        alpha = 0.1,
        color = "status")
dev.off()

# quickly summarize number of snps and number of loci
inty_wide %>% filter(status == "keep") %>%
  nrow()

inty_wide %>% filter(status == "keep") %>%
  select(CHROM) %>%
  unique() %>%
  nrow()

# Can also look at the before- and after-filtering of a single variable 
require(gridExtra)
a <- ggplot(inty_wide, aes(x = AF)) +
  geom_density( fill = "blue", alpha = 0.3) +
  geom_density(data = inty_wide %>% filter(status == "keep"), alpha = 0.3, fill = "red") 
b <- ggplot(inty_wide, aes(x = AB)) +
  geom_density( fill = "blue", alpha = 0.3) +
  geom_density(data = inty_wide %>% filter(status == "keep"), alpha = 0.3, fill = "red") 
c <- ggplot(inty_wide, aes(x = DP)) +
  geom_density( fill = "blue", alpha = 0.3) +
  geom_density(data = inty_wide %>% filter(status == "keep"), alpha = 0.3, fill = "red") 
d <- ggplot(inty_wide, aes(x = QUAL)) +
  geom_density( fill = "blue", alpha = 0.3) +
  geom_density(data = inty_wide %>% filter(status == "keep"), alpha = 0.3, fill = "red") 
grid.arrange(a, b, c, d, ncol=2)


####Calculate HWE####


##Call vcftools##
system("/usr/local/bin/vcftools --vcf inst/extdata/satro_144_filtered.vcf --hardy --max-missing 0.5 --out intermediates/satro_144_filtered2")

## A function to read vcf hardy output and put it into usable long format ##
long_hardy <- function(file) {
  x <- read.table(file, sep ="\t", header = TRUE, stringsAsFactors = FALSE) %>%
    tbl_df %>%
    separate(OBS.HOM1.HET.HOM2., c("obs_Hom1", "obs_Het", "obs_Hom2"), sep ="/") %>%
    separate(E.HOM1.HET.HOM2., c("exp_Hom1", "exp_Het", "exp_Hom2"), sep ="/")
  
  # now put it into long format
  obs <- x %>% select(CHR, POS, starts_with("obs_")) %>%
    gather(., obs_var, obs_cnt, starts_with("obs_")) %>%
    mutate(geno = str_replace(obs_var, "^obs_", "")) %>%
    mutate(obs_cnt = as.numeric(obs_cnt)) %>%
    select(-obs_var)
  
  # now add the total number of observed individuals in there
  ntot <- obs %>% group_by(CHR, POS) %>% summarise(nindiv = sum(obs_cnt))
  
  exp <- x %>% select(CHR, POS, starts_with("exp_")) %>%
    gather(., exp_var, exp_cnt, starts_with("exp_")) %>%
    mutate(geno = str_replace(exp_var, "^exp_", "")) %>%
    mutate(exp_cnt = as.numeric(exp_cnt)) %>%
    select(-exp_var)
  
  
  cnts <- inner_join(obs, exp) %>%
    mutate(geno = factor(geno, levels = c("Hom1", "Het", "Hom2"))) %>%
    inner_join(., ntot)
  
  list(cnts = cnts, p_etc = inner_join(x, ntot))
  
}

hwe <- long_hardy("intermediates/satro_144_filtered2.hwe")

# Let's see how close we are to equilibrium overall
hwplot <- ggplot(hwe$cnts, aes(x = exp_cnt, y = obs_cnt, colour = geno)) +
  geom_jitter(alpha = 0.75, size=5) +
  geom_abline(intercept = 0, slope = 1)

hwplot  

##Call plink##
# It appears that vcftools has a bug that prevents HWE being calculated for loci with any missing data
# vcftools also strips out triallelic snps - fortunately there are only three

system("/usr/local/bin/vcftools --vcf inst/extdata/")
system("/usr/local/bin/plink --noweb --file intermedaiates/satro_144 --hardy")

hwplink <- read.table(file="intermediates/plink.hwe", stringsAsFactors=F, header=T, sep="") %>% 
  tbl_df() %>% 
  select(-CHR) %>%
  filter(., TEST=="ALL") %>%
  mutate(status = ifelse(P < 0.0001, "toss", "keep"))

hwplplot <- ggplot(hwplink, aes(x = E.HET., y = O.HET., colour = status)) +
  geom_jitter(alpha = 0.75, size=5) +
  geom_abline(intercept = 0, slope = 1)

hwplplot  

