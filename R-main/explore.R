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
x <- readLines("inst/extdata/satro6_noNs_noMNP_noComplex.vcf", n = 1000)
header_line <- min(which(str_detect(x, "#CHROM")))
vcf <- read_tsv("inst/extdata/satro6_noNs_noMNP_noComplex.vcf", skip = header_line - 1) %>%
#vcf <- read_tsv("output/smar_filtered_ns32_ac8.vcf", skip = header_line - 1) %>%
  tbl_df
names(vcf)[1] <- "CHROM"


#### process the INFO fields ####


# pick CHROM POS and INFO
cpi <- vcf %>%
 select(CHROM, POS, INFO)

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

ggsave(g, filename = "histo_matrix.pdf", width = 20, height = 18)

# plot a single field
ggplot(filter(clean, info_field == "DP"), aes(x = info_numeric)) +
  geom_histogram(binwidth = 100, colour = "blue") + xlim(0,25000)


#### Make a plotmatrix of interesting info values ####
interesting <- c("DP", "AC", "AF", "AN", "NS", "MQM", "MQMR")


inty_wide <- clean %>%
  filter(info_field %in% interesting) %>%
  unite(CHROMPOS, CHROM, POS, remove = FALSE) %>%
  select(CHROMPOS, CHROM, POS, info_field, info_numeric) %>%
  spread(info_field, info_numeric) %>%
  mutate(status = ifelse(MQM > 30 & MQMR > 30 & NS >= 16 & AC >= 2 & DP >= 500 & DP < 10000, "keep", "toss"))


pdf(file = "vcf_matrix.pdf", width = 24, height = 20)
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
ggplot(inty_wide, aes(x = AF)) +
  geom_density( fill = "blue", alpha = 0.3) +
  geom_density(data = inty_wide %>% filter(status == "keep"), alpha = 0.3, fill = "red") 
 

