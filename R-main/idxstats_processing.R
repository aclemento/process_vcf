# Summarizing idxstats
install.packages("dplyr")
install.packages("ggplot2")
install.packages("RColorBrewer")
library(stringr)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)

# We have 96 files with idxstats that we need to read in, add an individual identifier column and then plot (raster?)

setwd("~/Genetics_Lab_Data/rproj/process_vcf/inst/idxdata/")

# read in the file list that I created with "ls -1 s* > files.txt"
files <- read.table("files.txt", stringsAsFactors=F)

# read in the data files each to a list component, grab the first part of the file name and name the list components
#stats <- lapply(files[,1], function(x) read.table(x, nrows=95))
stats <- lapply(files[,1], function(x) read.table(x, nrows=95, stringsAsFactors=F))
temp <- matrix(unlist(strsplit(files[,1], "_")), ncol=2, byrow=T)
names(stats) <- temp[,1]

# need to add a column to each list component that is the list name
idx <- list()
for (i in temp[,1]) {y <- rep(names(stats[i]),95)
                             idx[[i]]<- cbind(stats[[i]], y)}

# now merge into a single tidy data frame
idx_df <- tbl_df(bind_rows(idx))
names(idx_df) <- c("loc", "len", "reads", "unmapd", "ind")

# can we refactor/reorder the locs and/or the indivs? 
# first the locs
locs <- idx_df %>% group_by(loc) %>% summarise(sum(reads))
names(locs) <- c("loc", "total")
locs <- locs[order(locs$total),]
idx_df <- transform(idx_df, loc = factor(loc, levels = locs$loc))
# now individuals
inds <- idx_df %>% group_by(ind) %>% summarise(sum(reads))
names(inds) <- c("ind", "total")
inds <- inds[order(inds$total),]
idx_df <- transform(idx_df, ind = factor(ind, levels = inds$ind))

# and plot it up 
ggplot(idx_df, aes(loc, ind, fill = reads)) + geom_raster() + scale_fill_gradient(low="red", high="green", trans='log') 

# export the actual table
matrix_df <- select(idx_df, -len, -unmapd) %>% spread(., ind, reads)
write.table(matrix_df, file="gtseq6_idx_matrix.txt", quote=F, sep="\t")
