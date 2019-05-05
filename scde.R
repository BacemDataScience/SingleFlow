#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(scde)
require(reticulate)
source_python('pickle_reader.py')

file <- args[1]
df_name <- unlist(strsplit(file, "[.]"))[-2]
#print(df_name)

df <- read_pickle_file(file)

#load(file)
#df <- get(df_name)

print('data frame loaded')

df <- apply(df,2,function(x) {storage.mode(x) <- 'integer'; x})

print('data frame to integer done')

#write.table(df, file = paste0(df_name, '.csv'))

clusters <- substr(colnames(df), 1, 1)
levels <- c(clusters[1], clusters[length(clusters)])

sg <- factor(clusters, levels = levels)
#### the group factor should be named accordingly
names(sg) <- colnames(df)

o.ifm <- scde.error.models(counts = df, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
o.prior <- scde.expression.prior(models = o.ifm, counts = df, length.out = 400, show.plot = FALSE)
groups <- factor(clusters, levels = levels)
names(groups) <- row.names(o.ifm)

#### run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, df, o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)
#### ediff[order(ediff\$Z, decreasing  =  TRUE), ]
#### ediff[abs(ediff\$Z) > 1.96, ]
write.table(ediff, file=paste0(df_name, ".csv"), append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
