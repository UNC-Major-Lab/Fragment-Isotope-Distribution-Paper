#!/usr/bin/env Rscript
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

infile = args[1]
outfile = args[2]

data <- read.table(infile,header=F,sep="\t")

p <- ggplot(data=data, aes(x=V2, y=V3, group=V1, fill=factor(V1)))

pdf(outfile, width=9, height=6)

p + geom_bar(stat='identity', width=.01) + facet_wrap(~ V1, scale='free_y') + xlab("residual") + ylab("frequency") + ggtitle("Histogram of residuals (averagine - expected)")

dev.off()
