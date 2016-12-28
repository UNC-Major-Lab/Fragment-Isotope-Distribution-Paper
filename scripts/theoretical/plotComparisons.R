#!/usr/bin/env Rscript
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

infile = args[1]
outfile = args[2]

data <- read.table(infile,header=F,sep="\t")

p <- ggplot(data=data, aes(x=V2, y=V4, group=V3, fill=factor(V3)))

pdf(filename=outfile, width=600, height=400, units="px")

p 
+ geom_bar(stat='identity', width=.01) 
+ facet_wrap(~ V3, scale='free_y')
+ xlab("residual")
+ ylab("frequency")
+ ggtitle("Histogram of residuals (averagine fragment isotope probability - exact fragment isotope probability)")


dev.off()
