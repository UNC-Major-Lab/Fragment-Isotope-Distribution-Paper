#!/usr/bin/env Rscript
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

spline.infile <- args[1]
data.basedir <- args[2]
precursor <- args[3]
max.sulfur <- args[4]
outfile <- args[5]

#data <- read.table(infile,header=F,sep="\t")

scatter.infile <- paste(data.basedir, "/S", "0", "/data/Precursor", precursor, ".tab", sep="")
data <- read.table(scatter.infile, header=T, sep="\t")

setEPS()
postscript(outfile, width=9, height=6)

p <- ggplot(data, aes(x=precursor.mass, y=probability))

print(
  p
  + geom_point()
  )

dev.off()