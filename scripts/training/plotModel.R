#!/usr/bin/env Rscript
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

spline.infile <- args[1]
data.basedir <- args[2]
precursor <- args[3]
max.sulfur <- args[4]
outfile <- args[5]
max.mass <- as.numeric(args[6])

data <- data.frame()

if (max.sulfur == -1) {
  1
} else {
  for (sulfur in 0:max.sulfur) {
    scatter.infile <- paste(data.basedir, "S", toString(sulfur), "/data/Precursor", precursor, ".tab", sep="")
    data.tmp <- read.table(scatter.infile, header=T, sep="\t")
    data.tmp <- subset(data.tmp, data.tmp$precursor.mass <= max.mass)
    data.tmp$S <- sulfur
    data <- rbind(data, data.tmp)
  }
}
  
  
setEPS()
postscript(outfile, width=9, height=6)

p <- ggplot(data, aes(x=precursor.mass, y=probability, color=as.factor(S)))

print(
  p
  + geom_point(alpha=0.1)
  )

dev.off()