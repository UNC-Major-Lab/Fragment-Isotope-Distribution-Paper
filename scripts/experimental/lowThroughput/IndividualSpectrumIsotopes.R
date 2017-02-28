#!/usr/bin/env Rscript
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

infile1 <- args[1]
infile2 <- args[2]
infile3 <- args[3]
outfile <- args[4]

data <- read.table(infile1, sep="\t", header=T)
data.calc <- read.table(infile2, sep="\t", header=T)
data.scores <- read.table(infile3, sep="\t", header=T)

data <- subset(data, (data$ion.index == 67 | data$ion.index == 55) & data$isotope.range <= 5)
data.calc <- subset(data.calc, (data.calc$ion.index == 67 | data.calc$ion.index == 55) & data.calc$isotope.range <= 5)
data.scores <- subset(data.scores, (data.scores$ion.index == 67 | data.scores$ion.index == 55) & data.scores$isotope.range <= 5)

p <- ggplot(data=data, aes(x=mz,y=int))

pdf(outfile, width=8.5, height=11)

print(
  p
  + geom_line()
  + geom_point(data=data.calc, aes(group=method, color=method), size=3, shape=1)
  + facet_grid(isotope.range ~ ion.name, scale="free")
  + ylab("Intensity normalized to base peak")
  + xlab("m/z")
  + geom_text(data=data.scores, aes(group=method, color=method, x=x, y=y, label=label) size=1)
)

dev.off()