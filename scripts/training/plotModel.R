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
data.spline <- data.frame()

if (max.sulfur == -1) {
  scatter.infile <- paste(data.basedir, "Average_Spline", "/data/Precursor", precursor, ".tab", sep="")
  data <- read.table(scatter.infile, header=T, sep="\t")
  data <- subset(data, data$precursor.mass <= max.mass)
  data$S <- -1
  
  spline.infile <- paste(data.basedir, "Average_Spline", "/spline/eval/Precursor", precursor, ".tab", sep="")
  data.spline <- read.table(spline.infile, header=T, sep="\t")
  data.spline <- subset(data.spline, data.spline$precursor.mass <= max.mass)
  data.spline$S <- -1
} else {
  for (sulfur in 0:max.sulfur) {
    scatter.infile <- paste(data.basedir, "S", toString(sulfur), "/data/Precursor", precursor, ".tab", sep="")
    data.tmp <- read.table(scatter.infile, header=T, sep="\t")
    data.tmp <- subset(data.tmp, data.tmp$precursor.mass <= max.mass)
    data.tmp$S <- sulfur
    data <- rbind(data, data.tmp)
    
    spline.infile <- paste(data.basedir, "S", toString(sulfur), "/spline/eval/Precursor", precursor, ".tab", sep="")
    data.spline.tmp <- read.table(spline.infile, header=T, sep="\t")
    data.spline.tmp <- subset(data.spline/tmp, data.spline.tmp$precursor.mass <= max.mass)
    data.spline.tmp$S <- sulfur
    data.spline <- rbind(data.spline, data.spline.tmp)
  }
}
  
  
setEPS()
postscript(outfile, width=9, height=6)

p <- ggplot(data, aes(x=precursor.mass, y=probability, color=as.factor(S)))

print(
  p
  + geom_point(shape=1)
  + geom_line(data=data.spline, aes(x=precursor.mass, y=probability))
  )

dev.off()