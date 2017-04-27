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
data.spline.all <- data.frame()

spline.infile <- paste(data.basedir, "Average_Spline", "/spline/eval/Precursor", precursor, ".tab", sep="")
data.spline.all <- read.table(spline.infile, header=T, sep="\t")
data.spline.all <- subset(data.spline, data.spline$precursor.mass <= max.mass)
data.spline.all$S <- "All"

if (max.sulfur == -1) {
  scatter.infile <- paste(data.basedir, "Average_Spline", "/data/Precursor", precursor, ".tab", sep="")
  data <- read.table(scatter.infile, header=T, sep="\t")
  data <- subset(data, data$precursor.mass <= max.mass)
  data$S <- "All"
  data.spline.all$line.type <- "solid"
  data.spline.all$color <- "black"
} else {
  for (sulfur in 0:max.sulfur) {
    scatter.infile <- paste(data.basedir, "S", toString(sulfur), "/data/Precursor", precursor, ".tab", sep="")
    data.tmp <- read.table(scatter.infile, header=T, sep="\t")
    data.tmp <- subset(data.tmp, data.tmp$precursor.mass <= max.mass)
    data.tmp$S <- toString(sulfur)
    data <- rbind(data, data.tmp)
    
    spline.infile <- paste(data.basedir, "S", toString(sulfur), "/spline/eval/Precursor", precursor, ".tab", sep="")
    data.spline.tmp <- read.table(spline.infile, header=T, sep="\t")
    data.spline.tmp <- subset(data.spline.tmp, data.spline.tmp$precursor.mass <= max.mass)
    data.spline.tmp$S <- toString(sulfur)
    data.spline.tmp$line.type <- 'solid'
    data.spline <- rbind(data.spline, data.spline.tmp)
  }
  data.spline.all$line.type <- "dashed"
  data.spline.all$color <- 'red'
}
  
  
setEPS()
postscript(outfile, width=9, height=6)

p <- ggplot(data, aes(x=precursor.mass, y=probability, color=S))

print(
  p
  + geom_point(shape=1)
  + geom_line(data=data.spline, aes(x=precursor.mass, y=probability, group=S, linetype=line.type), color="black")
  + geom_line(data=data.spline.all, aes(x=precursor.mass, y=probability, linetype=line.type, color=color))
  + scale_linetype_identity()
  + scale_color_identity()
  )

dev.off()