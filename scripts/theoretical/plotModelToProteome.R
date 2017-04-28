#!/usr/bin/env Rscript
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

spline.infile <- args[1]
data.basedir <- args[2]
precursor <- args[3]
max.sulfur <- args[4]
outfile <- args[5]
max.mass <- as.numeric(args[6])


# Proteome theoretical data
scatter.infile <- paste(data.basedir, "training/proteome/Precursor", precursor, ".tab", sep="")
data <- read.table(scatter.infile, header=T, sep="\t")
data <- subset(data, data$precursor.mass <= max.mass)
data <- subset(data, as.numeric(data$sulfur) <= max.sulfur)
#data$sulfur[which(data$sulfur > max.sulfur)] <- paste(">",toString(max.sulfur),"")

# Proteome averagine eval
averagine.infile <- paste(data.basedir, "training/proteome/averagine/Precursor", precursor, ".tab", sep="")
data.averagine <- read.table(averagine.infile, header=T, sep="\t")
data.averagine <- subset(data.averagine, data.averagine$precursor.mass <= max.mass)

# Average spline
spline.infile <- paste(data.basedir, "Averagine/Average_Spline", "/spline/eval/Precursor", precursor, ".tab", sep="")
data.spline.avg <- read.table(spline.infile, header=T, sep="\t")
data.spline.avg <- subset(data.spline.avg, data.spline.avg$precursor.mass <= max.mass)

# Sulfur-specific splines
data.spline <- data.frame()
for (sulfur in 0:max.sulfur) {
  spline.infile <- paste(data.basedir, "Averagine/S", toString(sulfur), "/spline/eval/Precursor", precursor, ".tab", sep="")
  data.spline.tmp <- read.table(spline.infile, header=T, sep="\t")
  data.spline.tmp <- subset(data.spline.tmp, data.spline.tmp$precursor.mass <= max.mass)
  data.spline.tmp$S <- toString(sulfur)
  data.spline <- rbind(data.spline, data.spline.tmp)
}
  
  
setEPS()
postscript(outfile, width=9, height=6)

p <- ggplot()

print(
  p
  + geom_point(data=data, aes(x=precursor.mass, y=probability, color=as.factor(sulfur)), shape=1)
  # sulfur-specific splines
  + geom_line(data=data.spline, aes(x=precursor.mass, y=probability, group=S), color="black")
  # average spline
  + geom_line(data=data.spline.avg, aes(x=precursor.mass, y=probability), color="red", linetype="dashed")
  # averagine model
  + geom_line(data=data.averagine, aes(x=precursor.mass, y=probability), color="blue", linetype="dashed")
  )

dev.off()