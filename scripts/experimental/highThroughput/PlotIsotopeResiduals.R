library(ggplot2)
library(reshape2)
library(plyr)

args = c("/Users/dennisg/Downloads/isotopesScores.out",
         "/Users/dennisg/Downloads/")
#args = commandArgs(trailingOnly=TRUE)
infile1 <- args[1]
outPath <- args[2]



#load data from file
data <- read.table(infile1, header=T, sep="\t")

data.melted <- melt(data, id.vars=c("searchDepth","isSIM","isotope","monoMz","monoMass","intensity"), variable="method", value.name="residual")

p <- ggplot(data=data.melted, aes(x=as.factor(isotope), 
                                  y=residual, 
                                  color=as.factor(isotope)))

print(
  p
  + geom_boxplot(outlier.size=NA)
  + facet_grid(method ~ searchDepth + isSIM)
  + scale_y_continuous(limits=c(-0.25, 0.25))
)

#p <- ggplot(data=data, aes(x=monoMz, y=residualExactFragment))
#
#print(
#  p
#  + geom_point()
#  + facet_grid(isotope ~ searchDepth)
#)
