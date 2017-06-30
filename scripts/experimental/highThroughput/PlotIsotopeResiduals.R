library(ggplot2)
library(reshape2)
library(plyr)

args = c("/Users/dennisg/Documents/Manuscripts/FragmentIsotopicDistribution/data/LT/Angiotensin_z3_HCD_isotopes.tab",
                  "/Users/dennisg/Downloads/")
#args = c("/Users/dennisg/Downloads/isotopesScores.out",
#         "/Users/dennisg/Downloads/")
#args = commandArgs(trailingOnly=TRUE)
infile1 <- args[1]
outPath <- args[2]



#load data from file
data <- read.table(infile1, header=T, sep="\t")


#data <- subset(data, data$monoMass < 500)

data.melted <- melt(data, id.vars=c("searchDepth","scanDesc","isotope","monoMz","monoMass","precursorMass","intensity"), variable="method", value.name="residual")
data.melted <- subset(data.melted, data.melted$method == "residualCalibratedExactFragment")



p <- ggplot(data=data.melted, aes(x=as.factor(isotope),
                                  y=residual,
                                  color=as.factor(isotope)))

print(
  p
  + geom_point(position="jitter", alpha=0.2)
  + geom_boxplot(outlier.size=NA)
  + facet_grid(method ~ searchDepth + scanDesc)
  #+ scale_y_continuous(limits=c(-.15, .15))
)

#p <- ggplot(data=data, aes(x=precursorMass-monoMass, y=residualExactFragment))
#
#print(
#  p
#  + geom_point(alpha=0.05)
#  + facet_grid(isotope ~ searchDepth + isSIM)
#)

#p <- ggplot(data=data, aes(x=monoMass, y=residualExactFragment))

#print(
#  p
#  + geom_point(alpha=0.05)
#  + facet_grid(isotope ~ searchDepth + scanDesc)
#)


