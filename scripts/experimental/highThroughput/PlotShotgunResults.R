#load libraries
library(ggplot2)
library(reshape2)
library(plyr)

savePlot <- function(myplot, outfile) {
  setEPS()
  postscript(outfile)

  print(plotX2_atDepth)
  dev.off()
}

args = commandArgs(trailingOnly=TRUE)
infile1 <- args[1]
outPath <- args[2]
#load data from file
distData <- read.table(infile1, header=T, sep="\t")
#only distributions flagged as valid
distData <- distData[which(distData$distributionValid == 1), ]

#Chi-squared file headers
X2headers <- c("exactCondFragmentX2","approxFragmentFromWeightX2","approxFragmentFromWeightAndSX2")

#loop through each search depth
for (searchDepth in 2:5) {
  #truncate to distributions at search depth
  distAtDepth <- distData[which(distData$completeAtDepth == searchDepth),]
  #truncate to distributions at search depth but have complete distributions
  distAtDepthComplete <- distAtDepth[which(distAtDepth$completeFlag == 1),]
  distAtDepth <- distAtDepth[which(distAtDepth$completeFlag == 0),]
  distAtDepth <- distAtDepth[which(distAtDepth$searchDepth == searchDepth+1),]
  
  #melt for chi squared data
  meltedX2_AtDepth <- melt(distAtDepth, id.vars="ionID", measure.vars=X2headers)
  meltedX2_Complete <- melt(distAtDepthComplete, id.vars="ionID", measure.vars=X2headers)
  
  #rename melted factors for Chi Squared data
  meltedX2_AtDepth$variable <- as.character(meltedX2_AtDepth$variable)
  meltedX2_Complete$variable <- as.character(meltedX2_Complete$variable)
  meltedX2_AtDepth$variable[meltedX2_AtDepth$variable == "exactCondFragmentX2"] <- "ExactFragment"
  meltedX2_Complete$variable[meltedX2_Complete$variable == "exactCondFragmentX2"] <- "ExactFragment"
  meltedX2_AtDepth$variable[meltedX2_AtDepth$variable == "approxFragmentFromWeightX2"] <- "ApproxFragment"
  meltedX2_Complete$variable[meltedX2_Complete$variable == "approxFragmentFromWeightX2"] <- "ApproxFragment"
  meltedX2_AtDepth$variable[meltedX2_AtDepth$variable == "approxFragmentFromWeightAndSX2"] <- "ApproxFragmentSulf"
  meltedX2_Complete$variable[meltedX2_Complete$variable == "approxFragmentFromWeightAndSX2"] <- "ApproxFragmentSulf"
  meltedX2_AtDepth$variable <- as.factor(meltedX2_AtDepth$variable)
  meltedX2_Complete$variable <- as.factor(meltedX2_Complete$variable)
  #plot chi squared density
  plotX2_atDepth <- ggplot(data = meltedX2_AtDepth, mapping = aes(x=value, color=variable)) +
    geom_density() +
    scale_x_log10() +
    theme(legend.position = "right",
          legend.text = element_text(size = 10),
          axis.title = element_text(size=12),
          plot.title = element_text(size=16)) +
    ggtitle(paste("Chi-Squared at Depth: ", searchDepth, sep = "")) +
    ylab("Density") +
    xlab(expression(X^{2}))


  #plot chi squared density
  plotX2_Complete <- ggplot(data = meltedX2_Complete, mapping = aes(x=value, color=variable)) +
    geom_density() +
    scale_x_log10() +
    theme(legend.position = "right",
          legend.text = element_text(size = 10),
          axis.title = element_text(size=12),
          plot.title = element_text(size=16)) +
    ggtitle(paste("Chi-Squared Complete at Depth: ", searchDepth, sep = "")) +
    ylab("Density") +
    xlab(expression(X^{2}))

  savePlot(plotX2_atDepth, paste(outPath,sep = "", "/chi-squared_incomplete_",searchDepth, ".eps"))
  out.table <- ddply(meltedX2_AtDepth, "variable", summarize, median=median(value), min=min(value), max=max(value), sd=sd(value), count=length(value))
  write.table(out.table, sep="\t", file=paste(outPath, sep = "", "/chi-squared_incomplete_",searchDepth, ".tab"))
  savePlot(plotX2_Complete, paste(outPath, sep = "", "/chi-squared_complete_",searchDepth, ".eps"))
  out.table <- ddply(meltedX2_Complete, "variable", summarize, median=median(value), min=min(value), max=max(value), sd=sd(value), count=length(value))
  write.table(out.table, sep="\t", file=paste(outPath, sep = "", "/chi-squared_complete_",searchDepth, ".eps"))
}
