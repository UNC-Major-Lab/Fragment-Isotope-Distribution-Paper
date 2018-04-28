#load libraries
library(ggplot2)
library(reshape2)
library(plyr)
library(xtable)

savePlot <- function(myplot, outfile) {
  #setEPS()
  #postscript(outfile)
  pdf(outfile, width=11, height=7)
  print(myplot)
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
X2headers <- c("exactCondFragmentX2","approxFragmentFromWeightX2","approxFragmentFromWeightAndSX2","exactPrecursorX2","approxPrecursorX2","approxFragmentSplineFromWeightX2","approxFragmentSplineFromWeightAndSulfurX2")

#loop through each search depth
for (searchDepth in 2:4) {
  #truncate to distributions at search depth
  distAtDepth <- distData[which(distData$completeAtDepth == searchDepth),]
  #truncate to distributions at search depth but have complete distributions
  distAtDepthComplete <- distAtDepth[which(distAtDepth$completeFlag == 1),]
  #distAtDepth <- distAtDepth[which(distAtDepth$completeFlag == 0),]
  #distAtDepth <- distAtDepth[which(distAtDepth$searchDepth == searchDepth+1),]
  distAtDepth <- distData[which(distData$searchDepth == searchDepth),]

  #melt for chi squared data
  meltedX2_AtDepth <- melt(distAtDepth, id.vars=c("ionID","scanDesc","distributionMonoWeight","precursorMonoWeight"), measure.vars=X2headers)
  meltedX2_Complete <- melt(distAtDepthComplete, id.vars=c("ionID","scanDesc","distributionMonoWeight","precursorMonoWeight"), measure.vars=X2headers)

  #rename melted factors for Chi Squared data
  meltedX2_AtDepth$variable <- as.character(meltedX2_AtDepth$variable)
  meltedX2_Complete$variable <- as.character(meltedX2_Complete$variable)

  meltedX2_AtDepth$variable[meltedX2_AtDepth$variable == "exactCondFragmentX2"] <- "ExactFragment"
  meltedX2_Complete$variable[meltedX2_Complete$variable == "exactCondFragmentX2"] <- "ExactFragment"

  meltedX2_AtDepth$variable[meltedX2_AtDepth$variable == "approxFragmentFromWeightX2"] <- "ApproxFragment"
  meltedX2_Complete$variable[meltedX2_Complete$variable == "approxFragmentFromWeightX2"] <- "ApproxFragment"

  meltedX2_AtDepth$variable[meltedX2_AtDepth$variable == "approxFragmentFromWeightAndSX2"] <- "ApproxFragmentSulf"
  meltedX2_Complete$variable[meltedX2_Complete$variable == "approxFragmentFromWeightAndSX2"] <- "ApproxFragmentSulf"

  meltedX2_AtDepth$variable[meltedX2_AtDepth$variable == "exactPrecursorX2"] <- "ExactPrecursor"
  meltedX2_Complete$variable[meltedX2_Complete$variable == "exactPrecursorX2"] <- "ExactPrecursor"

  meltedX2_AtDepth$variable[meltedX2_AtDepth$variable == "approxPrecursorX2"] <- "ApproxPrecursor"
  meltedX2_Complete$variable[meltedX2_Complete$variable == "approxPrecursorX2"] <- "ApproxPrecursor"

  meltedX2_AtDepth$variable[meltedX2_AtDepth$variable == "approxFragmentSplineFromWeightX2"] <- "splineFragment"
  meltedX2_Complete$variable[meltedX2_Complete$variable == "approxFragmentSplineFromWeightX2"] <- "splineFragment"

  meltedX2_AtDepth$variable[meltedX2_AtDepth$variable == "approxFragmentSplineFromWeightAndSulfurX2"] <- "splineFragmentSulf"
  meltedX2_Complete$variable[meltedX2_Complete$variable == "approxFragmentSplineFromWeightAndSulfurX2"] <- "splineFragmentSulf"

  meltedX2_AtDepth$variable <- as.factor(meltedX2_AtDepth$variable)
  meltedX2_Complete$variable <- as.factor(meltedX2_Complete$variable)
  #plot chi squared density
  plotX2_atDepth <- ggplot(data = meltedX2_AtDepth, mapping = aes(x=value, color=variable)) +
    geom_density(adjust = 3) +
    #scale_x_log10(limits=c(1e-30,1e10)) +
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
    facet_wrap(~ scanDesc) +#, scales="free_y") +
    scale_x_log10() +
    theme(legend.position = "right",
          legend.text = element_text(size = 10),
          axis.title = element_text(size=12),
          plot.title = element_text(size=16)) +
    ggtitle(paste("Chi-Squared Complete at Depth: ", searchDepth, sep = "")) +
    ylab("Density") +
    xlab(expression(X^{2}))

  #plot chi squared density
  data.subset <- subset(meltedX2_Complete, meltedX2_Complete$variable == "ExactFragment")
  plotX2_Complete_vs_mass <- ggplot(data = data.subset, mapping = aes(y=value, x=precursorMonoWeight-distributionMonoWeight)) +
    #geom_point(alpha = 0.01) +
    geom_smooth(method="gam") +
    #scale_y_log10() +
    facet_wrap(~ scanDesc) +#, scales="free_y") +
    theme(legend.position = "right",
          legend.text = element_text(size = 10),
          axis.title = element_text(size=12),
          plot.title = element_text(size=16)) +
    ggtitle(paste("Chi-Squared Complete at Depth: ", searchDepth, sep = "")) +
    xlab("Mass") +
    ylab(expression(X^{2}))

  savePlot(plotX2_Complete_vs_mass, paste(outPath,sep = "", "/chi-squared_complete_vs_mass_",searchDepth, ".pdf"))

  savePlot(plotX2_atDepth, paste(outPath,sep = "", "/chi-squared_incomplete_",searchDepth, ".pdf"))
  out.table <- ddply(meltedX2_AtDepth, "variable", summarize, median=median(value), min=min(value), max=max(value), sd=sd(value), count=length(value))
  #write.table(out.table, sep="\t", file=paste(outPath, sep = "", "/chi-squared_incomplete_",searchDepth, ".tab"))
  print(xtable(out.table))

  savePlot(plotX2_Complete, paste(outPath, sep = "", "/chi-squared_complete_",searchDepth, ".pdf"))

  CID.data <- subset(meltedX2_AtDepth, meltedX2_AtDepth$scanDesc == "CID_25")
  #CID.data <- subset(meltedX2_Complete, meltedX2_Complete$scanDesc == "CID_25")
  out.table <- ddply(CID.data, "variable", summarize, median=median(value), mean=mean(value), count=length(value))
  #write.table(out.table, sep="\t", file=paste(outPath, sep = "", "/chi-squared_complete_",searchDepth, ".tab"))
  #print(xtable(out.table, digits=4))

  #HCD.data <- subset(meltedX2_Complete, meltedX2_Complete$scanDesc == "HCD")
  #out.table <- ddply(HCD.data, "variable", summarize, median=median(value), mean=mean(value), count=length(value))
  #write.table(out.table, sep="\t", file=paste(outPath, sep = "", "/chi-squared_complete_",searchDepth, ".tab"))
  #print(xtable(out.table, digits=4))


}
