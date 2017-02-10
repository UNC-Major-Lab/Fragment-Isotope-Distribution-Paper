#!/usr/bin/env Rscript
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

infile <- args[1]
outfile <- args[2]
binWidth <- as.numeric(args[3])
fragments <- args[4]

data <- read.table(infile,header=F,sep="\t")

if (fragments == "T") {

    p <- ggplot(data=data, aes(x=V2, y=V3, group=V1, fill=factor(V1)))

    pdf(outfile, width=9, height=6)

    p + geom_bar(stat='identity', width=binWidth) + facet_wrap(~ V1) + xlab("residual") + ylab("frequency") + ggtitle("Histogram of residuals (averagine - expected)")

    dev.off()

} else {

    p <- ggplot(data=data, aes(x=V2, y=V4, group=V1, fill=factor(V1)))

    pdf(outfile, width=9, height=6)

    p + geom_bar(stat='identity', width=binWidth) + facet_wrap(V3 ~ V1) + xlab("residual") + ylab("frequency") + ggtitle("Histogram of residuals (averagine - expected)")

    dev.off()

}
