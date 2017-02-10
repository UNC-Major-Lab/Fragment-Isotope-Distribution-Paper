#!/usr/bin/env Rscript
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

infile <- args[1]
outfile <- args[2]
binWidth <- as.numeric(args[3])
fragments <- args[4]

data <- read.table(infile,header=F,sep="\t")

if (fragments == "F") {

    p <- ggplot(data=data, aes(x=V2, y=V3, group=V1, fill=factor(V1)))

    pdf(outfile, width=9, height=6)

    print(p + geom_bar(stat='identity', width=binWidth) + facet_wrap(~ V1, scale="free") + xlab("residual") + ylab("frequency") + ggtitle("Histogram of residuals (averagine - expected)"))

    dev.off()

} else {

    p <- ggplot(data=data, aes(x=V2, y=V4, group=V1, fill=factor(V1)))

    pdf(outfile, width=9, height=20)

    print(p + geom_bar(stat='identity', width=binWidth) + facet_wrap(V3 ~ V1, scale="free") + xlab("residual") + ylab("frequency") + ggtitle("Histogram of residuals (averagine - expected)"))

    dev.off()

}
