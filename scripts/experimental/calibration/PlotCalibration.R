library(ggplot2)

savePlot <- function(myplot, outfile) {
  pdf(outfile, width=50, height=40)
  print(myplot)
  dev.off()
}

args = commandArgs(trailingOnly=TRUE)
infile1 <- args[1]
outPath <- args[2]

data <- read.table(infile1, header=T, sep="\t")

p <- ggplot(data=data, aes(x=offset, y=intensity)) +
  geom_point() +
  #+ geom_smooth(method='loess', aes(color=as.factor(target)))
  facet_wrap(width ~ target, scale="free")


savePlot(p, outPath)

