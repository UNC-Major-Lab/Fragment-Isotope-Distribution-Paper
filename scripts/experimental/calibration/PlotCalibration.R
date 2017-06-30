library(ggplot2)

savePlot <- function(myplot, outfile) {
  pdf(outfile, width=50, height=40)
  print(myplot)
  dev.off()
}

args = c("/Users/dennisg/Documents/Manuscripts/IsolationCalibration/OT_Quad.out",
         "/Users/dennisg/Documents/Manuscripts/IsolationCalibration/OT_Quad.pdf")
#args = commandArgs(trailingOnly=TRUE)
infile1 <- args[1]
outPath <- args[2]

data <- read.table(infile1, header=T, sep="\t")

p <- ggplot(data=data, aes(x=offset, y=intensity)) +
  geom_point() +
  geom_smooth(size = 2, method = "lm", formula = y ~ splines::bs(x, df=21), se = FALSE, aes(color=as.factor(target))) +
  geom_vline(size = 1.5, aes(xintercept = -width/2)) +
  geom_vline(size = 1.5, aes(xintercept = width/2)) +
  facet_wrap(width ~ target, scale="free")


savePlot(p, outPath)

