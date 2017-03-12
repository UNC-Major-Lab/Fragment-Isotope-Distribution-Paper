library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
infile <- args[1]
otufile <- args[2]

data <- read.table(infile, sep="\t", header=T)

setEPS()
postscript(outfile, width=9, height=6)

print (
  ggplot(data=data, aes(x=isotope, y=time, 
                        group=method, color=method,
                        facet=comparison))
  + geom_line()
  + geom_point()
  + facet_wrap(~ comparison, scale="free", ncol=1)
  + xlab("max_depth")
  + ylab("runtime (ms)")
  + scale_x_continuous(breaks=1:11)
  + scale_y_log10(breaks=c(100,200,400,800,1600,3200,6400,12800))
)

dev.off()