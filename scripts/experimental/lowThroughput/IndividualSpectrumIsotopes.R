library(ggplot2)

data <- read.table("/Users/dennisg/Documents/Manuscripts/FragmentIsotopicDistribution/data/LT/out04.tab", sep="\t", header=T)
data.calc <- read.table("/Users/dennisg/Documents/Manuscripts/FragmentIsotopicDistribution/data/LT/calc_out04.tab", sep="\t", header=T)

data <- subset(data, (data$ion.index == 67 | data$ion.index == 55) & data$isotope.range <= 5)
data.calc <- subset(data.calc, (data.calc$ion.index == 67 | data.calc$ion.index == 55) & data.calc$isotope.range <= 5)

p <- ggplot(data=data, aes(x=mz,y=int))

print(
  p
  + geom_line()
  + geom_point(data=data.calc, aes(group=method, color=method), size=3, shape=1)
  + facet_grid(isotope.range ~ ion.name, scale="free")
  + ylab("Intensity normalized to base peak")
  + xlab("m/z")
)