library("RColorBrewer")
library(gplots)

par(mar=c(5,6,4,1))
data_file <- read.table(file.choose(), sep="\t", header=TRUE, row.names=1)
data_matrix <- data.matrix(data_file)


bk <- seq(from=-2, to=2, by=0.05)
color.palette  <- colorRampPalette(c("firebrick4", "white", "dodgerblue4"))
png(filename='test.png',width=1000,height=1000, res=100)
heatmap.2(data_matrix, col=color.palette, srtCol=30, breaks=bk,cexCol=1, trace="none")
graphics.off()


par(mar=c(5,5,1,1))

plot(1:1000,1:1000,col=redblue(1000),pty=2)