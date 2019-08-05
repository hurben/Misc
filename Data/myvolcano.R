#install.packages("ggplot2")
#install.packages("gridExtra")
#install.packages("plotly")
#install.packages("ggrepel")

library("ggplot2")
library("gridExtra")
library("plotly")
library("ggrepel")

data <- read.table(file.choose(), header=TRUE, row.names=1, sep="\t")

#data[,3] == y axis
#data[,4] == x axis

y_axis = data[,3]
x_axis = data[,4]
gene = row.names(data)
bcl6_data = data["Bcl6",]
bcl6_logpval = bcl6_data[3]
bcl6_log2fc = bcl6_data[4]

png(file = "EAT.R.KOWT.png", res = 100)
#png(file = "EAT.WT.RF.png", res = 100)

plot(x_axis , y_axis, 
col = ifelse(
		y_axis >= 3, 

		ifelse(
			x_axis >= 0.2 , "firebrick3", "dodgerblue3"
			),

		 "black"
		),
pch=20,
text(bcl6_log2fc , bcl6_logpval, "Bcl6", col="red", cex=2) ,
xlim=c(-7,7), ylim=c(0,18), xlab="Log2FC", ylab="-log2(P-value)")
title(main = "EAT R (KO/WT)")
#title(main = "EAT WT (R/F)")
dev.off()



png(file = "EAT.R.KOWT.png", res = 100)
#png(file = "EAT.WT.RF.png", res = 100)

plot(x_axis , y_axis, 
pch=20,
text(bcl6_log2fc , bcl6_logpval, "Bcl6", col="red", cex=2) ,
xlim=c(-7,7), ylim=c(0,18), xlab="Log2FC", ylab="-log2(P-value)")
points( , 
title(main = "EAT R (KO/WT)")
#title(main = "EAT WT (R/F)")

dev.off()




data <- read.table(file.choose(), header=TRUE, row.names=1, sep="\t")
data_matrix <- data
data_matrix$bcl <- "NO"
data_matrix$bcl[rownames(data_matrix) %in% "Bcl6"] <- "YES"
data_matrix["Bcl6",]

p <- ggplot( data_matrix , aes( x = log2FC, y = logpval))
p
#p <- p + xlim(-10,10)
p <- p + xlim(-7,7)
p
p1 <- p + geom_point(color = "grey30")
p1
p2 <- p1 + geom_point( data = data["Bcl6",] , color = "red", size = 3 )
p2
p3 <- p2 + geom_label_repel(aes ( label= ifelse(data_matrix$bcl == "YES", "Bcl6", ""), color="red"), show.legend = FALSE)
p3
#p4 <- p3 + ggtitle("EAT R (KO/WT)")
p4 <- p3 + ggtitle("EAT WT (R/F)")
#png(file = "EAT.R.KOWT.png", res = 100)
png(file = "EAT.WT.RF.png", res = 100)
p4
dev.off()

