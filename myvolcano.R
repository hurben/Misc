library("ggplot2")
library("gridExtra")
library("plotly")
library("ggrepel")


data <- read.table(file.choose(), header=TRUE, row.names=1, sep="\t")
data_matrix <- data
#data_matrix$bcl6 <- "NO"
#data_matrix$bcl6[rownames(data_matrix) %in% "Bcl6"] <- "YES"

data_matrix$bcl6 <- "NO"
data_matrix[  which(  grepl("Bcl6_",rownames(data_matrix) )   ),]$bcl6 <- "YES"
data_matrix[which(  grepl("Bcl6_",rownames(data_matrix))),]


p <- ggplot( data_matrix , aes( x = LOG2FC, y = LOGPVAL))
p
#p <- p + xlim(-10,10)
p <- p + xlim(-2, 2) + ylim(0,10) + xlab("Log2 FoldChange") + ylab("-log10(PPEE)")
p
#p1 <- p + geom_point(color = "grey30")
p1 <- p + geom_point( data = data_matrix[data_matrix$LOGPVAL > 0,], color = "grey30")
p1
p1 <- p1 + geom_point( data = data_matrix[data_matrix$DEG == 1  & data$LOG2FC >= 0.1,], color="firebrick3")
p1 <- p1 + geom_point( data = data[data_matrix$DEG == 1 & data_matrix$LOG2FC <= -0.1,], color="dodgerblue3")
p1
p2 <- p1 + geom_point( data = data["Bcl6_",] , color = "orange2", size = 3 , shape=17)
p3 <- p2
#p3 <- p2 + geom_label_repel(aes ( label= ifelse(data_matrix$bcl == "YES", "Bcl6", ""), color="red"), show.legend = FALSE)
#p3 <- p2 + geom_label_repel(aes ( label= ifelse(data_matrix["Bcl6_",]$DEG == 1, "Bcl6", ""), color="red"), show.legend = FALSE)
#p3 <- p2 + geom_label_repel(aes ( label= ifelse(data_matrix$bcl6 == "YES", "Bcl6", ""), color="orange2"), show.legend = FALSE)
p3
p4 <- p3 + ggtitle("EAT R (KO/WT)")
p4 <- p3 + ggtitle("EAT WT (R/F)")
#png(file = "EAT.R.KOWT.png", res = 100)
#png(file = "EAT.WT.RF.png", res = 100)
p4
ggsave("EAT.R.KOWT.png")
ggsave("EAT.WT.RF.png")

