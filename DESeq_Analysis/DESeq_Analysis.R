# DESeq2 Run
#	DESeq is to detect differentially expressed genes (DEGs) from RNA-Seq data
#	Expression data must be in read counts (NOT Normalized data such as FPKM or RPKM)

# 1. Installation
#	source from bioconductor
source( "http://bioconductor.org/biocLite.R" )

# 1-1. DESeq
biocLite( "DESeq2" )
library( DESeq2 )

# 2. Data Import
EATCountTable	=	read.table( file.choose(), sep="\t", head=TRUE, row.names=1 )
				# WT_R vs WT_F:	EAT.WT.matrix.gene
				# KO_R vs KO_F:	EAT.KO.matrix.gene
				# HES vs NEA:	HEA.NEA.matrix.gene
temp			=	apply( EATCountTable, 2, as.integer )
rownames(temp)	=	rownames( EATCountTable )
EATCountTable	=	temp

EATDesign		=	data.frame(
					row.names	=	colnames( EATCountTable ),
					condition	=	rep( c('refeeding', 'feeding'), each=2)	# WT_R vs WT_F
#					condition	=	rep( c('refeeding', 'feeding'), each=2)	# KO_R vs KO_F
#					condition	=	rep( c('HEA', 'NEA'), each=3)			# HEA vs NEA
				)
EATDesign

countTable		=	EATCountTable
condition		=	EATDesign$condition

#	3-1. Check for the result
head( countTable )
condition

#	3-2. Instantiation of countDataSet (=cds)
#cds	=	newCountDataSet( countTable, condition )		# DESeq
cds	=	DESeqDataSetFromMatrix(	countData=countTable, 
						colData=EATDesign,
						design=~condition )	# DESeq2
#	3-3. Normalization (consideration of the library size)
#		To estimate the effective library size
cds	=	estimateSizeFactors( cds )
sizeFactors( cds )
#		Apply the size factor to the original counts data
head ( counts ( cds, normalized=TRUE ) )

# 4. Variance estimation
########################
#	Variance = (sample-to-sample variance) + (uncertainty in measurement)
#			(= dispersion)			(= shot/Poisson noise)
#		which dominates for...
#			highly expr genes			lowly expr genes
########################
#	4-1. To estimate the dispersions,
cds	=	estimateDispersions( cds )
#		For detailed information of estimateDispersions function,
plotDispEsts(cds)		# Plot

# 5. Inference: Calling differential expression
#	5-1. Standard comparison between two experimental conditions
cds	=	DESeq( cds )
res	=	results( cds, pAdjustMethod='fdr' )
#		Interpretation:
#			1) id				feature identifier
#			2) baseMean			mean normed cnts, av over all samples from both cdns
#			3) baseMeanA		mean normed cnts from cdn A
#			4) baseMeanB		mean normed cnts from cdn B
#			5) foldChange		fold change from condition A to B
#			6) log2FoldChange	the logarithm (to basis 2) of the fold change
#			7) pval				p value for the statistical significance of this change
#			8) padj				p.value adjusted for multiple testing to control FDR
#		To plot the result,
plotMA(res)
hist(res$pvalue, breaks=100, col="skyblue", border="slateblue", main="")



# Volcano plot
head(res)
print( paste( "The number of NA GENEs:", sum(is.na(res$padj)) ) )
res2	<-	res[!is.na(res$padj),]
dev.new(width=15,height=15)
par(mar=c(6,6,1,1))
plot( -log10(res2$pvalue) ~ res2$log2FoldChange, pch=16,
	xlim = c(-5,5), ylim = c(0,15), font.axis=2, font.lab=2, cex.lab = 2, cex.axis = 1.5,
	col = ifelse(-log10(res2$padj)>2, ifelse( res2$log2FoldChange>=0, "red", "blue"), "black"),
	cex = ifelse(-log10(res2$padj)>2, 0.5, 0.3),
	xlab="log2(FC)", ylab="-log(p-value)" )
box(lwd=4)

# Volcano Plot for log2FC < 0 
dev.new(width=8,height=15)
par(mar=c(5,5,1,1))
plot( -log10(res2$pvalue) ~ res2$log2FoldChange, pch=16,
	xlim = c(-5,0), ylim = c(0,15), font.axis=2, font.lab=2, cex.lab = 2, cex.axis = 1.5,
	col = ifelse(-log10(res2$padj)>2, ifelse( res2$log2FoldChange >=0, "red", "blue"), "black"),
	cex = ifelse(-log10(res2$padj)>2, 1, 0.3),
	xlab="log2(FC)", ylab="-log(p-value)" )
box(lwd=4)

# Volcano Plot for log2FC > 0
par(mar=c(5,3,1,1)) 
plot( -log10(res2$pvalue) ~ res2$log2FoldChange, pch=16,
	xlim = c(0,5), ylim = c(0,15), font.axis=2, font.lab=2, cex.lab = 2, cex.axis = 1.5,
	col = ifelse(-log10(res2$padj)>2, ifelse( res2$log2FoldChange>=0, "red", "blue"), "black"),
	cex = ifelse(-log10(res2$padj)>2, 1, 0.3),
	xlab="log2(FC)", ylab="" )
box(lwd=4)


#		To filter for significant genes,
isDEG		=	((res$padj < 0.01) & (res$padj!="NA"))
res.final	=	cbind( res, "DESeq_DEG" = isDEG )
resSig	=	res[which( isDEG == TRUE),]
#		and list the most significant DEGs
head( resSig[ order( resSig$pvalue), ])
#		or to look at the most strongly up/down-regulated of the significant genes,
head( resSig[ order( resSig$pvalue ), ] )
head( resSig[ order( -resSig$log2FoldChange, -resSig$baseMean ), ] );	sum(resSig$log2FoldChange>0);
head( resSig[ order( resSig$log2FoldChange, -resSig$baseMean ), ] );	sum(resSig$log2FoldChange<0);
#		To save the output,
write.csv( res.final, file.choose() )
