library(rtracklayer)
gtf1 <- rtracklayer::import('~/backup/Gencode/Release22/gencode.v22.long_noncoding_RNAs.gtf.gz')
gtf1.df <- as.data.frame(gtf1)
gene1.df  <- gtf1.df[gtf1.df$type=="gene",]
tmp.idx <- sort(c(grep("sense", gene1.df$gene_type), which(gene1.df$gene_type=="lincRNA"), which(gene1.df$gene_type=="TEC") ) )
LNC.df <- gene1.df[tmp.idx,]


cancers.v  <-  c("BLCA", "BRCA(ER+)", "BRCA(ER-)", "CESC", "CHOL", "COAD", "ESCA", "HNSC",
				 "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PAAD"  , "PRAD", "READ",  "THCA", "UCEC")


LncStat.m <- matrix(NA, nrow=length(cancers.v), ncol=5)
rownames(LncStat.m)  <-  cancers.v
colnames(LncStat.m)  <-  c("antisense", "lincRNA", "sense_intronic", "sense_overlapping", "TEC")

i <- 1
setwd("/mnt/b/lnc_DNAm/CancerRdata/LncExpStat/LncValExp0")
load(file=paste0(cancers.v[i], "_valFPKM.Rdata") )
ConsistentLnc.v <- rownames(LncValFPKM.m)
ConsistentLncAll.v <- rownames(LncValFPKM.m)

for(i in 1:length(cancers.v)){
	print(i)
	setwd("/mnt/b/lnc_DNAm/CancerRdata/LncExpStat/LncValExp0")
	load(file=paste0(cancers.v[i], "_valFPKM.Rdata") )
	tmp.idx <- match(rownames(LncValFPKM.m), LNC.df$gene_id)
	lncTypesStat.v <- table(LNC.df$gene_type[tmp.idx])
	tmp.idx <- match(names(lncTypesStat.v), colnames(LncStat.m))
	LncStat.m[i,] <- lncTypesStat.v[tmp.idx]
	ConsistentLnc.v <- intersect(ConsistentLnc.v , rownames(LncValFPKM.m))
	ConsistentLncAll.v <- union(ConsistentLncAll.v , rownames(LncValFPKM.m))
}

setwd("/mnt/b/lnc_DNAm/CancerRdata/LncExpStat/")
write.table(LncStat.m, file="LncTypeStat0.csv", sep="\t", quote=F, col.names=NA)
save(ConsistentLnc.v, file="ConsistentLnc0.Rdata")


cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
barplot(t(LncStat.m), ylim=c(1,6000), col=cbPalette[3:7])
LncStat.df <- as.data.frame(LncStat.m)

library("reshape2")
LncStat.df <- melt(LncStat.m )
colnames(LncStat.df) <- c("Cancer", "Type", "Num")
library("ggplot2")
library("extrafont")
font_import()
loadfonts()   
fonts()  
p <- ggplot(LncStat.df, aes(x = Cancer, y = Num, fill = Type)) + geom_bar(stat ="identity") + theme_bw() + theme(axis.text.x = element_text(angle=30, hjust=1,vjust=1), text=element_text(family="Arial"))
p 	
ggsave("LncTypeStat0.pdf",device=cairo_pdf, width=10,height=4)

library("RColorBrewer")
p <- ggplot(LncStat.df, aes(x = Cancer, y = Num, fill = Type)) + geom_bar(stat ="identity") + theme_bw() + theme(axis.text.x = element_text(angle=30, hjust=1,vjust=1), text=element_text(family="Arial")) + scale_fill_brewer(palette="Dark2") +  geom_text(aes(label = Num), size = 1, hjust = 0.2, vjust = 1.5, position = "stack") 
p
ggsave("LncTypeStat0_1_annNum.pdf",device=cairo_pdf, width=10,height=4)

###
p <- ggplot(LncStat.df, aes(x = Cancer, y = Num, fill = Type)) + geom_bar(stat ="identity") + theme_bw() + theme(axis.text.x = element_text(angle=30, hjust=1,vjust=1), text=element_text(family="Arial")) + scale_fill_brewer(palette="Accent")
p
ggsave("LncTypeStat0_2.pdf",device=cairo_pdf, width=10,height=4)
###
p <- ggplot(LncStat.df, aes(x = Cancer, y = Num, fill = Type)) + geom_bar(stat ="identity") + theme_bw() + theme(axis.text.x = element_text(angle=30, hjust=1,vjust=1), text=element_text(family="Arial")) + scale_fill_brewer(palette="Set1")
p
ggsave("LncTypeStat0_3.pdf",device=cairo_pdf, width=10,height=4)


ce <- arrange(LncStat.df, Cancer )



###########
tmp.idx <- match(ConsistentLncAll.v , LNC.df$gene_id)
lncTypesStat.df <- table(LNC.df$gene_type[tmp.idx])
barplot(t(lncTypesStat.df), col=cbPalette[3:7])

slicesExp <- as.data.frame(lncTypesStat.df)
#names(slicesExp) <- colnames(LncStat.m);
#lbls <- paste(slicesExp[,1], "\n(n = ", slicesExp[,2], ")", sep= "");
lbls <- paste( "n = ", slicesExp[,2], sep= "");

pdf("Figure_1C.pdf",width=3, height=3,family="Arial")
#par(family="Arial")
par(mar=c(0.28,0.05,0.1,0.52))
#pie(as.numeric(slicesExp[,2]), labels=lbls, col=brewer.pal(5,"Set1"), clockwise=TRUE, radius = 0.55, cex= 0.5, font = 2)
pie(as.numeric(slicesExp[,2]), labels=lbls, col=brewer.pal(5,"Set1"), clockwise=TRUE, radius = 0.55, cex= 0.5, font = 2)
legend(x=0.45, y=0.70, as.vector(slicesExp[,1]), pch=15, bty="n", col=brewer.pal(5,"Set1"), cex=0.5)
dev.off()


pdf("Figure_1C_1.pdf",width=3, height=3,family="Arial")
#par(family="Arial")
par(mar=c(0.28,0.05,0.1,0.52))
#pie(as.numeric(slicesExp[,2]), labels=lbls, col=brewer.pal(5,"Set1"), clockwise=TRUE, radius = 0.55, cex= 0.5, font = 2)
pie(as.numeric(slicesExp[,2]), labels=lbls, col=brewer.pal(5,"Dark2"), clockwise=TRUE, radius = 0.55, cex= 0.5, font = 2)
legend(x=0.45, y=0.70, as.vector(slicesExp[,1]), pch=15, bty="n", col=brewer.pal(5,"Dark2"), cex=0.5)
dev.off()


