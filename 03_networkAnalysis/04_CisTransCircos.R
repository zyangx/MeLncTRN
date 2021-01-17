
options(stringsAsFactors=FALSE)
library("OmicCircos")
library("naturalsort")
library(org.Hs.eg.db)
library(stringr) 

require("RColorBrewer")
display.brewer.all()
brewer.pal(9, "Purples")
colors.v <- colorRampPalette(brewer.pal(9, "Greens"))(2)
colors.v <- colorRampPalette(brewer.pal(9, "BuPu"))(2)
colors.v <- colorRampPalette(brewer.pal(9, "BuGn"))(2)
colors.v <- colorRampPalette(brewer.pal(9, "Blues"))(2)


## Load cytoband data. Cytoband data downloaded from UCSC for hg38.
hg38.cyto <- read.csv("./Circos/hg38_cytoband.txt", strip.white = T, sep="\t", header=T)
#names(hg38.cyto) <- c("chrom","chromStart","chromEnd","name","gieStain")
# Order the chromosomes in human readable format
hg38.cyto <- hg38.cyto[naturalorder(hg38.cyto$chrom),]


seg.name <- paste("chr",c(1:22, "X", "Y"),sep="");
db <- segAnglePo(hg38.cyto,seg=seg.name);
colors.v <- rainbow(24,alpha=0.5);

cancers.v  <-  c("BLCA", "BRCA(ER+)", "BRCA(ER-)", "CESC", "CHOL", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC",
				"LUAD", "LUSC", "PAAD" , "PRAD", "READ", "THCA", "UCEC")

load("./Circos/gencode.Rdata")

EnsemblIDs.v <- str_split(PCG.df$gene_id,'[.]',simplify = T)[,1]
EntrezIDs.v <- mapIds(org.Hs.eg.db, keys=EnsemblIDs.v, column= 'ENTREZID', keytype="ENSEMBL")
PCG.df$EntrezID <- EntrezIDs.v 


pdf("Lnc_PCG_interact_Circos.pdf", width=10,height=13)
par(mfrow=c(5,4));

for(c in 1:length(cancers.v)) {
#rm(list=ls())
print(cancers.v[c])

load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
load(paste("./InteractNetworkStat1/", cancers.v[c], "_InteractNetStat.Rdata", sep="") )

LncPos.df <- as.data.frame(matrix(NA, nrow=nrow(InteractLG.m), ncol=3))
colnames(LncPos.df) <- c("chr", "po", "NAME")
LncPos.df[,3] <-  as.vector(rownames(InteractLG.m))
tmp.idx <- match(LncPos.df[,3], LNC.df$gene_id)
LncPos.df[,1] <- LNC.df$seqnames[tmp.idx]
LncPos.df[,2] <- LNC.df$start[tmp.idx]

PCGPos.df <- as.data.frame(matrix(NA, nrow=ncol(InteractLG.m), ncol=3))
colnames(PCGPos.df) <- c("chr", "po", "NAME")
PCGPos.df[,3] <-  as.vector(colnames(InteractLG.m))
tmp.idx <- match(PCGPos.df[,3], PCG.df$EntrezID)
PCGPos.df[,1] <- PCG.df$seqnames[tmp.idx]
PCGPos.df[,2] <- PCG.df$start[tmp.idx]

LncPCGLink.df <- as.data.frame(matrix(NA, nrow=nrow(validLncTargetStat.df), ncol=6))
colnames(LncPCGLink.df) <- c("Chr1", "Pos1", "NAME1", "Chr2", "Pos2", "NAME2")
LncPCGLink.df[,3] <-  validLncTargetStat.df[,1]
tmp.idx <- match(LncPCGLink.df[,3], LNC.df$gene_id)
LncPCGLink.df[,1] <- LNC.df$seqnames[tmp.idx]
LncPCGLink.df[,2] <- LNC.df$start[tmp.idx]

LncPCGLink.df[,6] <-  validLncTargetStat.df[,2]
tmp.idx <- match(LncPCGLink.df[,6], PCG.df$EntrezID)
LncPCGLink.df[,4] <- PCG.df$seqnames[tmp.idx]
LncPCGLink.df[,5] <- PCG.df$start[tmp.idx]

par(mar=c(1,1,1,1));
plot(c(1,800),c(1,800),type="n",axes=FALSE,xlab="",ylab="",main=cancers.v[c]);
par(mar=c(0.5, 0.5, 0.5, 0.5));
circos(R=400, cir=db, type="chr", col=colors.v, print.chr.lab=TRUE, W=10, scale=TRUE);
circos(R=360,cir=db,W=40,mapping=LncPos.df,type="b3",B=FALSE, col="red",lwd=1);
circos(R=320,cir=db,W=40,mapping=PCGPos.df,type="b3",B=FALSE,col="blue",lwd=1);
circos(R=300, cir=db, W=20, mapping=LncPCGLink.df, type="link", lwd=0.1, col="gray");


}

dev.off()

################
################
library (plotrix)

AllLncPCGLink.df  <- as.data.frame(matrix(NA, nrow=1, ncol=6))
colnames(AllLncPCGLink.df) <- c("Chr1", "Pos1", "NAME1", "Chr2", "Pos2", "NAME2")

#pcom.lv <- list()

pdf("All_Interact_CisTrans_barplot.pdf", width=10,height=13);
par(mfrow=c(5,4))
par(mar=c(4.5,4,3,1))

for(c in 1:length(cancers.v)) {
#rm(list=ls())
print(cancers.v[c])

load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
load(paste("./InteractNetworkStat1/", cancers.v[c], "_InteractNetStat.Rdata", sep="") )


LncPCGLink.df <- as.data.frame(matrix(NA, nrow=nrow(validLncTargetStat.df), ncol=6))
colnames(LncPCGLink.df) <- c("Chr1", "Pos1", "NAME1", "Chr2", "Pos2", "NAME2")
LncPCGLink.df[,3] <-  as.vector(validLncTargetStat.df[,1])
tmp.idx <- match(LncPCGLink.df[,3], LNC.df$gene_id)
LncPCGLink.df[,1] <- as.vector(LNC.df$seqnames[tmp.idx])
LncPCGLink.df[,2] <- LNC.df$start[tmp.idx]

LncPCGLink.df[,6] <-  as.vector(validLncTargetStat.df[,2])
tmp.idx <- match(LncPCGLink.df[,6], PCG.df$EntrezID)
LncPCGLink.df[,4] <- as.vector(PCG.df$seqnames[tmp.idx])
LncPCGLink.df[,5] <- PCG.df$start[tmp.idx]

################
AllLncPCGLink.df  <-  rbind(AllLncPCGLink.df, LncPCGLink.df)

LncPCGLink.df$SameChromo <- "No"
tmp.idx <-  which(LncPCGLink.df$Chr1 == LncPCGLink.df$Chr2)
LncPCGLink.df$SameChromo[tmp.idx] <- "Yes"

LncPCGLink.df$ChromoDist <-  abs(LncPCGLink.df$Pos1 - LncPCGLink.df$Pos2)
LncPCGLink.df$ChromoDist[which(LncPCGLink.df$SameChromo == "No")] <- NA


LncPCGLink.df$DistMark  <-  "Different Chromosome"
LncPCGLink.df$DistMark[which(LncPCGLink.df$ChromoDist <= 100000)] <- "< 100kb"
LncPCGLink.df$DistMark[intersect(which(LncPCGLink.df$ChromoDist > 100000), which(LncPCGLink.df$ChromoDist <= 10000000))] <- "100kb ~ 10Mb"
LncPCGLink.df$DistMark[intersect(which(LncPCGLink.df$ChromoDist > 10000000), which(LncPCGLink.df$ChromoDist <= 100000000))] <- "10Mb ~ 100Mb"
LncPCGLink.df$DistMark[which(LncPCGLink.df$ChromoDist > 100000000)] <- "> 100Mb"

colors.v <- colorRampPalette(brewer.pal(9, "BuGn"))(12)[c(7,5,4,3,2)]

# barplot(table(LncPCGLink.df$DistMark)[c(5,2,4,3,1)], las=2)
gaps.v <- c(table(LncPCGLink.df$DistMark)[c(5,2,4,3,1)][1] * 0.1, table(LncPCGLink.df$DistMark)[c(5,2,4,3,1)][1] * 0.9)
print(gaps.v )
yaxlab.v <- c(table(LncPCGLink.df$DistMark)[c(5,2,4,3,1)][1] * 0.2 *0.1, table(LncPCGLink.df$DistMark)[c(5,2,4,3,1)][1] * 0.2 *0.4, table(LncPCGLink.df$DistMark)[c(5,2,4,3,1)][1] *0.95,  table(LncPCGLink.df$DistMark)[c(5,2,4,3,1)][1] *0.99)
#yaxlab.v <- c(table(LncPCGLink.df$DistMark)[c(5,2,4,3,1)][1] *0.1, table(LncPCGLink.df$DistMark)[c(5,2,4,3,1)][1]  *0.4, table(LncPCGLink.df$DistMark)[c(5,2,4,3,1)][1] *0.6,  table(LncPCGLink.df$DistMark)[c(5,2,4,3,1)][1] *0.8)

print(as.numeric(round(yaxlab.v , -3)))
xlables.v <- c("Diff Chromo", "> 100Mb" ,  "10Mb ~ 100Mb" ,"100kb ~ 10Mb", "< 100kb" )
gap.barplot(table(LncPCGLink.df$DistMark)[c(5,2,4,3,1)], gap=gaps.v, ylim=c(0, table(LncPCGLink.df$DistMark)[c(5,2,4,3,1)][1]*0.2), main=cancers.v[c], xaxlab= xlables.v,
 col=colors.v, xlab="", ylab="", ytics= as.numeric(round(yaxlab.v , -3)) , xaxt="n", las=1)
minpos <- table(LncPCGLink.df$DistMark)[c(5,2,4,3,1)][1] * -0.015
axis(1,at =1:5, label=xlables.v , tick = TRUE, srt = 45, las=2, cex.axis =0.6,  font.axis = 2)
#text(y= -minpos, x =1:5, label=xlables.v , xpd=T, adj=0.5, srt = 45)

}
dev.off()



AllLncPCGLink.df <- AllLncPCGLink.df[-1,]
library(dplyr)

AllLncPCGLink.df <- distinct(AllLncPCGLink.df)

AllLncPCGLink.df$SameChromo <- "No"
tmp.idx <-  which(AllLncPCGLink.df$Chr1 == AllLncPCGLink.df$Chr2)
AllLncPCGLink.df$SameChromo[tmp.idx] <- "Yes"

AllLncPCGLink.df$ChromoDist <-  abs(AllLncPCGLink.df$Pos1 - AllLncPCGLink.df$Pos2)
AllLncPCGLink.df$ChromoDist[which(AllLncPCGLink.df$SameChromo == "No")] <- NA


AllLncPCGLink.df$DistMark  <-  NA
AllLncPCGLink.df$DistMark[which(AllLncPCGLink.df$ChromoDist <= 100000)] <- "< 100kb"
AllLncPCGLink.df$DistMark[intersect(which(AllLncPCGLink.df$ChromoDist > 100000), which(AllLncPCGLink.df$ChromoDist <= 10000000))] <- "100kb ~ 10Mb"
AllLncPCGLink.df$DistMark[intersect(which(AllLncPCGLink.df$ChromoDist > 10000000), which(AllLncPCGLink.df$ChromoDist <= 100000000))] <- "10Mb ~ 100Mb"
AllLncPCGLink.df$DistMark[which(AllLncPCGLink.df$ChromoDist > 100000000)] <- "> 100Mb"



slices1 <- round(as.numeric(table(AllLncPCGLink.df$SameChromo)) /nrow(AllLncPCGLink.df)*100, digits = 2)
names(slices1) <- c("Different\nChromosome\n", "Same\nChromosome\n");
lbls1 <- paste(names(slices1), " (", slices1, "%)", sep= "");

slices2 <- round(as.numeric(table(AllLncPCGLink.df$DistMark)) /sum(AllLncPCGLink.df$SameChromo == "Yes")*100, digits = 2)
names(slices2) <- names(table(AllLncPCGLink.df$DistMark))
lbls2 <- paste(names(slices2), "\n(", slices2, "%)", sep= "");
slices2 <- slices2[c(1,3,4,2)]
lbls2 <- lbls2[c(1,3,4,2)]

pdf("All_Interact_CisTrans_Pie.pdf", width=6.5,height=4)
par(family="ArialMT")
par(mar=c(0.28,0.2,0.1,0))
par(mfrow=c(1,2))
colors.v <- colorRampPalette(brewer.pal(9, "BuGn"))(12)[c(5,6)]
pie(slices1, labels="", col=colors.v,  init.angle = 10,  radius = 0.65, cex= 0.68, font = 2  )
text(0.4, 0.15, labels = lbls1[2], xpd = TRUE, cex = 0.8, font=2)
text(-0.2, -0.35, labels = lbls1[1], xpd = TRUE, cex = 0.8, font=2)
segments(0.65, 0.12, 1, 0.185, col= 'gray51', lty = 2)
segments(0.65, -0.13, 1, -0.20, col= 'gray51', lty = 2)

par(mar=c(0.28,0,0.1,0.5))
#text(0,1, labels = "Exp", xpd = TRUE, cex = 1.15, font=2)
colors.v <- colorRampPalette(brewer.pal(9, "BuGn"))(12)[c(2,3,4,5)]
pie(slices2, labels="", col=colors.v,  init.angle = 12,  radius = 0.45, cex= 0.68, font = 2  )
text(0.56, 0.13, labels = lbls2[1], xpd = TRUE, cex = 0.5, font=2)
text(0.28, 0.2, labels = lbls2[2], xpd = TRUE, cex = 0.5, font=2)
text(-0.2, -0.2, labels = lbls2[3], xpd = TRUE, cex = 0.5, font=2)
text(0.27, -0.12, labels = lbls2[4], xpd = TRUE, cex = 0.5, font=2)
dev.off()


