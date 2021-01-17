
cancers.v  <-  c("BLCA", "BRCA(ER+)", "BRCA(ER-)", "CESC", "CHOL", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC",
				"LUAD", "LUSC", "PAAD" , "PRAD", "READ", "THCA", "UCEC")

SampleStat.m <- matrix(NA, nrow=18, ncol=4)
rownames(SampleStat.m) <- cancers.v
colnames(SampleStat.m) <- c("GeneExp", "DNAm", "CNV", "Intersect")


ValidGeneStat.m <- matrix(NA, nrow=length(cancers.v), ncol=1)
rownames(ValidGeneStat.m)  <-  cancers.v
colnames(ValidGeneStat.m)[1]  <-  c("Num")


for(c in 1:length(cancers.v)) {
print(cancers.v[c])

#source("Dolimma.R")

load(paste("./Exp_gene/", cancers.v[c], "_expGenelevel.Rdata", sep=""))
#PhenoTypesExp.lv <- PhenoTypes.lv
#DEG.m <- Dolimma(PCGexpFPKM.m, as.vector(PhenoTypesExp.lv$Cancer),  "1" , "0")
SampleStat.m[c,1] <- ncol(PCGexpFPKM.m)

load(paste("./DNAm_gene/", cancers.v[c], "_DNAmGenelevel.Rdata", sep=""))
#PhenoTypesDNAm.lv <- PhenoTypes.lv
#DMG.m <- Dolimma(avbeta.m, as.vector(PhenoTypesDNAm.lv$Cancer),  "1" , "0")
SampleStat.m[c,2] <- ncol(avbeta.m)

#load(paste("./LncValExp0/", cancers.v[c], "_valFPKM.Rdata", sep=""))
#PhenoTypesLnc.lv <- PhenoTypes.lv
#DEL.m <- Dolimma(LncValFPKM.m, as.vector(PhenoTypesLnc.lv$Cancer),  "1" 

load(paste("./CNV_gene/", cancers.v[c], "_cnvGenelevel.Rdata", sep=""))
colnames(CNVscore.m)  <-  paste0(colnames(CNVscore.m), "A") 
SampleStat.m[c,3] <- ncol(CNVscore.m)

interGenes.v <- intersect(intersect(rownames(PCGexpFPKM.m), rownames(avbeta.m)), rownames(CNVscore.m))
interSamples.v <- intersect(intersect(colnames(PCGexpFPKM.m), colnames(avbeta.m)), colnames(CNVscore.m)  )
SampleStat.m[c,4] <- length(interSamples.v)


load(paste("./MultiOmicsData/", cancers.v[c], "_MultiOmicsProcess.Rdata", sep=""))
ValidGeneStat.m[c,1] <- nrow(statLMG.lm[[1]])

}


write.table(SampleStat.m, file="SampleStat_Exp_DNAm_CNV.txt", quote=F, sep="\t", col.names=NA)


library(RColorBrewer)
library(extrafont)
font_import()
library(ggsci)
library(scales)

slices1 <- as.numeric(SampleStat.m[,1])
names(slices1) <- cancers.v;
lbls1 <- paste(names(slices1), " (", slices1, ")", sep= "");

slices2 <- as.numeric(SampleStat.m[,2])
names(slices2) <- cancers.v;
lbls2 <- paste(names(slices2), " (", slices2, ")", sep= "");

slices3 <- as.numeric(SampleStat.m[,3])
names(slices3) <- cancers.v;
lbls3 <- paste(names(slices3), " (", slices3, ")", sep= "");

slices4 <- as.numeric(SampleStat.m[,4])
names(slices4) <- cancers.v;
lbls4 <- paste(names(slices4), " (", slices4, ")", sep= "");





pdf("SampleStat_Exp_DNAm_CNV.pdf", width=8, height=8)
par(family="ArialMT")
par(mar=c(0.28,0.2,0.1,0.52))
par(mfrow=c(2,2))
pie(slices1, labels=lbls1, col=pal_igv("default")(length(lbls1)), clockwise=TRUE, radius = 0.65, cex= 0.68, font = 2  )
text(-1,0.85, labels = "A", xpd = TRUE, cex = 1.2, font=2)
#text(0,1, labels = "Exp", xpd = TRUE, cex = 1.15, font=2)

pie(slices2, labels=lbls2, col=pal_igv("default")(length(lbls2)), clockwise=TRUE, radius = 0.65, cex= 0.68, font = 2  )
text(-1,0.85, labels = "B", xpd = TRUE, cex = 1.2, font=2)
#text(0,1, labels = "Exp", xpd = TRUE, cex = 1.15, font=2)

pie(slices3, labels=lbls3, col=pal_igv("default")(length(lbls3)), clockwise=TRUE, radius = 0.65, cex= 0.68, font = 2  )
text(-1,0.85, labels = "C", xpd = TRUE, cex = 1.2, font=2)
#text(0,1, labels = "Exp", xpd = TRUE, cex = 1.15, font=2)

pie(slices4, labels=lbls4, col=pal_igv("default")(length(lbls4)), clockwise=TRUE, radius = 0.65, cex= 0.68, font = 2  )
text(-1,0.85, labels = "D", xpd = TRUE, cex = 1.2, font=2)
#text(0,1, labels = "Exp", xpd = TRUE, cex = 1.15, font=2)

dev.off()



mycolors.v <- pal_d3("category20")(18)


library(ggplot2)
ValidGeneStat.df <- as.data.frame(ValidGeneStat.m)
ValidGeneStat.df$Cancer <- rownames(ValidGeneStat.df )
p <- ggplot(ValidGeneStat.df, aes(x = Cancer, y = Num, fill = Cancer)) + geom_bar(stat ="identity", fill='#1B9E77') + theme_bw() + theme(axis.text.x = element_text(angle=30, hjust=1,vjust=1), text=element_text(family="Arial"), legend.position='none')
#p <- p + scale_fill_manual(values =mycolors.v )
p 

ggsave("DNAmDrivenGenesStat.pdf", device=cairo_pdf, width=10,height=4)

