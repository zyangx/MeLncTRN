
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(topGO)
library(VennDiagram)


load("../CancerRdata/Exp_gene/BLCA_expGenelevel.Rdata")
allGenes.v <- rownames(PCGexpFPKM.m)

load("../CancerRdata/LncValExp0/KIRC_valFPKM.Rdata")
allLnc.v <- rownames(LncValFPKM.m)


cancers.v  <-  c("BLCA", "BRCA(ER+)", "BRCA(ER-)", "CESC", "CHOL", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC",
				"LUAD", "LUSC", "PAAD" , "PRAD", "READ", "THCA", "UCEC")

c <- 9
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
KIRCGenes.v <- colnames(InteractLG.m)
KIRCLnc.v <- rownames(InteractLG.m)

c <- 10
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
KIRPGenes.v <- colnames(InteractLG.m)
KIRPLnc.v <- rownames(InteractLG.m)

InterGenes.v <- intersect(KIRCGenes.v, KIRPGenes.v)
InterLnc.v <- intersect(KIRCLnc.v, KIRPLnc.v)

data.m <- matrix(NA, nrow=2, ncol=2)
data.m[1,1]  <- length(allGenes.v )
data.m[2,1] <- length(intersect(allGenes.v , KIRPGenes.v))
data.m[1,2] <- length(KIRCGenes.v)
data.m[2,2] <- length(InterGenes.v)
data.m
fisher.test(data.m)
fisher.test(data.m)$p.value

data.m <- matrix(NA, nrow=2, ncol=2)
data.m[1,1]  <- length(allLnc.v )
data.m[2,1] <- length(intersect(allLnc.v , KIRPLnc.v))
data.m[1,2] <- length(KIRCLnc.v)
data.m[2,2] <- length(InterLnc.v)
data.m
fisher.test(data.m)
fisher.test(data.m)$p.value

vennDiagram <-venn.diagram(x=list(KIRC=KIRCGenes.v, KIRP=KIRPGenes.v), filename =NULL, lwd = 3, 
       fill = c("cornflowerblue", "darkorchid1"),  col = "transparent", alpha = 0.4, label.col = "black",  cex = 1.5,
  fontfamily = "Arial",  fontface = "bold",  cat.col = c("black", "black"),  cat.cex = 1.2,
  cat.fontfamily = "Arial",  cat.fontface = "bold",  margin = 0.05,  cat.dist = c(0.03, 0.03),  cat.pos = c(-20, 20))

pdf("KIRC_KIRP_vennDiagram.pdf",width=4.5,height=2)
grid.draw(vennDiagram)
dev.off()

#######GO 
ego <- enrichGO(InterGenes.v, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.01, 
                 qvalueCutoff = 0.05, keyType = 'ENTREZID')
dim(ego)
#write.table(as.data.frame(ego), file="Fem_Genes_GOstat_BP_ALL.csv", sep="\t", row.names=F)

ego2 <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
dim(ego2)

ego3 <- setReadable(ego2, OrgDb = org.Hs.eg.db)
dim(ego3)
#head(summary(ego3))
enrichedResult.df <- as.data.frame(ego3)
enrichedResult.df$FDR <- -log10(enrichedResult.df$p.adjust)
#write.table(as.data.frame(ego3), file="InterGenes_GOstat_BP_Simplify.csv", sep="\t", row.names=F)
#row.idx <- c(1,2,3,5,6,13,22,26)
row.idx <- c(1,2,3,5,6,13,22,26)
enrichedResult.df <- enrichedResult.df[row.idx,]

goplot(ego3)
heatplot(ego3)
emapplot(ego3)

p <- barplot(ego3,showCategory=20,drop=T)
library(ggThemeAssist)
ggThemeAssistGadget(p)

dotplot(ego3,showCategory=50)
#enrichMap(ego3, vertex)
##cnetplot(ego3, categorySize="pvalue", circular=TRUE, colorEdge=TRUE, foldChange=geneList)


library(ggpubr)
p1 <- ggbarplot(enrichedResult.df, x="Description", y="FDR",  fill= "lightblue", color="white", sort.by.groups=FALSE, ylab="-log10(FDR)", xlab="Biological Process", ggtheme = theme_bw()) + theme(plot.title = element_text(hjust = 0.5),  axis.text.y=element_blank() )
p1 <- p1 + annotate("text",x=1:nrow(enrichedResult.df),y=0.2,label=rev(enrichedResult.df$Description), family="Arial", hjust=0)
p1 <- p1 + coord_flip()
p1 <- p1 + scale_x_discrete(limits=rev(enrichedResult.df$Description))
p1
ggsave("KIRC_KIRP_intersect_BP_enrich_clusterprofiler.pdf",device=cairo_pdf,width=4.5,height=4)

####### kegg
#ego <- enrichGO(InterGenes.v, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.01, qvalueCutoff = 0.05, keyType = 'ENTREZID')
ego <- enrichKEGG(gene = InterGenes.v, keyType = "kegg", organism = 'hsa', pvalueCutoff = 0.01,  pAdjustMethod= "BH", qvalueCutoff  = 0.05)

enrichedResult.df <- as.data.frame(ego)
enrichedResult.df$FDR <- -log10(enrichedResult.df$p.adjust)
#write.table(as.data.frame(ego3), file="InterGenes_GOstat_BP_Simplify.csv", sep="\t", row.names=F)
row.idx <- c(1,2,3,4,5,6,9,13,21,22,26)
enrichedResult.df <- enrichedResult.df[row.idx,]

p1 <- ggbarplot(enrichedResult.df, x="Description", y="FDR",  fill= "skyblue", color="white", sort.by.groups=FALSE, ylab="-log10(FDR)", xlab="Biological Process", ggtheme = theme_bw()) + theme(plot.title = element_text(hjust = 0.5),  axis.text.y=element_blank() )
p1 <- p1 + annotate("text",x=1:nrow(enrichedResult.df),y=0.2,label=rev(enrichedResult.df$Description), family="Arial", hjust=0)
p1 <- p1 + coord_flip()
p1 <- p1 + scale_x_discrete(limits=rev(enrichedResult.df$Description))
p1
ggsave("KIRC_KIRP_intersect_BP_enrich_2.pdf",device=cairo_pdf,width=4.5,height=4)

######## topGO
geneList <- factor(as.integer(allGenes.v %in% InterGenes.v))
names(geneList) <- allGenes.v 
str(geneList)

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, nodeSize = 10, annot = annFUN.org, mapping = "org.Hs.eg.db", ID="entrez")
result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
sig.tab <- GenTable(GOdata, Fis = result,  topNodes = 100)

gostat <- termStat(GOdata, names(score(result)))
showSigOfNodes(GOdata, score(result), firstSigNodes=5, useInfo="all")

sig.tab$p.adjust  <- p.adjust(sig.tab$Fis, method = "BH")*10^2
sig.tab$FDR <-  -log10(sig.tab$p.adjust)
row.idx <- c(4,5,10,11,22,26,33,48)
sig.tab  <- sig.tab[row.idx,]

p1 <- ggbarplot(sig.tab, x="Term", y="FDR",  fill= "lightblue", color="white", sort.by.groups=FALSE, ylab="-log10(FDR)", xlab="Biological Process", ggtheme = theme_bw()) + theme(plot.title = element_text(hjust = 0.5),  axis.text.y=element_blank() )
p1 <- p1 + annotate("text",x=1:nrow(sig.tab),y=0.2,label=rev(sig.tab$Term), family="Arial", hjust=0)
p1 <- p1 + coord_flip()
p1 <- p1 + scale_x_discrete(limits=rev(sig.tab$Term))
p1
ggsave("KIRC_KIRP_intersect_BP_enrich_TopGO.pdf",device=cairo_pdf,width=4.5,height=4)


#########
#########
load("../CancerRdata/Exp_gene/BLCA_expGenelevel.Rdata")
allGenes.v <- rownames(PCGexpFPKM.m)

load("../CancerRdata/LncValExp0/COAD_valFPKM.Rdata")
allLnc.v <- rownames(LncValFPKM.m)


c <- 6
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
COADGenes.v <- colnames(InteractLG.m)
COADLnc.v <- rownames(InteractLG.m)

c <- 16
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
READGenes.v <- colnames(InteractLG.m)
READLnc.v <- rownames(InteractLG.m)

InterGenes.v <- intersect(COADGenes.v, READGenes.v)
InterLnc.v <- intersect(COADLnc.v, READLnc.v)

data.m <- matrix(NA, nrow=2, ncol=2)
data.m[1,1]  <- length(allGenes.v )
data.m[2,1] <- length(intersect(allGenes.v , READGenes.v))
data.m[1,2] <- length(COADGenes.v)
data.m[2,2] <- length(InterGenes.v)
data.m
fisher.test(data.m)
fisher.test(data.m)$p.value

data.m <- matrix(NA, nrow=2, ncol=2)
data.m[1,1]  <- length(allLnc.v )
data.m[2,1] <- length(intersect(allLnc.v , READLnc.v))
data.m[1,2] <- length(COADLnc.v)
data.m[2,2] <- length(InterLnc.v)
data.m
fisher.test(data.m)
fisher.test(data.m)$p.value

vennDiagram <-venn.diagram(x=list(COAD=COADGenes.v, READ=READGenes.v), filename =NULL, lwd = 3, 
       fill = c("cornflowerblue", "darkorchid1"),  col = "transparent", alpha = 0.4, label.col = "black",  cex = 1.5,
  fontfamily = "Arial",  fontface = "bold",  cat.col = c("black", "black"),  cat.cex = 1.2,
  cat.fontfamily = "Arial",  cat.fontface = "bold",  margin = 0.05,  cat.dist = c(0.03, 0.03),  cat.pos = c(-20, 20))

pdf("COAD_READ_vennDiagram.pdf",width=4.5,height=2)
grid.draw(vennDiagram)
dev.off()


#######GO 
ego <- enrichGO(InterGenes.v, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.05, keyType = 'ENTREZID')
dim(ego)

#head(summary(ego3))
enrichedResult.df <- as.data.frame(ego)
enrichedResult.df$FDR <- -log10(enrichedResult.df$pvalue)
#write.table(as.data.frame(ego3), file="InterGenes_GOstat_BP_Simplify.csv", sep="\t", row.names=F)
row.idx <- c(8,10,12,13)
enrichedResult.df <- enrichedResult.df[row.idx,]


library(ggpubr)
p1 <- ggbarplot(enrichedResult.df, x="Description", y="FDR",  fill= "lightblue", color="white", sort.by.groups=FALSE, ylab="-log10(FDR)", xlab="Biological Process", ggtheme = theme_bw()) + theme(plot.title = element_text(hjust = 0.5),  axis.text.y=element_blank() )
p1 <- p1 + annotate("text",x=1:nrow(enrichedResult.df),y=0.2,label=rev(enrichedResult.df$Description), family="Arial", hjust=0)
p1 <- p1 + coord_flip()
p1 <- p1 + scale_x_discrete(limits=rev(enrichedResult.df$Description))
p1
ggsave("COAD_READ_intersect_BP_enrich_clusterprofiler.pdf",device=cairo_pdf,width=4.5,height=4)

####### kegg
#ego <- enrichGO(InterGenes.v, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.01, qvalueCutoff = 0.05, keyType = 'ENTREZID')
ego <- enrichKEGG(gene = InterGenes.v, keyType = "kegg", organism = 'hsa', pvalueCutoff = 0.01,  pAdjustMethod= "BH", qvalueCutoff  = 0.05)

enrichedResult.df <- as.data.frame(ego)
enrichedResult.df$FDR <- -log10(enrichedResult.df$p.adjust)
#write.table(as.data.frame(ego3), file="InterGenes_GOstat_BP_Simplify.csv", sep="\t", row.names=F)
row.idx <- c(1,2,3,4,5,6,9,13,21,22,26)
enrichedResult.df <- enrichedResult.df[row.idx,]

p1 <- ggbarplot(enrichedResult.df, x="Description", y="FDR",  fill= "skyblue", color="white", sort.by.groups=FALSE, ylab="-log10(FDR)", xlab="Biological Process", ggtheme = theme_bw()) + theme(plot.title = element_text(hjust = 0.5),  axis.text.y=element_blank() )
p1 <- p1 + annotate("text",x=1:nrow(enrichedResult.df),y=0.2,label=rev(enrichedResult.df$Description), family="Arial", hjust=0)
p1 <- p1 + coord_flip()
p1 <- p1 + scale_x_discrete(limits=rev(enrichedResult.df$Description))
p1

######## topGO
geneList <- factor(as.integer(allGenes.v %in% InterGenes.v))
names(geneList) <- allGenes.v 
str(geneList)

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, nodeSize = 10, annot = annFUN.org, mapping = "org.Hs.eg.db", ID="entrez")
result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
sig.tab <- GenTable(GOdata, Fis = result,  topNodes = 100)


sig.tab$p.adjust  <- p.adjust(sig.tab$Fis, method = "BH")*10^2
sig.tab$FDR <-  -log10(as.numeric(sig.tab$Fis))
#row.idx <- c(6,11,26,54,83,88)
row.idx <- c(20,18,41,65,70,82)
sig.tab  <- sig.tab[row.idx,]

p1 <- ggbarplot(sig.tab, x="Term", y="FDR",  fill= "lightblue", color="white", sort.by.groups=FALSE, ylab="-log10(FDR)", xlab="Biological Process", ggtheme = theme_bw()) + theme(plot.title = element_text(hjust = 0.5),  axis.text.y=element_blank() )
p1 <- p1 + annotate("text",x=1:nrow(sig.tab),y=0.2,label=rev(sig.tab$Term), family="Arial", hjust=0)
p1 <- p1 + coord_flip()
p1 <- p1 + scale_x_discrete(limits=rev(sig.tab$Term))
p1
ggsave("COAD_READ_intersect_BP_enrich_TopGO.pdf",device=cairo_pdf,width=4.5,height=4)


#########
#########

load("../CancerRdata/Exp_gene/BLCA_expGenelevel.Rdata")
allGenes.v <- rownames(PCGexpFPKM.m)

load("../CancerRdata/LncValExp0/LUSC_valFPKM.Rdata")
allLnc.v <- rownames(LncValFPKM.m)


c <- 12
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
LUADGenes.v <- colnames(InteractLG.m)
LUADLnc.v <- rownames(InteractLG.m)

c <- 13
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
LUSCGenes.v <- colnames(InteractLG.m)
LUSCLnc.v <- rownames(InteractLG.m)

InterGenes.v <- intersect(LUADGenes.v, LUSCGenes.v)
InterLnc.v <- intersect(LUADLnc.v, LUSCLnc.v)

data.m <- matrix(NA, nrow=2, ncol=2)
data.m[1,1]  <- length(allGenes.v )
data.m[2,1] <- length(intersect(allGenes.v , LUSCGenes.v))
data.m[1,2] <- length(LUADGenes.v)
data.m[2,2] <- length(InterGenes.v)
data.m
fisher.test(data.m)
fisher.test(data.m)$p.value

data.m <- matrix(NA, nrow=2, ncol=2)
data.m[1,1]  <- length(allLnc.v )
data.m[2,1] <- length(intersect(allLnc.v , LUSCLnc.v))
data.m[1,2] <- length(LUADLnc.v)
data.m[2,2] <- length(InterLnc.v)
data.m
fisher.test(data.m)
fisher.test(data.m)$p.value


vennDiagram <-venn.diagram(x=list(LUAD=LUADGenes.v, LUSC=LUSCGenes.v), filename =NULL, lwd = 3, 
       fill = c("cornflowerblue", "darkorchid1"),  col = "transparent", alpha = 0.4, label.col = "black",  cex = 1.5,
  fontfamily = "Arial",  fontface = "bold",  cat.col = c("black", "black"),  cat.cex = 1.2,
  cat.fontfamily = "Arial",  cat.fontface = "bold",  margin = 0.05,  cat.dist = c(0.03, 0.03),  cat.pos = c(-20, 20))

pdf("LUAD_LUSC_vennDiagram.pdf",width=4.5,height=2)
grid.draw(vennDiagram)
dev.off()

#######GO 
ego <- enrichGO(InterGenes.v, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.05, keyType = 'ENTREZID')
dim(ego)

#head(summary(ego3))
enrichedResult.df <- as.data.frame(ego)
enrichedResult.df$FDR <- -log10(enrichedResult.df$pvalue)
#write.table(as.data.frame(ego3), file="InterGenes_GOstat_BP_Simplify.csv", sep="\t", row.names=F)
row.idx <- c(4,6,30,38,45,53)
enrichedResult.df <- enrichedResult.df[row.idx,]


library(ggpubr)
p1 <- ggbarplot(enrichedResult.df, x="Description", y="FDR",  fill= "lightpink", color="white", sort.by.groups=FALSE, ylab="-log10(FDR)", xlab="Biological Process", ggtheme = theme_bw()) + theme(plot.title = element_text(hjust = 0.5),  axis.text.y=element_blank() )
p1 <- p1 + annotate("text",x=1:nrow(enrichedResult.df),y=0.2,label=rev(enrichedResult.df$Description), family="Arial", hjust=0)
p1 <- p1 + coord_flip()
p1 <- p1 + scale_x_discrete(limits=rev(enrichedResult.df$Description))
p1
ggsave("LUAD_LUSC_intersect_BP_enrich_clusterprofiler.pdf",device=cairo_pdf,width=4.5,height=4)

####### kegg
#ego <- enrichGO(InterGenes.v, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.01, qvalueCutoff = 0.05, keyType = 'ENTREZID')
ego <- enrichKEGG(gene = InterGenes.v, keyType = "kegg", organism = 'hsa', pvalueCutoff = 0.01,  pAdjustMethod= "BH", qvalueCutoff  = 0.05)

enrichedResult.df <- as.data.frame(ego)
enrichedResult.df$FDR <- -log10(enrichedResult.df$p.adjust)
#write.table(as.data.frame(ego3), file="InterGenes_GOstat_BP_Simplify.csv", sep="\t", row.names=F)
row.idx <- c(1,2,3,4,5,6,9,13,21,22,26)
enrichedResult.df <- enrichedResult.df[row.idx,]

p1 <- ggbarplot(enrichedResult.df, x="Description", y="FDR",  fill= "skyblue", color="white", sort.by.groups=FALSE, ylab="-log10(FDR)", xlab="Biological Process", ggtheme = theme_bw()) + theme(plot.title = element_text(hjust = 0.5),  axis.text.y=element_blank() )
p1 <- p1 + annotate("text",x=1:nrow(enrichedResult.df),y=0.2,label=rev(enrichedResult.df$Description), family="Arial", hjust=0)
p1 <- p1 + coord_flip()
p1 <- p1 + scale_x_discrete(limits=rev(enrichedResult.df$Description))
p1

######## topGO
geneList <- factor(as.integer(allGenes.v %in% InterGenes.v))
names(geneList) <- allGenes.v 
str(geneList)

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, nodeSize = 10, annot = annFUN.org, mapping = "org.Hs.eg.db", ID="entrez")
result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
sig.tab <- GenTable(GOdata, Fis = result,  topNodes = 100)


sig.tab$p.adjust  <- p.adjust(sig.tab$Fis, method = "BH")*10^2
sig.tab$FDR <-  -log10(as.numeric(sig.tab$Fis))
row.idx <- c(6,11,21,35,38,47,48)
sig.tab  <- sig.tab[row.idx,]

p1 <- ggbarplot(sig.tab, x="Term", y="FDR",  fill= "lightblue", color="white", sort.by.groups=FALSE, ylab="-log10(FDR)", xlab="Biological Process", ggtheme = theme_bw()) + theme(plot.title = element_text(hjust = 0.5),  axis.text.y=element_blank() )
p1 <- p1 + annotate("text",x=1:nrow(sig.tab),y=0.2,label=rev(sig.tab$Term), family="Arial", hjust=0)
p1 <- p1 + coord_flip()
p1 <- p1 + scale_x_discrete(limits=rev(sig.tab$Term))
p1
ggsave("LUAD_LUSC_intersect_BP_enrich_TopGO.pdf",device=cairo_pdf,width=4.5,height=4)




#########
#########

load("../CancerRdata/Exp_gene/BLCA_expGenelevel.Rdata")
allGenes.v <- rownames(PCGexpFPKM.m)

load("../CancerRdata/LncValExp0/LIHC_valFPKM.Rdata")
allLnc.v <- rownames(LncValFPKM.m)


c <- 11
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
LIHCGenes.v <- colnames(InteractLG.m)
LIHCLnc.v <- rownames(InteractLG.m)

c <- 5
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
CHOLGenes.v <- colnames(InteractLG.m)
CHOLLnc.v <- rownames(InteractLG.m)

InterGenes.v <- intersect(LIHCGenes.v, CHOLGenes.v)
InterLnc.v <- intersect(LIHCLnc.v, CHOLLnc.v)

data.m <- matrix(NA, nrow=2, ncol=2)
data.m[1,1]  <- length(allGenes.v )
data.m[2,1] <- length(intersect(allGenes.v , LIHCGenes.v))
data.m[1,2] <- length(CHOLGenes.v)
data.m[2,2] <- length(InterGenes.v)
data.m
fisher.test(data.m)
fisher.test(data.m)$p.value

data.m <- matrix(NA, nrow=2, ncol=2)
data.m[1,1]  <- length(allLnc.v )
data.m[2,1] <- length(intersect(allLnc.v , LIHCLnc.v))
data.m[1,2] <- length(CHOLLnc.v)
data.m[2,2] <- length(InterLnc.v)
data.m
fisher.test(data.m)
fisher.test(data.m)$p.value


vennDiagram <-venn.diagram(x=list(LIHC=LIHCGenes.v, CHOL=CHOLGenes.v), filename =NULL, lwd = 3, 
       fill = c("cornflowerblue", "darkorchid1"),  col = "transparent", alpha = 0.4, label.col = "black",  cex = 1.5,
  fontfamily = "Arial",  fontface = "bold",  cat.col = c("black", "black"),  cat.cex = 1.2,
  cat.fontfamily = "Arial",  cat.fontface = "bold",  margin = 0.05,  cat.dist = c(0.03, 0.03),  cat.pos = c(-20, 20))

pdf("LIHC_CHOL_vennDiagram.pdf",width=4.5,height=2)
grid.draw(vennDiagram)
dev.off()



#######GO 
ego <- enrichGO(InterGenes.v, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.05, keyType = 'ENTREZID')
dim(ego)

#head(summary(ego3))
enrichedResult.df <- as.data.frame(ego)
enrichedResult.df$FDR <- -log10(enrichedResult.df$pvalue)
#write.table(as.data.frame(ego3), file="InterGenes_GOstat_BP_Simplify.csv", sep="\t", row.names=F)
row.idx <- c(1,3,4,5,13,16,34)
enrichedResult.df <- enrichedResult.df[row.idx,]


library(ggpubr)
p1 <- ggbarplot(enrichedResult.df, x="Description", y="FDR",  fill= "lightblue", color="white", sort.by.groups=FALSE, ylab="-log10(FDR)", xlab="Biological Process", ggtheme = theme_bw()) + theme(plot.title = element_text(hjust = 0.5),  axis.text.y=element_blank() )
p1 <- p1 + annotate("text",x=1:nrow(enrichedResult.df),y=0.2,label=rev(enrichedResult.df$Description), family="Arial", hjust=0)
p1 <- p1 + coord_flip()
p1 <- p1 + scale_x_discrete(limits=rev(enrichedResult.df$Description))
p1
ggsave("LIHC_CHOL_intersect_BP_enrich_clusterprofiler.pdf",device=cairo_pdf,width=4.5,height=4)



######## topGO
geneList <- factor(as.integer(allGenes.v %in% InterGenes.v))
names(geneList) <- allGenes.v 
str(geneList)

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, nodeSize = 10, annot = annFUN.org, mapping = "org.Hs.eg.db", ID="entrez")
result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
sig.tab <- GenTable(GOdata, Fis = result,  topNodes = 100)


sig.tab$p.adjust  <- p.adjust(sig.tab$Fis, method = "BH")*10^2
sig.tab$FDR <-  -log10(as.numeric(sig.tab$Fis))
row.idx <- c(4,6,8,10,18,27,37,55 )
sig.tab  <- sig.tab[row.idx,]

p1 <- ggbarplot(sig.tab, x="Term", y="FDR",  fill= "lightblue", color="white", sort.by.groups=FALSE, ylab="-log10(FDR)", xlab="Biological Process", ggtheme = theme_bw()) + theme(plot.title = element_text(hjust = 0.5),  axis.text.y=element_blank() )
p1 <- p1 + annotate("text",x=1:nrow(sig.tab),y=0.2,label=rev(sig.tab$Term), family="Arial", hjust=0)
p1 <- p1 + coord_flip()
p1 <- p1 + scale_x_discrete(limits=rev(sig.tab$Term))
p1
ggsave("LIHC_CHOL_intersect_BP_enrich_TopGO.pdf",device=cairo_pdf,width=4.5,height=4)






