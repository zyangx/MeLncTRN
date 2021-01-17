
options(stringsAsFactors=FALSE)
library("naturalsort")
library(org.Hs.eg.db)
library(stringr) 

cancers.v  <-  c("BLCA", "BRCA(ER+)", "BRCA(ER-)", "CESC", "CHOL", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC",
				"LUAD", "LUSC", "PAAD" , "PRAD", "READ", "THCA", "UCEC")


load("./Circos/gencode.Rdata")

EnsemblIDs.v <- str_split(PCG.df$gene_id,'[.]',simplify = T)[,1]
EntrezIDs.v <- mapIds(org.Hs.eg.db, keys=EnsemblIDs.v, column= 'ENTREZID', keytype="ENSEMBL")
PCG.df$EntrezID <- EntrezIDs.v 



AllLncPCGLinkClass.df  <- matrix(NA, nrow=18, ncol=4)
colnames(AllLncPCGLinkClass.df) <- c("Enhance_Activate", "Enhance_Inhibit", "Invert_Activate", "Invert_Inhibit")
rownames(AllLncPCGLinkClass.df) <- cancers.v

for(c in 1:length(cancers.v)) {
#rm(list=ls())
print(cancers.v[c])

load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
load(paste("./InteractNetworkStat1/", cancers.v[c], "_InteractNetStat.Rdata", sep="") )
load(paste("./MultiOmicsData/", cancers.v[c], "_MultiOmicsProcess.Rdata", sep=""))

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
LncPCGLink.df$RegulateType <- "Invert"
LncPCGLink.df[which(validLncTargetStat.df[,7] > 0), 7] <- "Enhance"

tmp.idx <- match(LncPCGLink.df[,6], rownames(DEG.m))
LncPCGLink.df$RegulateClass <- "Inhibit"
LncPCGLink.df[which(DEG.m[tmp.idx,1] > 0), 8] <- "Activate"

AllLncPCGLinkClass.df[c,] <-  table(paste(LncPCGLink.df[,7],  LncPCGLink.df[,8], sep="_"))

save(LncPCGLink.df, file=paste("./InteractNetworkStat1/", cancers.v[c], "_InteractNetStatClass.Rdata", sep="") )

}

tmp.df <- AllLncPCGLinkClass.df
tmp.df <- t( tmp.df /rowSums(tmp.df) )

Name <- rownames(tmp.df)
tmp.df <- data.frame(Name,tmp.df)
colnames(tmp.df) <- c("Name", cancers.v)
rownames(tmp.df) <- LETTERS[1:4]
tmp.df

library(ggradar)
ggradar(tmp.df,legend.position="right", group.point.size=0.1, group.line.width=0.8)

pdf("LncPCG_Interact_Class_Radar.pdf", width=8, height=5)
ggradar(tmp.df,legend.position="right", group.point.size=0.1, group.line.width=0.8)
dev.off()



