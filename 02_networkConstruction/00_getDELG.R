options(stringsAsFactors=FALSE)

cancers.v  <-  c("BLCA", "BRCA(ER+)", "BRCA(ER-)", "CESC", "CHOL", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC",
				"LUAD", "LUSC", "PAAD" , "PRAD", "READ", "THCA", "UCEC")

source("Dolimma.R")

LncDEInfo.df  <- as.data.frame(matrix(NA, nrow=1, ncol=9))
colnames(LncDEInfo.df) <-  c("CancerType", "LncID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val","fiveNumControl", "fiveNumCase")

for(c in 1:length(cancers.v)) {
print(cancers.v[c])
load(paste("./LncValExp0/", cancers.v[c], "_valFPKM.Rdata", sep=""))
DEL.m <- Dolimma(LncValFPKM.m, as.vector(PhenoTypes.lv$Cancer),  "1" , "0")

fiveNum.m <- matrix(NA, nrow=nrow(LncValFPKM.m), ncol=2)
rownames(fiveNum.m) <- rownames(LncValFPKM.m)
colnames(fiveNum.m) <- c("fiveNumControl", "fiveNumCase")
for(i in 1:nrow(fiveNum.m)){
fiveNum.m[i,1] <-  paste(format(quantile(LncValFPKM.m[i, which(PhenoTypes.lv$Cancer == 0)]),digits=3), collapse =", ")
fiveNum.m[i,2] <-  paste(format(quantile(LncValFPKM.m[i, which(PhenoTypes.lv$Cancer == 1)]),digits=3), collapse =", ")
}

tmp.idx <- match(rownames(DEL.m), rownames(fiveNum.m))
diffExpInfo.m <- cbind( cbind( cbind(cancers.v[c], rownames(DEL.m)), DEL.m[,-6]), fiveNum.m[tmp.idx,])
colnames(diffExpInfo.m)[1:2] <- c("CancerType", "LncID")

LncDEInfo.df  <-  rbind(LncDEInfo.df, diffExpInfo.m)

}

LncDEInfo.df <- LncDEInfo.df[-1,]
write.table(LncDEInfo.df , file="LncDEInfo.txt", sep="\t", quote=F, row.names=F)


##########

PCGDEInfo.df  <- as.data.frame(matrix(NA, nrow=1, ncol=9))
colnames(PCGDEInfo.df) <-  c("CancerType", "GeneID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val","fiveNumControl", "fiveNumCase")

for(c in 1:length(cancers.v)) {
print(cancers.v[c])
load(paste("./Exp_gene/", cancers.v[c], "_expGenelevel.Rdata", sep=""))
DEG.m <- Dolimma(PCGexpFPKM.m, as.vector(PhenoTypes.lv$Cancer),  "1" , "0")

fiveNum.m <- matrix(NA, nrow=nrow(PCGexpFPKM.m), ncol=2)
rownames(fiveNum.m) <- rownames(PCGexpFPKM.m)
colnames(fiveNum.m) <- c("fiveNumControl", "fiveNumCase")
for(i in 1:nrow(fiveNum.m)){
fiveNum.m[i,1] <-  paste(format(quantile(PCGexpFPKM.m[i, which(PhenoTypes.lv$Cancer == 0)]),digits=3), collapse =", ")
fiveNum.m[i,2] <-  paste(format(quantile(PCGexpFPKM.m[i, which(PhenoTypes.lv$Cancer == 1)]),digits=3), collapse =", ")
}

tmp.idx <- match(rownames(DEG.m), rownames(fiveNum.m))
diffExpInfo.m <- cbind( cbind( cbind(cancers.v[c], rownames(DEG.m)), DEG.m[,-6]), fiveNum.m[tmp.idx,])
colnames(diffExpInfo.m)[1:2] <- c("CancerType", "GeneID")

PCGDEInfo.df  <-  rbind(PCGDEInfo.df, diffExpInfo.m)

}

PCGDEInfo.df <- PCGDEInfo.df[-1,]
write.table(PCGDEInfo.df , file="GeneDEInfo.txt", sep="\t", quote=F, row.names=F)


#########

PCGDMInfo.df  <- as.data.frame(matrix(NA, nrow=1, ncol=9))
colnames(PCGDMInfo.df) <-  c("CancerType", "GeneID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val","fiveNumControl", "fiveNumCase")

for(c in 1:length(cancers.v)) {
print(cancers.v[c])
load(paste("./DNAm_gene/", cancers.v[c], "_DNAmGenelevel.Rdata", sep=""))
DMG.m <- Dolimma(avbeta.m, as.vector(PhenoTypes.lv$Cancer),  "1" , "0")

fiveNum.m <- matrix(NA, nrow=nrow(avbeta.m), ncol=2)
rownames(fiveNum.m) <- rownames(avbeta.m)
colnames(fiveNum.m) <- c("fiveNumControl", "fiveNumCase")
for(i in 1:nrow(fiveNum.m)){
fiveNum.m[i,1] <-  paste(format(quantile(avbeta.m[i, which(PhenoTypes.lv$Cancer == 0)]),digits=3), collapse =", ")
fiveNum.m[i,2] <-  paste(format(quantile(avbeta.m[i, which(PhenoTypes.lv$Cancer == 1)]),digits=3), collapse =", ")
}

tmp.idx <- match(rownames(DMG.m), rownames(fiveNum.m))
diffDNAmInfo.m <- cbind( cbind( cbind(cancers.v[c], rownames(DMG.m)), DMG.m[,-6]), fiveNum.m[tmp.idx,])
colnames(diffDNAmInfo.m)[1:2] <- c("CancerType", "GeneID")

PCGDMInfo.df  <-  rbind(PCGDMInfo.df, diffDNAmInfo.m)

}

PCGDMInfo.df <- PCGDMInfo.df[-1,]
write.table(PCGDMInfo.df , file="GeneDMInfo.txt", sep="\t", quote=F, row.names=F)



