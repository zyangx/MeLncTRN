load("probe450kfemanno.rda");

cancers.v  <-  c("BLCA", "BRCA(ER+)", "BRCA(ER-)", "CESC", "CHOL", "COAD", "ESCA", "HNSC",
				 "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PAAD"  , "PRAD", "READ",  "THCA", "UCEC")

#load("../DNAm_450k/BLCA_bmiqbeta.Rdata")
for(i in 1:length(cancers.v)){
load(paste("../DNAm_450k/", cancers.v[i], "_bmiqbeta.Rdata", sep=""))
print(cancers.v[i])
rm(avbeta.m)
dnaM.m <- bmiqbeta.m

extractFn <- function(tmp.v, ext.idx) {
	 return(tmp.v[ext.idx])
}

map.idx <- match(rownames(dnaM.m), probe450kfemanno$probeID);
probeInfo.lv <- lapply(probe450kfemanno, extractFn, map.idx)
beta.lm <- list()


for (g in 1:6) {
    group.idx <- which(probeInfo.lv[[3]] == g)
    tmp.m <- dnaM.m[group.idx, ]
    rownames(tmp.m) <- probeInfo.lv$eid[group.idx]; ## genes with region(group) g
    sel.idx <- which(is.na(rownames(tmp.m)) == FALSE);
    tmp.m <- tmp.m[sel.idx,];
    nL <- length(factor(rownames(tmp.m)));
    nspg.v <- summary(factor(rownames(tmp.m)),maxsum=nL);
    beta.lm[[g]] <- rowsum(tmp.m,group=rownames(tmp.m))/nspg.v;  ### aggregate??
    print(paste("Done for regional gene group ", g, sep = ""))
}

unqEID.v <- unique(c(rownames(beta.lm[[2]]), rownames(beta.lm[[4]]), rownames(beta.lm[[1]])))

avbeta.m <- matrix(nrow = length(unqEID.v), ncol = ncol(dnaM.m))
colnames(avbeta.m) <- colnames(dnaM.m)
rownames(avbeta.m) <- unqEID.v
for (gr in c(1, 4, 2)) { ### order is important TSS1500, 1stExon, TSS200
	avbeta.m[match(rownames(beta.lm[[gr]]), rownames(avbeta.m)), ] <- beta.lm[[gr]]
}

save(avbeta.m, PhenoTypes.lv, file=paste(cancers.v[i], "_DNAmGenelevel.Rdata", sep="")	)
}


