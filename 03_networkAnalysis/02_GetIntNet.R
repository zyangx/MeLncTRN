options(stringsAsFactors=FALSE)
cancers.v  <-  c("BLCA", "BRCA(ER+)", "BRCA(ER-)", "CESC", "CHOL", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC",
				"LUAD", "LUSC", "PAAD" , "PRAD", "READ", "THCA", "UCEC")

for(c in 1:length(cancers.v)) {
#rm(list=ls())
#load("BLCA_MultiOmicsProcess.Rdata")
print(cancers.v[c])
load(paste("./MultiOmicsData/", cancers.v[c], "_MultiOmicsProcess.Rdata", sep="") )

#interaction.v <- c()
#tmp.idx <- which(abs(DEL.m[,1]) > 0.7)
#validLnc.v <-  rownames(DEL.m)[tmp.idx]


tmp.idx <- which(abs(DEG.m[,1]) > 0.0) #### exp fold change screen
DEMGenes.v <- intersect(rownames(DMG.m), rownames(DEG.m)[tmp.idx])
DNAmDrivenGenes.v <- intersect(DEMGenes.v,  rownames(statEMC.m))

InteractLG.m <- matrix(0, nrow=length(statLMG.lm), ncol=nrow(statEMC.m))
rownames(InteractLG.m) <- names(statLMG.lm)
colnames(InteractLG.m) <- rownames(statEMC.m)


for (i in 1:length(statLMG.lm)){
        print(i)

	stat.m <- statLMG.lm[[i]]
	pnTarget.idx <- intersect( intersect(which(stat.m[,1] > 0.3), which(stat.m[,4] < 0.05)), intersect(which(stat.m[,5] < -0.3), which(stat.m[,8] < 0.05)) )
	npTarget.idx <- intersect( intersect(which(stat.m[,1] < -0.3), which(stat.m[,4] < 0.05)), intersect(which(stat.m[,5] > 0.3), which(stat.m[,8] < 0.05)) )
	Lnctargets.idx <- sort(c(pnTarget.idx, npTarget.idx ))
    Lnctargets.v <- rownames(stat.m)[Lnctargets.idx]
	validLncTargets.v <- intersect(Lnctargets.v , DNAmDrivenGenes.v )


	tmp.idx <- match(validLncTargets.v, colnames(InteractLG.m))
	InteractLG.m[i, tmp.idx] <- 1
	#interaction.v <- c(interaction.v, length(validLncTargets.v))
}

row.idx <- which(rowSums(InteractLG.m) > 0)
col.idx <- which(colSums(InteractLG.m) > 0)

InteractLG.m <- InteractLG.m[row.idx, col.idx]
print(dim(InteractLG.m))

save(InteractLG.m, file=paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )

}



