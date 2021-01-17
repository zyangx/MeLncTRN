options(stringsAsFactors=FALSE)

cancers.v  <-  c("BLCA", "BRCA(ER+)", "BRCA(ER-)", "CESC", "CHOL", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC",
				"LUAD", "LUSC", "PAAD" , "PRAD", "READ", "THCA", "UCEC")

#mycolors.v <- pal_d3("category20")(18)

for(c in 1:length(cancers.v)) {

print(cancers.v[c])
load(paste("./MultiOmicsData/", cancers.v[c], "_MultiOmicsProcess.Rdata", sep="") )


tmp.idx <- which(abs(DEG.m[,1]) > 0.0) #### exp fold change screen
DEMGenes.v <- intersect(rownames(DMG.m), rownames(DEG.m)[tmp.idx])
DNAmDrivenGenes.v <- intersect(DEMGenes.v,  rownames(statEMC.m))


validLncTargetStat.df <- data.frame("corM" = 0, "tM"  = 0,  "PM" = 0 ,  "FDRM" = 0, "corG" = 0, "tG" = 0 ,  "PG" = 0,   "FDRG" = 0, "LncID" = 0, stringsAsFactors=F)
validLncTargetStat.df <- validLncTargetStat.df[-1,]


for (i in 1:length(statLMG.lm)){
        print(i)
	
	stat.m <- statLMG.lm[[i]]
	pnTarget.idx <- intersect( intersect(which(stat.m[,1] > 0.3), which(stat.m[,4] < 0.05)), intersect(which(stat.m[,5] < -0.3), which(stat.m[,8] < 0.05)) )
	npTarget.idx <- intersect( intersect(which(stat.m[,1] < -0.3), which(stat.m[,4] < 0.05)), intersect(which(stat.m[,5] > 0.3), which(stat.m[,8] < 0.05)) )
	Lnctargets.idx <- sort(c(pnTarget.idx, npTarget.idx ))
        Lnctargets.v <- rownames(stat.m)[Lnctargets.idx]
	validLncTargets.v <- intersect(Lnctargets.v , DNAmDrivenGenes.v )

	
	if(length(validLncTargets.v) > 0){
		validLnctargets.idx <- match(validLncTargets.v, rownames(stat.m))
		if(length(validLncTargets.v) > 1){
		    validLncStat.df <- as.data.frame(stat.m[validLnctargets.idx,], stringsAsFactors=F)
		    validLncStat.df$LncID <- as.character(names(statLMG.lm)[i])
		    validLncStat.df$GeneID <-  as.character(rownames(validLncStat.df))
		    validLncTargetStat.df <- rbind(validLncTargetStat.df, validLncStat.df )
		}else{
		    validLncStat.df <- as.data.frame(stat.m[validLnctargets.idx,], stringsAsFactors=F)
		    validLncStat.df <- t(validLncStat.df)
		    rownames(validLncStat.df) <- validLncTargets.v
                 
		    validLncStat.df <- cbind(validLncStat.df, as.character(names(statLMG.lm)[i]))
		    colnames(validLncStat.df)[9] <- "LncID"
		    validLncStat.df <- cbind(validLncStat.df, rownames(validLncStat.df))
		    colnames(validLncStat.df)[10] <- "GeneID"
		    validLncTargetStat.df <- rbind(validLncTargetStat.df, validLncStat.df )
		}
	}
	#print(validLncTargets.v )

}

rownames(validLncTargetStat.df) <- NULL
validLncTargetStat.df <- validLncTargetStat.df[,c(9,10,1:8)]
save(validLncTargetStat.df, file=paste("./InteractNetworkStat3/", cancers.v[c], "_InteractNetStat.Rdata", sep="") )

}



