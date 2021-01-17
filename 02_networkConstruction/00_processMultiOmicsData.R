
cancers.v  <-  c("BLCA", "BRCA(ER+)", "BRCA(ER-)", "CESC", "CHOL", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC",
				"LUAD", "LUSC", "PAAD" , "PRAD", "READ", "THCA", "UCEC")


for(c in 1:length(cancers.v)) {
print(cancers.v[c])

source("Dolimma.R")

load(paste("./Exp_gene/", cancers.v[c], "_expGenelevel.Rdata", sep=""))
PhenoTypesExp.lv <- PhenoTypes.lv
DEG.m <- Dolimma(PCGexpFPKM.m, as.vector(PhenoTypesExp.lv$Cancer),  "1" , "0")

load(paste("./DNAm_gene/", cancers.v[c], "_DNAmGenelevel.Rdata", sep=""))
PhenoTypesDNAm.lv <- PhenoTypes.lv
DMG.m <- Dolimma(avbeta.m, as.vector(PhenoTypesDNAm.lv$Cancer),  "1" , "0")

load(paste("./LncValExp0/", cancers.v[c], "_valFPKM.Rdata", sep=""))
PhenoTypesLnc.lv <- PhenoTypes.lv
DEL.m <- Dolimma(LncValFPKM.m, as.vector(PhenoTypesLnc.lv$Cancer),  "1" , "0")

#load("./CNV_gene/COAD_cnvGenelevel.Rdata")
load(paste("./CNV_gene/", cancers.v[c], "_cnvGenelevel.Rdata", sep=""))
colnames(CNVscore.m)  <-  paste0(colnames(CNVscore.m), "A") 
interGenes.v <- intersect(intersect(rownames(PCGexpFPKM.m), rownames(avbeta.m)), rownames(CNVscore.m))
interSamples.v <- intersect(intersect(colnames(PCGexpFPKM.m), colnames(avbeta.m)), colnames(CNVscore.m)  )

statEMC.m <- matrix(NA, nrow=length(interGenes.v),ncol=3);
colnames(statEMC.m) <- c("t","P","FDR");
rownames(statEMC.m) <- interGenes.v;

statECM.m <- matrix(NA, nrow=length(interGenes.v),ncol=3);
colnames(statECM.m) <- c("t","P","FDR");
rownames(statECM.m) <- interGenes.v;

for(i in 1:length(interGenes.v) ){
	lm.o <- lm(PCGexpFPKM.m[interGenes.v[i], interSamples.v ] ~ avbeta.m[interGenes.v[i], interSamples.v ] + CNVscore.m[interGenes.v[i], interSamples.v ]);
	statEMC.m[i,1:2] <- summary(lm.o)$coeff[2,3:4];
	statECM.m[i,1:2]  <- summary(lm.o)$coeff[3,3:4];
	#print(i);
}

statEMC.m[,3] <- p.adjust(statEMC.m[,2], method="fdr")
statEMC.m <- statEMC.m[intersect(which(statEMC.m[,3] < 0.05), which(statEMC.m[,1] < 0)) ,]

statECM.m[,3] <- p.adjust(statECM.m[,2], method="fdr")
statECM.m <- statECM.m[intersect(which(statECM.m[,3] < 0.05), which(statECM.m[,1] > 0)) ,]

##
#validGenes.v <- intersect(rownames(DEG.m), rownames(statEMC.m))
validGenes.v <-  rownames(statEMC.m) ## may be consider DEG later
validavBeta.m <- avbeta.m[validGenes.v,]
validExpFPKM.m <- PCGexpFPKM.m[validGenes.v,]

###################

statLMG.lm  <- list()
interSample.v <-  intersect(colnames(LncValFPKM.m), colnames(validavBeta.m))

for (i in 1:nrow(LncValFPKM.m)){
	print(i)
	statLM.m <- matrix(NA, nrow=nrow(validavBeta.m), ncol=4)
	rownames(statLM.m) <- rownames(validavBeta.m)
	colnames(statLM.m) <- c("corM","tM","PM","FDRM")
	statLM.m[,1] <- apply(validavBeta.m[,interSample.v], 1, function(x) cor(x , LncValFPKM.m[i, interSample.v]) )
	statLM.m[,2:3] <- t(apply(validavBeta.m[,interSample.v], 1, function(x) summary(lm(x ~ LncValFPKM.m[i, interSample.v]))$coeff[2,3:4]))
	statLM.m[,4] <- p.adjust(statLM.m[,3], method="fdr")
	
	statLG.m <- matrix(NA, nrow=nrow(validExpFPKM.m), ncol=4)
	rownames(statLG.m) <- rownames(validExpFPKM.m)
	colnames(statLG.m) <- c("corG","tG","PG","FDRG")
	statLG.m[,1] <- apply(validExpFPKM.m, 1, function(x) cor(x , LncValFPKM.m[i,]) )
	statLG.m[,2:3] <- t(apply(validExpFPKM.m, 1, function(x) summary(lm(x ~ LncValFPKM.m[i,]))$coeff[2,3:4]))
	statLG.m[,4] <- p.adjust(statLG.m[,3], method="fdr")

	statLMG.lm[[i]] <- cbind(statLM.m, statLG.m) 
}

names(statLMG.lm) <- rownames(LncValFPKM.m)
save(DEG.m, DMG.m, DEL.m, statEMC.m, statECM.m, statLMG.lm, file=paste("./MultiOmicsData/", cancers.v[c], "_MultiOmicsProcess.Rdata", sep="") )

}


