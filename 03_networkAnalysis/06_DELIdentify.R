

cancers.v  <-  c("BLCA", "BRCA(ER+)", "BRCA(ER-)", "CESC", "CHOL", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC",
				"LUAD", "LUSC", "PAAD" , "PRAD", "READ", "THCA", "UCEC")


for(c in 1:length(cancers.v)) {
print(cancers.v[c])

load(paste("./LncValExp0/", cancers.v[c], "_valFPKM.Rdata", sep=""))
PhenoTypesLnc.lv <- PhenoTypes.lv

DEL.m <- matrix(NA, nrow=nrow(LncValFPKM.m), ncol=6 )
rownames(DEL.m) <- rownames(LncValFPKM.m)
colnames(DEL.m) <- c("meanC","meanN","logFC","T-stat","P_value1","FDR1")

LncSwitchStat.m <- matrix(NA, nrow=nrow(LncValFPKM.m), ncol=5 )
rownames(LncSwitchStat.m) <- rownames(LncValFPKM.m)
colnames(LncSwitchStat.m) <- c("FreqC","FreqN","FreqFold","P_value2","FDR2")

for(i in 1:nrow(LncValFPKM.m)){

t.o <- t.test(LncValFPKM.m[i,] ~ as.vector(PhenoTypesLnc.lv$Cancer))
DEL.m[i,1] <-  t.o$estimate[2]
DEL.m[i,2] <-  t.o$estimate[1]
DEL.m[i,3] <-  mean(LncValFPKM.m[i, which(PhenoTypesLnc.lv$Cancer == "1")] ) - mean(LncValFPKM.m[i, which(PhenoTypesLnc.lv$Cancer == "0")] ) 
DEL.m[i,4] <-  t.o$statistic
DEL.m[i,5] <-  t.o$p.value

status.v <- ifelse(LncValFPKM.m[i,] > 0, 1, 0)
cancer.v <- PhenoTypesLnc.lv$Cancer
if(length(which(status.v == 0)) > 0 ){
tmp.m <- table(status.v, cancer.v )
fisher.o <- fisher.test(tmp.m)
LncSwitchStat.m[i,1] <-  tmp.m[1,2] / tmp.m[2,2]
LncSwitchStat.m[i,2] <-  tmp.m[1,1] / tmp.m[2,1]
LncSwitchStat.m[i,3] <-  LncSwitchStat.m[i,2] / LncSwitchStat.m[i,1] 
LncSwitchStat.m[i,4] <-  fisher.o$p.value
}


}
DEL.m[,6] <- p.adjust(DEL.m[,5], method="BH")
LncSwitchStat.m[,5] <- p.adjust(LncSwitchStat.m[,4], method="BH")

LncDEStatus.m <- cbind(DEL.m, LncSwitchStat.m)
##########
save(LncDEStatus.m, file=paste("./LncDEstatus/", cancers.v[c], "_LncDEStatus.Rdata", sep="") )

}

break;
tmp1.idx <-  intersect( sort(c(which(DEL.m[,3] > 1), which(DEL.m[,3] < -1))), which(DEL.m[,6] <0.05) )
tmp2.idx <-  intersect( sort(c(which(LncSwitchStat.m[,3] > 2), which(LncSwitchStat.m[,3] < 0.5))), which(LncSwitchStat.m[,5] <0.05) )


