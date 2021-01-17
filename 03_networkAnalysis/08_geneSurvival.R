options(stringsAsFactors=FALSE)

library(survival)
library(forestplot)
library(survminer)


cancers.v  <-  c("BLCA", "BRCA(ER+)", "BRCA(ER-)", "CESC", "CHOL", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC",
				"LUAD", "LUSC", "PAAD" , "PRAD", "READ", "THCA", "UCEC")

ComputeHR <- function(cox.o) {tmp.v <- c("Coeff" = summary(cox.o)$coeff[1,1],
					"Zscore" = summary(cox.o)$coeff[1,4],
					"HR"=summary(cox.o)$coeff[1,2],
					"lower95"=summary(cox.o)$conf[1,3],
					"upper95"=summary(cox.o)$conf[1,4],
					summary(cox.o)$sctest[3]) }


for(c in 1:length(cancers.v)) {
print(cancers.v[c])

load(paste("./Exp_gene/", cancers.v[c], "_expGenelevel.Rdata", sep=""))
PhenoTypesExp.lv <- PhenoTypes.lv

load(paste("./LncValExp0/", cancers.v[c], "_valFPKM.Rdata", sep=""))
PhenoTypesLnc.lv <- PhenoTypes.lv

survival.df <- data.frame("Time" = as.numeric(PhenoTypes.lv$Overall_Survival / 365.25 ), "Event" = as.numeric(as.vector(PhenoTypes.lv$Vital_status) ) )
survival.df$Time[which(survival.df$Time > 10)] <-  10
rownames(survival.df) <- colnames(PCGexpFPKM.m)

cancer.idx <- which(PhenoTypes.lv[[1]] == 1)
PCGexpFPKM.m <- PCGexpFPKM.m[,cancer.idx]
LncValFPKM.m <- LncValFPKM.m[,cancer.idx]
survival.df <- survival.df[cancer.idx,]

train.idx <- sample(x=1:ncol(PCGexpFPKM.m), size=ceiling(ncol(PCGexpFPKM.m)/2), replace=F)
test.idx <- setdiff(1:ncol(PCGexpFPKM.m), train.idx)
#save(train.idx, test.idx, file="randomSample.Rdata")

PCGexpFPKMTrain.m <- PCGexpFPKM.m[,train.idx]
PCGexpFPKMTest.m <- PCGexpFPKM.m[,test.idx]

LncValFPKMTrain.m <- LncValFPKM.m[,train.idx]
LncValFPKMTest.m <- LncValFPKM.m[,test.idx]

survivalTrain.df <- survival.df[train.idx,]
survivalTest.df <- survival.df[test.idx,]

##### normalization input matrix ######
tmp.m <- PCGexpFPKMTrain.m-rowMeans(PCGexpFPKMTrain.m);
tmp.m <- tmp.m/sqrt(apply(tmp.m,1,var));

cox_results.m <- matrix(NA, nrow=nrow(PCGexpFPKMTrain.m), ncol=6)
rownames(cox_results.m) <- rownames(PCGexpFPKMTrain.m)
colnames(cox_results.m) <- c("Coeff", "Zscore",  "HR",  "HR.95L",  "HR.95H",  "P")

for(i in 1:nrow(tmp.m)){
  gene= tmp.m[i,]
  if(length(which(is.nan(gene))) == 0){
    datatmp <- cbind(survivalTrain.df[,1:2], gene)
    datatmp$group <- ifelse(gene > median(gene),'high','low') 
    #sfit.o <- survfit(Surv( Time, Event ) ~ group, data=datatmp)
    #sdiff.o <- survdiff(Surv(Time, Event) ~ group, data = datatmp)
    #cox.o <- coxph(Surv( Time, Event ) ~ group, data=datatmp)
    cox.o <- coxph(Surv( Time, Event ) ~ gene, data=datatmp)
    coeff.v <- ComputeHR(cox.o)
    cox_results.m[i,] <- round(coeff.v, digits=3)
    #p.val = 1 - pchisq(sdiff.o$chisq, length(sdiff.o$n) - 1)
    #cox_results.m[i,6] <- p.val
  }
}
PCG_cox_results.m <- cox_results.m

##### normalization input matrix ######
tmp.m <- LncValFPKMTrain.m-rowMeans(LncValFPKMTrain.m);
tmp.m <- tmp.m/sqrt(apply(tmp.m,1,var));

cox_results.m <- matrix(NA, nrow=nrow(LncValFPKMTrain.m), ncol=6)
rownames(cox_results.m) <- rownames(LncValFPKMTrain.m)
colnames(cox_results.m) <- c("Coeff", "Zscore",  "HR",  "HR.95L",  "HR.95H",  "P")

for(i in 1:nrow(tmp.m)){
  gene= tmp.m[i,]
  if(length(which(is.nan(gene))) == 0){
    datatmp <- cbind(survivalTrain.df[,1:2], gene)
    datatmp$group <- ifelse(gene > median(gene),'high','low') 

    #sfit.o <- survfit(Surv( Time, Event ) ~ group, data=datatmp)
    #sdiff.o <- survdiff(Surv(Time, Event) ~ group, data = datatmp)
    #cox.o <- coxph(Surv( Time, Event ) ~ group, data=datatmp)
    cox.o <- coxph(Surv( Time, Event ) ~ gene, data=datatmp)
    coeff.v <- ComputeHR(cox.o)
    cox_results.m[i,] <- round(coeff.v, digits=3)
    #p.val = 1 - pchisq(sdiff.o$chisq, length(sdiff.o$n) - 1)
    #cox_results.m[i,6] <- p.val
  }
}
Lnc_cox_results.m <- cox_results.m

save(PCGexpFPKMTrain.m, PCGexpFPKMTest.m, LncValFPKMTrain.m, LncValFPKMTest.m, survivalTrain.df, survivalTest.df, PCG_cox_results.m, Lnc_cox_results.m, train.idx, test.idx, file=paste("./Survival_Lnc_Gene/", cancers.v[c], "_Survival.Rdata", sep=""))

}


