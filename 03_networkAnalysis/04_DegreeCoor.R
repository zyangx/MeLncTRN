options(stringsAsFactors=FALSE)

cancers.v  <-  c("BLCA", "BRCA(ER+)", "BRCA(ER-)", "CESC", "CHOL", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC",
				"LUAD", "LUSC", "PAAD" , "PRAD", "READ", "THCA", "UCEC")

#mycolors.v <- pal_d3("category20")(18)

##############
Lab.palette <- colorRampPalette(c("blue", "skyblue", "green", "yellow", "orange", "red"), space = "Lab");



pdf("LncDegree_DNAm_cor.pdf", width=10,height=13);
par(mfrow=c(5,4));
par(mar=c(4,4,3,1));

for(c in 1:length(cancers.v)) {

print(cancers.v[c])
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
load(paste("./InteractNetworkStat1/", cancers.v[c], "_InteractNetStat.Rdata", sep="") )

LncDegree.v <- rowSums(InteractLG.m)
tmp.idx <- match(validLncTargetStat.df$LncID, names(LncDegree.v))
#plot(LncDegree.v[tmp.idx], abs(as.numeric(validLncTargetStat.df$corM)), col="skyblue", main="cancer", xlab="Degree", ylab="abs(PCC)", cex.lab=1.2, las=1)
smoothScatter(abs(as.numeric(validLncTargetStat.df$corM)) ~ LncDegree.v[tmp.idx], colramp = Lab.palette, nrpoints=0, main=cancers.v[c], xlab="Degree", ylab="abs(PCC)", cex.lab=1.2, las=1)

lm.o <- lm(abs(as.numeric(validLncTargetStat.df$corM)) ~ LncDegree.v[tmp.idx]);
rsquare=parse(text=paste('R^2==',round(summary(lm.o)$r.sq,2),sep=''));
cor.o <- cor.test(LncDegree.v[tmp.idx], abs(as.numeric(validLncTargetStat.df$corM)) , method="spearman");
#text(x=max(LncDegree.v)*0.65, y=max(abs(as.numeric(validLncTargetStat.df$corM)))*0.95, adj= 0, cex=0.9, font=2, col="#FF0000", rsquare);
text(x=max(LncDegree.v)*0.65, y=max(abs(as.numeric(validLncTargetStat.df$corM)))*0.9, adj= 0, cex=0.9, font=2, col="#FF0000", label= paste('SCC = ',round(cor.o$estimate,digits=2),sep=''));
if(cor.o$p.value == 0){
text(x=max(LncDegree.v)*0.65, y=max(abs(as.numeric(validLncTargetStat.df$corM)))*0.85, adj= 0, cex=0.9, font=2, col="#FF0000", label= 'P < 1.0e-50');
}else{
text(x=max(LncDegree.v)*0.65, y=max(abs(as.numeric(validLncTargetStat.df$corM)))*0.85, adj= 0, cex=0.9, font=2, col="#FF0000", label= paste('P = ',format(cor.o$p.value,digits=3),sep=''));
}

}
dev.off()

######
c <- 1

pdf("LncDegree_DNAm_cor_BLCA.pdf", width=4.5,height=4.5);
par(mar=c(4,4,2,1));

print(cancers.v[c])
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
load(paste("./InteractNetworkStat1/", cancers.v[c], "_InteractNetStat.Rdata", sep="") )

LncDegree.v <- rowSums(InteractLG.m)
tmp.idx <- match(validLncTargetStat.df$LncID, names(LncDegree.v))
#plot(LncDegree.v[tmp.idx], abs(as.numeric(validLncTargetStat.df$corM)), col="skyblue", main="cancer", xlab="Degree", ylab="abs(PCC)", cex.lab=1.2, las=1)
smoothScatter(abs(as.numeric(validLncTargetStat.df$corM)) ~ LncDegree.v[tmp.idx], colramp = Lab.palette, nrpoints=0, xlab="Degree", ylab="abs(PCC)", cex.lab=1.2, las=1)

lm.o <- lm(abs(as.numeric(validLncTargetStat.df$corM)) ~ LncDegree.v[tmp.idx]);
rsquare=parse(text=paste('R^2==',round(summary(lm.o)$r.sq,2),sep=''));
cor.o <- cor.test(LncDegree.v[tmp.idx], abs(as.numeric(validLncTargetStat.df$corM)) , method="spearman");
#text(x=max(LncDegree.v)*0.65, y=max(abs(as.numeric(validLncTargetStat.df$corM)))*0.95, adj= 0, cex=0.9, font=2, col="#FF0000", rsquare);
text(x=max(LncDegree.v)*0.65, y=max(abs(as.numeric(validLncTargetStat.df$corM)))*0.9, adj= 0, cex=0.9, font=2, col="#FF0000", label= paste('SCC = ',round(cor.o$estimate,digits=2),sep=''));
if(cor.o$p.value == 0){
text(x=max(LncDegree.v)*0.65, y=max(abs(as.numeric(validLncTargetStat.df$corM)))*0.85, adj= 0, cex=0.9, font=2, col="#FF0000", label= 'P < 1.0e-50');
}else{
text(x=max(LncDegree.v)*0.65, y=max(abs(as.numeric(validLncTargetStat.df$corM)))*0.85, adj= 0, cex=0.9, font=2, col="#FF0000", label= paste('P = ',format(cor.o$p.value,digits=3),sep=''));
}
dev.off()


##############

pdf("GeneDegree_DNAm_cor.pdf", width=10,height=13);
par(mfrow=c(5,4));
par(mar=c(4,4,3,1));

for(c in 1:length(cancers.v)) {

print(cancers.v[c])
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
load(paste("./InteractNetworkStat1/", cancers.v[c], "_InteractNetStat.Rdata", sep="") )

GeneDegree.v <- colSums(InteractLG.m)
tmp.idx <- match(validLncTargetStat.df$GeneID, names(GeneDegree.v))
#plot(GeneDegree.v[tmp.idx], abs(as.numeric(validLncTargetStat.df$corM)))
smoothScatter(abs(as.numeric(validLncTargetStat.df$corM)) ~ GeneDegree.v[tmp.idx], colramp = Lab.palette, nrpoints=0, main=cancers.v[c], xlab="Degree", ylab="coffecent", cex.lab=1.2)

lm.o <- lm(abs(as.numeric(validLncTargetStat.df$corM)) ~ GeneDegree.v[tmp.idx]);
rsquare=parse(text=paste('R^2==',round(summary(lm.o)$r.sq,2),sep=''));
cor.o <- cor.test(GeneDegree.v[tmp.idx], abs(as.numeric(validLncTargetStat.df$corM)) , method="spearman");
#text(x=max(GeneDegree.v)*0.65, y=max(abs(as.numeric(validLncTargetStat.df$corM)))*0.95, adj= 0, cex=0.9, font=2, col="#FF0000", rsquare);
text(x=max(GeneDegree.v)*0.65, y=max(abs(as.numeric(validLncTargetStat.df$corM)))*0.9, adj= 0, cex=0.9, font=2, col="#FF0000", label= paste('SCC = ',round(cor.o$estimate,digits=2),sep=''));
if(cor.o$p.value == 0){
text(x=max(GeneDegree.v)*0.65, y=max(abs(as.numeric(validLncTargetStat.df$corM)))*0.85, adj= 0, cex=0.9, font=2, col="#FF0000", label= 'P < 1.0e-50');
}else{
text(x=max(GeneDegree.v)*0.65, y=max(abs(as.numeric(validLncTargetStat.df$corM)))*0.85, adj= 0, cex=0.9, font=2, col="#FF0000", label= paste('P = ',format(cor.o$p.value,digits=3),sep=''));
}

}
dev.off()


##############

pdf("LncDegree_Exp_cor.pdf", width=10,height=13);
par(mfrow=c(5,4));
par(mar=c(4,4,3,1));

for(c in 1:length(cancers.v)) {

print(cancers.v[c])
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
load(paste("./InteractNetworkStat1/", cancers.v[c], "_InteractNetStat.Rdata", sep="") )

LncDegree.v <- rowSums(InteractLG.m)
tmp.idx <- match(validLncTargetStat.df$LncID, names(LncDegree.v))
#plot(LncDegree.v[tmp.idx], abs(as.numeric(validLncTargetStat.df$corG)))
smoothScatter(abs(as.numeric(validLncTargetStat.df$corG)) ~ LncDegree.v[tmp.idx], colramp = Lab.palette, nrpoints=0, main=cancers.v[c], xlab="Degree", ylab="abs(PCC)", cex.lab=1.2, las=1)

lm.o <- lm(abs(as.numeric(validLncTargetStat.df$corG)) ~ LncDegree.v[tmp.idx]);
rsquare=parse(text=paste('R^2==',round(summary(lm.o)$r.sq,2),sep=''));
cor.o <- cor.test(LncDegree.v[tmp.idx], abs(as.numeric(validLncTargetStat.df$corG)) , method="spearman");
#text(x=max(LncDegree.v)*0.65, y=max(abs(as.numeric(validLncTargetStat.df$corG)))*0.95, adj= 0, cex=0.9, font=2, col="#FF0000", rsquare);
text(x=max(LncDegree.v)*0.65, y=max(abs(as.numeric(validLncTargetStat.df$corG)))*0.9, adj= 0, cex=0.9, font=2, col="#FF0000", label= paste('SCC = ',round(cor.o$estimate,digits=2),sep=''));
if(cor.o$p.value ==0){
text(x=max(LncDegree.v)*0.65, y=max(abs(as.numeric(validLncTargetStat.df$corG)))*0.85, adj= 0, cex=0.9, font=2, col="#FF0000", label= 'P < 1.0e-50');
}else{
text(x=max(LncDegree.v)*0.65, y=max(abs(as.numeric(validLncTargetStat.df$corG)))*0.85, adj= 0, cex=0.9, font=2, col="#FF0000", label= paste('P = ',format(cor.o$p.value,digits=3),sep=''));
}

}
dev.off()

##
library(plotrix)

c <- 1
pdf("LncDegree_Exp_cor_BLCA_colorlegend.pdf", width=5,height=4.5);
par(mar=c(4,4,2,3.5));

print(cancers.v[c])
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
load(paste("./InteractNetworkStat1/", cancers.v[c], "_InteractNetStat.Rdata", sep="") )

LncDegree.v <- rowSums(InteractLG.m)
tmp.idx <- match(validLncTargetStat.df$LncID, names(LncDegree.v))
#plot(LncDegree.v[tmp.idx], abs(as.numeric(validLncTargetStat.df$corG)))
smoothScatter(abs(as.numeric(validLncTargetStat.df$corG)) ~ LncDegree.v[tmp.idx], colramp = Lab.palette, nrpoints=0,  xlab="Degree", ylab="abs(PCC)", cex.lab=1.2, las=1)

lm.o <- lm(abs(as.numeric(validLncTargetStat.df$corG)) ~ LncDegree.v[tmp.idx]);
rsquare=parse(text=paste('R^2==',round(summary(lm.o)$r.sq,2),sep=''));
cor.o <- cor.test(LncDegree.v[tmp.idx], abs(as.numeric(validLncTargetStat.df$corG)) , method="spearman");
#text(x=max(LncDegree.v)*0.65, y=max(abs(as.numeric(validLncTargetStat.df$corG)))*0.95, adj= 0, cex=0.9, font=2, col="#FF0000", rsquare);
text(x=max(LncDegree.v)*0.65, y=max(abs(as.numeric(validLncTargetStat.df$corG)))*0.9, adj= 0, cex=0.9, font=2, col="#FF0000", label= paste('SCC = ',round(cor.o$estimate,digits=2),sep=''));
if(cor.o$p.value ==0){
text(x=max(LncDegree.v)*0.65, y=max(abs(as.numeric(validLncTargetStat.df$corG)))*0.85, adj= 0, cex=0.9, font=2, col="#FF0000", label= 'P < 1.0e-50');
}else{
text(x=max(LncDegree.v)*0.65, y=max(abs(as.numeric(validLncTargetStat.df$corG)))*0.85, adj= 0, cex=0.9, font=2, col="#FF0000", label= paste('P = ',format(cor.o$p.value,digits=3),sep=''));
}
color.legend(148, 0.4, 154, 0.8, rect.col=Lab.palette(200), legend=c("0.00", "0.01", "0.02", "0.03", "0.04", "0.05"), align="rb",gradient="y")
dev.off()



############
pdf("GeneDegree_Exp_cor.pdf", width=10,height=13);
par(mfrow=c(5,4));
par(mar=c(4,4,3,1));

for(c in 1:length(cancers.v)) {

print(cancers.v[c])
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
load(paste("./InteractNetworkStat1/", cancers.v[c], "_InteractNetStat.Rdata", sep="") )

GeneDegree.v <- colSums(InteractLG.m)
tmp.idx <- match(validLncTargetStat.df$GeneID, names(GeneDegree.v))
#plot(GeneDegree.v[tmp.idx], abs(as.numeric(validLncTargetStat.df$corG)))
smoothScatter(abs(as.numeric(validLncTargetStat.df$corG)) ~ GeneDegree.v[tmp.idx], colramp = Lab.palette, nrpoints=0, main=cancers.v[c], xlab="Degree", ylab="coffecent", cex.lab=1.2)

lm.o <- lm(abs(as.numeric(validLncTargetStat.df$corG)) ~ GeneDegree.v[tmp.idx]);
rsquare=parse(text=paste('R^2==',round(summary(lm.o)$r.sq,2),sep=''));
cor.o <- cor.test(GeneDegree.v[tmp.idx], abs(as.numeric(validLncTargetStat.df$corG)) , method="spearman");
#text(x=max(GeneDegree.v)*0.65, y=max(abs(as.numeric(validLncTargetStat.df$corG)))*0.95, adj= 0, cex=0.9, font=2, col="#FF0000", rsquare);
text(x=max(GeneDegree.v)*0.65, y=max(abs(as.numeric(validLncTargetStat.df$corG)))*0.9, adj= 0, cex=0.9, font=2, col="#FF0000", label= paste('SCC = ',round(cor.o$estimate,digits=2),sep=''));
if(cor.o$p.value == 0){
text(x=max(GeneDegree.v)*0.65, y=max(abs(as.numeric(validLncTargetStat.df$corG)))*0.85, adj= 0, cex=0.9, font=2, col="#FF0000", label= 'P < 1.0e-50');
}else{
text(x=max(GeneDegree.v)*0.65, y=max(abs(as.numeric(validLncTargetStat.df$corG)))*0.85, adj= 0, cex=0.9, font=2, col="#FF0000", label= paste('P = ',format(cor.o$p.value,digits=3),sep=''));
}

}
dev.off()



##############
##############
Lab.palette <- colorRampPalette(c("blue", "skyblue", "green", "yellow", "orange", "red"), space = "Lab");

pdf("LncGeneDegree_Exp_cor.pdf", width=10,height=13);
par(mfrow=c(5,4));
par(mar=c(4,4,3,1));

for(c in 1:length(cancers.v)) {

print(cancers.v[c])
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
load(paste("./InteractNetworkStat1/", cancers.v[c], "_InteractNetStat.Rdata", sep="") )

LncDegree.v <- rowSums(InteractLG.m)
tmp.idx <- match(validLncTargetStat.df$LncID, names(LncDegree.v))
#plot(LncDegree.v[tmp.idx], abs(as.numeric(validLncTargetStat.df$corG)), col="skyblue", main="cancer", xlab="Degree", ylab="abs(PCC)", cex.lab=1.2, las=1)
tmp1.m <- cbind(LncDegree.v[tmp.idx], abs(as.numeric(validLncTargetStat.df$corG)))


GeneDegree.v <- colSums(InteractLG.m)
tmp.idx <- match(validLncTargetStat.df$GeneID, names(GeneDegree.v))
#plot(GeneDegree.v[tmp.idx], abs(as.numeric(validLncTargetStat.df$corG)))
tmp2.m <- cbind(GeneDegree.v[tmp.idx], abs(as.numeric(validLncTargetStat.df$corG )))


data.m <- rbind(tmp1.m, tmp2.m)

smoothScatter(data.m[,2] ~ data.m[,1], colramp = Lab.palette, nrpoints=0, main=cancers.v[c], xlab="Degree", ylab="abs(PCC)", cex.lab=1.2, las=1)


lm.o <- lm(data.m[,2] ~ data.m[,1]);
rsquare=parse(text=paste('R^2==',round(summary(lm.o)$r.sq,2),sep=''));
cor.o <- cor.test(data.m[,2], data.m[,1], method="spearman");
#text(x=max(data.m[,1])*0.65, y=max(data.m[,2])*0.95, adj= 0, cex=0.9, font=2, col="#FF0000", rsquare);
text(x=max(data.m[,1])*0.65, y=max(data.m[,2])*0.9, adj= 0, cex=0.9, font=2, col="#FF0000", label= paste('SCC = ',round(cor.o$estimate,digits=2),sep=''));
if(cor.o$p.value  == 0){
text(x=max(data.m[,1])*0.65, y=max(data.m[,2])*0.85, adj= 0, cex=0.9, font=2, col="#FF0000", label= 'P < 1.0e-50');
}else{
text(x=max(data.m[,1])*0.65, y=max(data.m[,2])*0.85, adj= 0, cex=0.9, font=2, col="#FF0000", label= paste('P = ',format(cor.o$p.value,digits=3),sep=''));
}

}
dev.off()



##############
##############
#pdf("LncGeneDegree_Exp_Poscor.pdf", width=10,height=13);
pdf("LncGeneDegree_Exp_Negcor.pdf", width=10,height=13);
par(mfrow=c(5,4));
par(mar=c(4,4,3,1));

for(c in 1:length(cancers.v)) {

print(cancers.v[c])
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
load(paste("./InteractNetworkStat1/", cancers.v[c], "_InteractNetStat.Rdata", sep="") )

LncDegree.v <- rowSums(InteractLG.m)
tmp.idx <- match(validLncTargetStat.df$LncID, names(LncDegree.v))
#plot(LncDegree.v[tmp.idx], as.numeric(validLncTargetStat.df$corG), col="skyblue", main="cancer", xlab="Degree", ylab="PCC", cex.lab=1.2, las=1)
tmp1.m <- cbind(LncDegree.v[tmp.idx], as.numeric(validLncTargetStat.df$corG))

GeneDegree.v <- colSums(InteractLG.m)
tmp.idx <- match(validLncTargetStat.df$GeneID, names(GeneDegree.v))
#plot(GeneDegree.v[tmp.idx], as.numeric(validLncTargetStat.df$corG))
tmp2.m <- cbind(GeneDegree.v[tmp.idx], as.numeric(validLncTargetStat.df$corG ))

data.m <- rbind(tmp1.m, tmp2.m)
#data.m <- data.m[which(data.m[,2] > 0),]
data.m <- data.m[which(data.m[,2] < 0),]

smoothScatter(data.m[,2] ~ data.m[,1], colramp = Lab.palette, nrpoints=0, main=cancers.v[c], xlab="Degree", ylab="PCC", cex.lab=1.2, las=1)


lm.o <- lm(data.m[,2] ~ data.m[,1]);
rsquare=parse(text=paste('R^2==',round(summary(lm.o)$r.sq,2),sep=''));
cor.o <- cor.test(data.m[,2], data.m[,1], method="spearman");
#text(x=max(data.m[,1])*0.65, y=max(data.m[,2])*0.95, adj= 0, cex=0.9, font=2, col="#FF0000", rsquare);
text(x=max(data.m[,1])*0.65, y=min(data.m[,2])*0.85, adj= 0, cex=0.9, font=2, col="#FF0000", label= paste('SCC = ',round(cor.o$estimate,digits=2),sep=''));
if(cor.o$p.value  == 0){
text(x=max(data.m[,1])*0.65, y=min(data.m[,2])*0.9, adj= 0, cex=0.9, font=2, col="#FF0000", label= 'P < 1.0e-50');
}else{
text(x=max(data.m[,1])*0.65, y=min(data.m[,2])*0.9, adj= 0, cex=0.9, font=2, col="#FF0000", label= paste('P = ',format(cor.o$p.value,digits=3),sep=''));
}

}
dev.off()






