####### DoSVD and draw correlation p-value heatmap ######
source("DoSVD.R")


load("BRCA_PCG_expFPKM_PhenoSupp.Rdata")
svdCorPv.o <- GetSVDcorPv(PCGexpFPKM.m, PhenoTypesBRCA.lv);
f.v <- svdCorPv.o$fraction
svdCorPv.m <- log10(svdCorPv.o$svdCorPv.m)
DrawPvHeatmap(f.v , svdCorPv.m)


load("BRCA_expFPKM_PhenoSupp.Rdata")
svdCorPv.o <- GetSVDcorPv(expFPKM.m, PhenoTypesBRCA.lv);
f.v <- svdCorPv.o$fraction
svdCorPv.m <- log10(svdCorPv.o$svdCorPv.m)
DrawPvHeatmap(f.v , svdCorPv.m)


load("BRCA_bmiqbeta_PhenoSupp.Rdata")
DNAmsvdCorPv.o <- GetSVDcorPv(bmiqbeta.m, PhenoTypesBRCA.lv);
f.v <- DNAmsvdCorPv.o$fraction
svdCorPv.m <- log10(DNAmsvdCorPv.o$svdCorPv.m)
DrawPvHeatmap(f.v , svdCorPv.m)



### image heatmap
#library(marray);
Palette.v <- c("darkred","red","orange","pink","white");
breaks.v <- c(-200,-10,-5,-2,log10(0.05),0);

pdf("SVDheatmapLNC.pdf",width=6,height=7);
layout(matrix(c(1,2,2,2), nr=4, byrow=TRUE))
par(mar=c(5,10,2,1));

plot(f.v,ylab="fVAR",xlab="PC",pch=23,type="b",lwd=2,col="red",cex.axis=1.5,cex.lab=1.5);
text(par("usr")[1]-6.2 , 0.6, adj = 0, labels = "A)", xpd = TRUE, cex=1.8)

image(x=1:nrow(svdCorPv.m),y=1:ncol(svdCorPv.m),z=svdCorPv.m,col=Palette.v,breaks=breaks.v,xlab="Pheno",ylab="",axes=FALSE);
axis(1,at=1:nrow(svdCorPv.m),labels=rownames(svdCorPv.m),las=2);
axis(2,at=1:ncol(svdCorPv.m),labels=colnames(svdCorPv.m),las=2.5,cex.axis=1.1);
text(par("usr")[1] - 6 , 25, adj = 0, labels = "B)", xpd = TRUE, cex=1.8)
dev.off();


library(reshape2)
library(ggplot2)
library(scales)

#Palette.v <- c("darkred","red","orange","pink","white");
#breaks.v <- c(-100,-10,-5,-2,log10(0.05),0);

Palette.v <- c("darkred","red","orange","pink","white", "white");
breaks.v <- c(-100,-10,-5,-2,log10(0.05),log10(0.0499999),0);


f.v <- LNCsvdCorPv.o$fraction
svdCorPv.df <- as.data.frame(t(log10(LNCsvdCorPv.o$svdCorPv.m)))
svdCorPv.df$Pheno <- rownames(svdCorPv.df)
data.df <- melt(svdCorPv.df, id.vars=c("Pheno"))
ggheatmap <- ggplot(data.df, aes(x=variable, y=Pheno)) + geom_tile(aes(fill = value))
ggheatmap <- ggheatmap + scale_fill_gradientn(colours=Palette.v, values=rescale(breaks.v))
ggheatmap <- ggheatmap + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.position="top")
ggheatmap <- ggheatmap + guides(fill=FALSE)
ggsave(ggheatmap, filename="SVDheatmapLNC1.pdf", width=6, height=6, units=c("cm"),colormodel="srgb")


breaks.v <- c(0, 1E-10, 1E-5, 1E-2, 0.05,0.05000001,1);
svdCorPv.df <- as.data.frame(t(LNCsvdCorPv.o$svdCorPv.m))
svdCorPv.df$Pheno <- rownames(svdCorPv.df)
data.df <- melt(svdCorPv.df, id.vars=c("Pheno"))
ggheatmap <- ggplot(data.df, aes(x=variable, y=Pheno)) + geom_tile(aes(fill = value))
ggheatmap <- ggheatmap + scale_fill_gradientn(colours=Palette.v, values=breaks.v, minor_breaks=NULL)
ggheatmap <- ggheatmap + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.position="top")
ggheatmap <- ggheatmap + guides(fill=FALSE)
ggsave(ggheatmap, filename="SVDheatmapLNC2.pdf", width=6, height=6, units=c("cm"),colormodel="srgb")





