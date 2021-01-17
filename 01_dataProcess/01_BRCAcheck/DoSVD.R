### DoSVD.R
GetSVDcorPv <- function(input.m, PhenoTypes.lv){
	tmp.m <- input.m - rowMeans(input.m);
	svd.o <- svd(tmp.m);

	library(isva);
	rmt.o <- EstDimRMT(tmp.m);
	topPCA <- rmt.o$dim;

	selcat.idx <- 1:length(PhenoTypes.lv);

	corPv.m <- matrix(nrow=topPCA, ncol=length(selcat.idx));
	colnames(corPv.m) <- names(PhenoTypes.lv)[selcat.idx];
	rownames(corPv.m) <- paste("PC-",1:nrow(corPv.m),sep="")

	factor.log <- c();
	for(i in selcat.idx){
		factor.log <- c(factor.log, is.factor(PhenoTypes.lv[[i]]));
	}

	for(c in 1:topPCA){
		tmp.v <- svd.o$v[,c];
		for(f in 1:length(selcat.idx)){
			if(factor.log[f]){
				corPv.m[c,f] <- kruskal.test(tmp.v ~ as.factor(PhenoTypes.lv[[selcat.idx[f]]]))$p.value;
			}
			else{
				corPv.m[c,f] <- summary(lm(tmp.v ~ PhenoTypes.lv[[selcat.idx[f]]]))$coeff[2,4];
			}
		}
		#print(c);
	}

	tmp.v <- svd.o$d[1:topPCA];
	f.v <- tmp.v^2/sum(svd.o$d^2);
	names(f.v) <- paste("PC-",1:nrow(corPv.m),sep="")

	svdCorPv.lv <- list(corPv.m, f.v)
	names(svdCorPv.lv) <- c("svdCorPv.m", "fraction")

	return(svdCorPv.lv)
}


DrawPvHeatmap  <- function(fraction.v,svdCorPv.m){
	Palette.v <- c("darkred","red","orange","pink","white");
	breaks.v <- c(-200,-10,-5,-2,log10(0.05),0);

	pdf("SVDheatmap.pdf",width=6,height=7 );

	layout(matrix(c(1,2,2,2), nr=4, byrow=TRUE))
	par(family="ArialMT")
	par(mar=c(5,10,2,1));


	plot(fraction.v,ylab="fVAR",xlab="PC",pch=23,type="b",lwd=2,col="red",cex.axis=1.5,cex.lab=1.5);
	#text(par("usr")[1]-6.2 , 0.6, adj = 0, labels = "A)", xpd = TRUE, cex=1.8)

	image(x=1:nrow(svdCorPv.m),y=1:ncol(svdCorPv.m),z=svdCorPv.m,col=Palette.v,breaks=breaks.v,xlab="",ylab="",axes=FALSE);
	#axis(1,at=1:nrow(svdCorPv.m),labels=rownames(svdCorPv.m),las=2);
	axis(1,at=seq(1, nrow(svdCorPv.m), by=5),labels=rownames(svdCorPv.m)[seq(1, nrow(svdCorPv.m), by=5)],las=2);
	axis(2,at=1:ncol(svdCorPv.m),labels=colnames(svdCorPv.m),las=2.5,cex.axis=1.1);
	#text(par("usr")[1] - 1 , 25, adj = 0, labels = "B)", xpd = TRUE, cex=1.8)

	dev.off();
}


DrawPvHeatmap1  <- function(fraction.v,svdCorPv.m){
	library(reshape2)
	library(ggplot2)
	library(scales)

	#Palette.v <- c("darkred","red","orange","pink","white");
	#breaks.v <- c(-100,-10,-5,-2,log10(0.05),0);
	Palette.v <- c("darkred","red","orange","pink","white", "white");
	#breaks.v <- c(-100,-10,-5,-2,log10(0.05),log10(0.0499999),0);
	breaks.v <- c(0, 1E-10, 1E-5, 1E-2, 0.05,0.05000001, 1);

	#svdCorPv.df <- as.data.frame(t(log10(LNCsvdCorPv.o$svdCorPv.m)))
	svdCorPv.df <- as.data.frame(t(LNCsvdCorPv.o$svdCorPv.m))
	svdCorPv.df$Pheno <- rownames(svdCorPv.df)
	data.df <- melt(svdCorPv.df, id.vars=c("Pheno"))

	#fraction.v 

	ggheatmap <- ggplot(data.df, aes(x=variable, y=Pheno)) + geom_tile(aes(fill = value))
	#ggheatmap <- ggheatmap + scale_fill_gradientn(colours=Palette.v, values=rescale(breaks.v))
	ggheatmap <- ggheatmap + scale_fill_gradientn(colours=Palette.v, values=breaks.v, minor_breaks=NULL)
	ggheatmap <- ggheatmap + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.position="top")
	ggheatmap <- ggheatmap + guides(fill=FALSE)
	ggsave(ggheatmap, filename="SVDheatmapLNC1.pdf", width=6, height=6, units=c("cm"),colormodel="srgb")

}

