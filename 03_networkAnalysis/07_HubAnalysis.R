options(stringsAsFactors=FALSE)

load("../CancerRdata/exp_FPKM/BLCA_LNC_expFPKM.Rdata")
load("../CancerRdata/Exp_gene/BLCA_expGenelevel.Rdata")


cancers.v  <-  c("BLCA", "BRCA(ER+)", "BRCA(ER-)", "CESC", "CHOL", "COAD", "ESCA", "HNSC",
				 "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PAAD"  , "PRAD", "READ",  "THCA", "UCEC")



LncHubCount.m <- matrix(0, nrow=nrow(LNCexpFPKM.m), ncol=18)
rownames(LncHubCount.m) <- rownames(LNCexpFPKM.m)
colnames(LncHubCount.m) <- cancers.v

PCGHubCount.m <- matrix(0, nrow=nrow(PCGexpFPKM.m), ncol=18)
rownames(PCGHubCount.m) <- rownames(PCGexpFPKM.m)
colnames(PCGHubCount.m) <- cancers.v

LncHubMarker.m <- matrix(0, nrow=nrow(LNCexpFPKM.m), ncol=18)
rownames(LncHubMarker.m) <- rownames(LNCexpFPKM.m)
colnames(LncHubMarker.m) <- cancers.v

PCGHubMarker.m <- matrix(0, nrow=nrow(PCGexpFPKM.m), ncol=18)
rownames(PCGHubMarker.m) <- rownames(PCGexpFPKM.m)
colnames(PCGHubMarker.m) <- cancers.v

for(c in 1:length(cancers.v)) {

	print(cancers.v[c])
	load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
	load(paste("./InteractNetworkStat1/", cancers.v[c], "_InteractNetStat.Rdata", sep="") )

	LncHubCount.v <- rowSums(InteractLG.m) 
	tmp.idx <- match(names(LncHubCount.v), rownames(LncHubCount.m))
	LncHubCount.m[tmp.idx,c] <- LncHubCount.v 

	LncHubs.v <- sort(LncHubCount.v,decreasing=T)[1:round(length(LncHubCount.v) * 0.1)]
	tmp.idx <- match(names(LncHubs.v), rownames(LncHubMarker.m))
	LncHubMarker.m[tmp.idx,c] <- 1 


	PCGHubCount.v <- colSums(InteractLG.m) 
	tmp.idx <- match(names(PCGHubCount.v), rownames(PCGHubCount.m))
	PCGHubCount.m[tmp.idx,c] <- PCGHubCount.v

	PCGHubs.v <- sort(PCGHubCount.v,decreasing=T)[1:round(length(PCGHubCount.v) * 0.1)]
	tmp.idx <- match(names(PCGHubs.v), rownames(PCGHubMarker.m))
	PCGHubMarker.m[tmp.idx,c] <- 1 

}



library(RColorBrewer)
library(extrafont)
#font_import()
library(ggsci)
library(scales)
mycolors.v <- pal_d3("category20")(18)


pdf("Lnc_Degree_Cumulative_Distribution_plot_0.1.pdf", width=4.5, height=4.5)
par(mar=c(4,4,2,1))
plot(ecdf(LncHubCount.m[,1]), ylim=c(0.5,1), main="", xlim=c(1,600), ylab="CDF", xlab="Degree",lwd = 2 , las=1, cex.axis=0.75, cex.lab=0.75, col=mycolors.v[1], do.points=FALSE, verticals=TRUE)
#plot(ecdf(setdiff(LncHubCount.m[,1], 0)), ylim=c(0,1), main="", xlim=c(1,1200), ylab="CDF", xlab="Degree",lwd = 2 , las=1, cex.axis=0.75, cex.lab=0.75, col=mycolors.v[1], do.points=FALSE, verticals=TRUE)

grid(nx=20,ny=20)
for(c in 2:ncol(LncHubCount.m)){
lines(ecdf(LncHubCount.m[,c]),col =mycolors.v[c], lwd = 2 , do.points=FALSE, verticals=TRUE)
#lines(ecdf(setdiff(LncHubCount.m[,c], 0)), col =mycolors.v[c], lwd = 2 , do.points=FALSE, verticals=TRUE)
}
legend("bottomright", lty=1, legend=cancers.v, col=mycolors.v, cex=0.5)
dev.off()


pv1 <- ks.test(LncHubCount.m[,6], LncHubCount.m[,16])$p.v
pv2 <- ks.test(LncHubCount.m[,9], LncHubCount.m[,10])$p.v
pv3 <- ks.test(LncHubCount.m[,12], LncHubCount.m[,13])$p.v




pdf("PCG_Degree_Cumulative_Distribution_plot_0.1.pdf", width=4.5, height=4.5)
par(mar=c(4,4,2,1))
plot(ecdf(PCGHubCount.m[,1]), ylim=c(0.5,1), main="", xlim=c(1,600), ylab="CDF", xlab="Degree",las=1, cex.axis=0.75, cex.lab=0.75, col=mycolors.v[1], do.points=FALSE, verticals=TRUE)
#plot(ecdf(setdiff(PCGHubCount.m[,1],0)), ylim=c(0,1), main="", xlim=c(1,1200), ylab="CDF", xlab="Degree",las=1, cex.axis=0.75, cex.lab=0.75, col=mycolors.v[1], do.points=FALSE, verticals=TRUE)

grid(nx=20,ny=20)
for(c in 2:ncol(PCGHubCount.m)){
lines(ecdf(PCGHubCount.m[,c]),col =mycolors.v[c], lwd = 2 , do.points=FALSE, verticals=TRUE)
#lines(ecdf(setdiff(PCGHubCount.m[,c],0)),col =mycolors.v[c], lwd = 2 , do.points=FALSE, verticals=TRUE)

}
legend("bottomright", lty=1, legend=cancers.v, col=mycolors.v, cex=0.5)
dev.off()




#######
#######
tmp.idx <- which(rowSums(LncHubMarker.m) > 0)
LncHubMarker.m <- LncHubMarker.m[tmp.idx,]
LncHubCount.m <- LncHubCount.m[tmp.idx,]


order.idx <- c()
for(c in 1:length(cancers.v)) {

	print(cancers.v[c])
	tmp.idx <- which(LncHubMarker.m[,c] > 0)
	tmp.idx <- setdiff(tmp.idx, order.idx)
	order.idx <- c(order.idx, tmp.idx)

}
order.idx <- rev(order.idx)

LncHubMarker.m <- LncHubMarker.m[order.idx,]

pdf("Lnc_Hub_Pan_cancer_consistent_plot_0.1.pdf", width=6, height=6.5)
layout(matrix(c(1,1,1,1,2), nr=1, byrow=TRUE))
par(mar=c(5,3,1.8,0.5))
image(x=1:ncol(LncHubMarker.m),y=1:nrow(LncHubMarker.m),z=t(LncHubMarker.m),col=c("white", "blue"), xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", xlab="",ylab="",axes=T,  bty="o")
grid(nx=18,ny=12)
axis(1,at=1:ncol(LncHubMarker.m),labels=colnames(LncHubMarker.m),las=2,cex.axis=0.9, srt = 45);
axis(2, at=NULL, labels=FALSE, tick=FALSE, lwd = 50)
axis(3, at=NULL, labels=FALSE, tick=FALSE, lwd = 2)
axis(4, at=NULL, labels=FALSE, tick=FALSE, lwd = 5)
#text(seq(1, 18, by=1)-0.5, par("usr")[3] - 70, labels = cancers.v, srt = 90, pos = 1, adj =1.1, xpd = TRUE,cex=0.9,font=1)
mtext("Hub Lnc", side=2, line=1, cex=1,font=1)
#axis(2,at=1:nrow(LncHubMarker.m),labels=rownames(LncHubMarker.m),las=2);
legend(x =par("usr")[1]+3 , y = 500, legend = c("Hub", "Non Hub"), fill = c("blue", "white"),horiz=F, bty = "n", xpd = NA)
par(mar=c(3.1,0.1,0,1.2))
barplot(as.numeric(rowSums(LncHubMarker.m)), horiz=T, yaxt="n", col="green", xlab=c("Number of Cancer Types"))
dev.off()


#######
tmp.idx <- which(rowSums(PCGHubMarker.m) > 0)
PCGHubMarker.m <- PCGHubMarker.m[tmp.idx,]
PCGHubCount.m <- PCGHubCount.m[tmp.idx,]

#PCGHubMarker <- PCGHubMarker.m[tmp.idx,]
order.idx <- c()
for(c in 1:length(cancers.v)) {

	print(cancers.v[c])
	tmp.idx <- which(PCGHubMarker.m[,c] > 0)
	tmp.idx <- setdiff(tmp.idx, order.idx)
	order.idx <- c(order.idx, tmp.idx)

}
order.idx <- rev(order.idx)

PCGHubMarker.m <- PCGHubMarker.m[order.idx,]

pdf("PCG_Hub_Pan_cancer_consistent_plot_0.1.pdf", width=6, height=6.5)
layout(matrix(c(1,1,1,1,2), nr=1, byrow=TRUE))
par(mar=c(5,3,1.8,0.5))
image(x=1:ncol(PCGHubMarker.m),y=1:nrow(PCGHubMarker.m),z=t(PCGHubMarker.m),col=c("white", "red"), xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", xlab="",ylab="",axes=T,  bty="o")
grid(nx=18,ny=12)
axis(1,at=1:ncol(PCGHubMarker.m),labels=colnames(PCGHubMarker.m),las=2,cex.axis=0.9, srt = 45);
axis(2, at=NULL, labels=FALSE, tick=FALSE, lwd = 50)
axis(3, at=NULL, labels=FALSE, tick=FALSE, lwd = 2)
axis(4, at=NULL, labels=FALSE, tick=FALSE, lwd = 5)
#text(seq(1, 18, by=1)-0.5, par("usr")[3] - 70, labels = cancers.v, srt = 90, pos = 1, adj =1.1, xpd = TRUE,cex=0.9,font=1)
mtext("Hub Gene", side=2, line=1, cex=1,font=1)
#axis(2,at=1:nrow(PCGHubMarker.m),labels=rownames(PCGHubMarker.m),las=2);
legend(x =par("usr")[1]+3 , y = 150, legend = c("Hub", "Non Hub"), fill = c("red", "white"),horiz=F, bty = "n", xpd = NA)
par(mar=c(3.1,0.1,0,1.2))
barplot(as.numeric(rowSums(PCGHubMarker.m)), horiz=T, yaxt="n", col="green", xlab=c("Number of Cancer Types"))
dev.off()



#######################
#######################

dataEL.df <- data.frame(Exp=0, Type=0, Cancer=0)
dataEL.df <- dataEL.df[-1,]

dataEP.df <- data.frame(Exp=0, Type=0, Cancer=0)
dataEP.df <- dataEP.df[-1,]

LncHubCount.m <- matrix(NA, nrow=nrow(LNCexpFPKM.m), ncol=18)
rownames(LncHubCount.m) <- rownames(LNCexpFPKM.m)
colnames(LncHubCount.m) <- cancers.v

PCGHubCount.m <- matrix(NA, nrow=nrow(PCGexpFPKM.m), ncol=18)
rownames(PCGHubCount.m) <- rownames(PCGexpFPKM.m)
colnames(PCGHubCount.m) <- cancers.v

LncHubMarker.m <- matrix(NA, nrow=nrow(LNCexpFPKM.m), ncol=18)
rownames(LncHubMarker.m) <- rownames(LNCexpFPKM.m)
colnames(LncHubMarker.m) <- cancers.v

PCGHubMarker.m <- matrix(NA, nrow=nrow(PCGexpFPKM.m), ncol=18)
rownames(PCGHubMarker.m) <- rownames(PCGexpFPKM.m)
colnames(PCGHubMarker.m) <- cancers.v

for(c in 1:length(cancers.v)) {

	print(cancers.v[c])
	load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
	load(paste("./InteractNetworkStat1/", cancers.v[c], "_InteractNetStat.Rdata", sep="") )
	load(paste("../CancerRdata/exp_FPKM/", cancers.v[c], "_LNC_expFPKM.Rdata", sep=""))
	load(paste("../CancerRdata/Exp_gene/", cancers.v[c], "_expGenelevel.Rdata", sep=""))

	LncHubCount.v <- rowSums(InteractLG.m) 
	tmp.idx <- match(names(LncHubCount.v), rownames(LncHubCount.m))
	LncHubCount.m[tmp.idx,c] <- LncHubCount.v 

	threshold <- quantile(LncHubCount.v, probs = c(90)/100)
	LncHubMarker.v <- ifelse(LncHubCount.v > threshold, "Hub", "Non Hub")

	#LncHubs.v <- LncHubMarker.v[which(LncHubMarker.v == 1)]
	tmp.idx <- match(names(LncHubMarker.v), rownames(LncHubMarker.m))
	LncHubMarker.m[tmp.idx,c] <- LncHubMarker.v

	LncExp.v <-  rowMeans(LNCexpFPKM.m)
	tmp.idx <- match(names(LncHubMarker.v),  names(LncExp.v))

	data.df <- as.data.frame(cbind(LncExp.v[tmp.idx], LncHubMarker.v ))
	colnames(data.df) <- c("Exp", "Type")
	data.df$Exp <- as.numeric(data.df$Exp)
	data.df <- within(data.df, Type <- factor(Type, levels = c("Hub", "Non Hub")))
	data.df$Cancer <- cancers.v[c]
	dataEL.df <- rbind(dataEL.df, data.df)

        ######
	PCGHubCount.v <- colSums(InteractLG.m) 
	tmp.idx <- match(names(PCGHubCount.v), rownames(PCGHubCount.m))
	PCGHubCount.m[tmp.idx,c] <- PCGHubCount.v

	threshold <- quantile(PCGHubCount.v, probs = c(95)/100)
	PCGHubMarker.v <- ifelse(PCGHubCount.v > threshold, "Hub", "Non Hub")

	#PCGHubs.v <- PCGHubMarker.v[which(PCGHubMarker.v == 1)]
	tmp.idx <- match(names(PCGHubMarker.v), rownames(PCGHubMarker.m))
	PCGHubMarker.m[tmp.idx,c] <- PCGHubMarker.v 

	PCGExp.v <-  rowMeans(PCGexpFPKM.m)
	tmp.idx <- match(names(PCGHubMarker.v),  names(PCGExp.v))

	data.df <- as.data.frame(cbind(PCGExp.v[tmp.idx], PCGHubMarker.v ))
	colnames(data.df) <- c("Exp", "Type")
	data.df$Exp <- as.numeric(data.df$Exp)
	data.df <- within(data.df, Type <- factor(Type, levels = c("Hub", "Non Hub")))
	data.df$Cancer <- cancers.v[c]
	dataEP.df <- rbind(dataEP.df, data.df)

}
library(ggpubr)

p <- ggboxplot(dataEL.df, x="Cancer", y="Exp", fill = "Cancer", color="Type", palette = mycolors.v, add = "none", outlier.size=0.2,  ylab="Exp",  xlab="Cancer Types", ggtheme = theme_bw()) + theme(plot.title = element_text(hjust = 0.5),  axis.text.x=element_blank())  
p <- p+ stat_compare_means(aes(group=Type),label = "p.signif")
p <- p+ scale_colour_manual(values=c("darkred", "black"))
p <- p+ annotate("text", x=15, y=7.5, label="****  P < 0.0001") 
#p <- p+ annotate("text", x=15, y=7, label="***  P < 0.001") 
p <- p+ annotate("text", x=15, y=7, label="**  P < 0.01") 
p
p1 <- p
ggsave("LncExp_HubCompare_by_PanCancer1.pdf", width=15,height=6)

p <- ggboxplot(dataEP.df, x="Cancer", y="Exp", fill = "Cancer", color="Type", palette = mycolors.v, add = "none", outlier.size=0.2, ylab="Exp", xlab="Cancer Types", ggtheme = theme_bw()) + theme(plot.title = element_text(hjust = 0.5),  axis.text.x=element_blank())
p <- p+ stat_compare_means(aes(group=Type),label = "p.signif")
p <- p+ scale_colour_manual(values=c("darkred", "black"))
#p <- p+ annotate("text", x=3, y=600, label="****  P < 0.0001") 
p

p2 <- p
ggsave("PCGExp_HubCompare_by_PanCancer.pdf", width=12,height=6)


ggarrange(p1, p2, nrow=2, ncol=1, legend="right", common.legend=TRUE)
ggsave("LNC_PCGExp_HubCompare_by_PanCancer.pdf", width=12,height=6)


############################

dataEL.df <- data.frame(Exp=0, Type=0, Cancer=0)
dataEL.df <- dataE.df[-1,]

dataEP.df <- data.frame(Exp=0, Type=0, Cancer=0)
dataEP.df <- dataEP.df[-1,]


LncHubCount.m <- matrix(0, nrow=nrow(LNCexpFPKM.m), ncol=18)
rownames(LncHubCount.m) <- rownames(LNCexpFPKM.m)
colnames(LncHubCount.m) <- cancers.v

PCGHubCount.m <- matrix(0, nrow=nrow(PCGexpFPKM.m), ncol=18)
rownames(PCGHubCount.m) <- rownames(PCGexpFPKM.m)
colnames(PCGHubCount.m) <- cancers.v

LncHubMarker.m <- matrix(0, nrow=nrow(LNCexpFPKM.m), ncol=18)
rownames(LncHubMarker.m) <- rownames(LNCexpFPKM.m)
colnames(LncHubMarker.m) <- cancers.v

PCGHubMarker.m <- matrix(0, nrow=nrow(PCGexpFPKM.m), ncol=18)
rownames(PCGHubMarker.m) <- rownames(PCGexpFPKM.m)
colnames(PCGHubMarker.m) <- cancers.v

for(c in 1:length(cancers.v)) {

	print(cancers.v[c])
	load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
	load(paste("./InteractNetworkStat1/", cancers.v[c], "_InteractNetStat.Rdata", sep="") )
	load(paste("../CancerRdata/exp_FPKM/", cancers.v[c], "_LNC_expFPKM.Rdata", sep=""))
	load(paste("../CancerRdata/Exp_gene/", cancers.v[c], "_expGenelevel.Rdata", sep=""))

	LncHubCount.v <- rowSums(InteractLG.m) 
	tmp.idx <- match(names(LncHubCount.v), rownames(LncHubCount.m))
	LncHubCount.m[tmp.idx,c] <- LncHubCount.v 

	threshold <- quantile(LncHubCount.v, probs = c(95)/100)
	LncHubMarker.v <- ifelse(LncHubCount.v > threshold, 1, 0)

	#LncHubs.v <- LncHubMarker.v[which(LncHubMarker.v == 1)]
	tmp.idx <- match(names(LncHubMarker.v), rownames(LncHubMarker.m))
	LncHubMarker.m[tmp.idx,c] <- LncHubMarker.v

	LncExp.v <-  rowMeans(LNCexpFPKM.m)
	tmp.idx <- match(names(LncHubMarker.v),  names(LncExp.v))

	data.df <- as.data.frame(cbind(LncExp.v[tmp.idx], LncHubMarker.v ))
	colnames(data.df) <- c("Exp", "Type")
	data.df$Exp <- as.numeric(data.df$Exp)
	data.df <- within(data.df, Type <- factor(Type, levels = c(1, 0)))
	data.df$Cancer <- cancers.v[c]
	dataEL.df <- rbind(dataEL.df, data.df)

        ######
	PCGHubCount.v <- colSums(InteractLG.m) 
	tmp.idx <- match(names(PCGHubCount.v), rownames(PCGHubCount.m))
	PCGHubCount.m[tmp.idx,c] <- PCGHubCount.v

	threshold <- quantile(PCGHubCount.v, probs = c(95)/100)
	PCGHubMarker.v <- ifelse(PCGHubCount.v > threshold, 1, 0)

	#PCGHubs.v <- PCGHubMarker.v[which(PCGHubMarker.v == 1)]
	tmp.idx <- match(names(PCGHubMarker.v), rownames(PCGHubMarker.m))
	PCGHubMarker.m[tmp.idx,c] <- PCGHubMarker.v 

	PCGExp.v <-  rowMeans(PCGexpFPKM.m)
	tmp.idx <- match(names(PCGHubMarker.v),  names(PCGExp.v))

	data.df <- as.data.frame(cbind(PCGExp.v[tmp.idx], PCGHubMarker.v ))
	colnames(data.df) <- c("Exp", "Type")
	data.df$Exp <- as.numeric(data.df$Exp)
	data.df <- within(data.df, Type <- factor(Type, levels = c(1, 0)))
	data.df$Cancer <- cancers.v[c]
	dataEP.df <- rbind(dataEP.df, data.df)

}

sort(rowSums(LncHubMarker.m))


##############
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
keytypes(org.Hs.eg.db)

########ENSG00000241684
tmp.idx <-  which(rownames(LncHubMarker.m) == "ENSG00000241684.4")
LncHubMarker.m[tmp.idx,]
validCancer.idx <- which(LncHubMarker.m[tmp.idx,] == 1)
validCancer.idx 


c <- 1 ## blca
sample_id <- names(validCancer.idx)
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )

tmp.idx <- which(rownames(InteractLG.m) == "ENSG00000241684.4")
#InteractLG.m[tmp.idx, ]
hubtarget.idx <- which(InteractLG.m[tmp.idx,] == 1)
hubtarget.v <- names(hubtarget.idx)

geneID.v <- hubtarget.v
core_Gene_ID.v <- geneID.v
gene_num.v <- length(geneID.v)

allHub.v <- c()
for(c in as.vector(validCancer.idx)[-1] ) { ### exclude blca

	print(cancers.v[c])
	load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )

	tmp.idx <- which(rownames(InteractLG.m) == "ENSG00000241684.4")
	#InteractLG.m[tmp.idx, ]
	hubtarget.idx <- which(InteractLG.m[tmp.idx,] == 1)
	hubtarget.v <- names(hubtarget.idx)
	allHub.v <- c(allHub.v, hubtarget.v)

	geneID.v <- hubtarget.v
	core_Gene_ID.v <- intersect(core_Gene_ID.v, geneID.v )
	gene_num.v <- c(gene_num.v, length(geneID.v ))
	##ego <- enrichGO(hubtarget.v, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.01, qvalueCutoff = 0.05, keyType = 'ENTREZID')
	#as.data.frame(ego)

}

sort(table(allHub.v))
core_num <- length(core_Gene_ID.v)

library(plotrix)

ellipse_col <- c('#6181BD4E','#F348004E','#64A10E4E','#9300264E','#464E044E','#049a0b4E','#4E0C664E','#D000004E','#FF6C004E','#FF00FF4E','#c7475b4E','#00F5FF4E','#BDA5004E','#A5CFED4E','#f0301c4E','#2B8BC34E','#FDA1004E','#54adf54E','#CDD7E24E','#9295C14E')

ellipse_col <-  mycolors.v[as.vector(validCancer.idx)]

flower_plot <- function(sample, gene_num, core_gene, start, a, b, r, ellipse_col, circle_col) {
    par( bty = 'n', ann = F, xaxt = 'n', yaxt = 'n', mar = c(1,1,1,1))
    plot(c(0,10),c(0,10),type='n')
    n   <- length(sample)
    deg <- 360 / n
    res <- lapply(1:n, function(t){
        draw.ellipse(x = 5 + cos((start + deg * (t - 1)) * pi / 180), 
            y = 5 + sin((start + deg * (t - 1)) * pi / 180), 
            col = ellipse_col[t],
            border = ellipse_col[t],
            a = a, b = b, angle = deg * (t - 1))
        text(x = 5 + 2.5 * cos((start + deg * (t - 1)) * pi / 180),
            y = 5 + 2.5 * sin((start + deg * (t - 1)) * pi / 180),
            gene_num[t])
        
        if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
            text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
                y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
                sample[t],
                srt = deg * (t - 1) - start,
                adj = 1,
                cex = 1
                )
        } else {
            text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
                y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
                sample[t],
                srt = deg * (t - 1) + start,
                adj = 0,
                cex = 1
                )
        }            
    })
    draw.circle(x = 5, y = 5, r = r, col = circle_col, border = NA)
    text(x = 5, y = 5, paste('Core:', core_gene))
}

pdf('Hub_ENSG00000241684_target_flower1.pdf', width = 6, height = 6)
flower_plot(sample = sample_id, gene_num = gene_num.v, core_gene = core_num, 
    start = 90, a = 0.5, b = 2, r = 1, ellipse_col = ellipse_col, circle_col = 'white')
dev.off()
### core gene:222698

> validCancer.idx 
     BLCA BRCA(ER+) BRCA(ER-)      CESC      COAD      ESCA      HNSC      PRAD 
        1         2         3         4         6         7         8        15 
     READ      UCEC 
       16        18 

2,  8,  16,  18
library(stringr)

c <- 1
print(cancers.v[c])
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
tmp.idx <- which(rownames(InteractLG.m) == "ENSG00000241684.4")
#InteractLG.m[tmp.idx, ]
hubtarget.idx <- which(InteractLG.m[tmp.idx,] == 1)
hubtarget.v <- names(hubtarget.idx)
ego <- enrichGO(hubtarget.v, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.01, qvalueCutoff = 0.05, keyType = 'ENTREZID')
as.data.frame(ego)

#ego2 <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
ego3 <- setReadable(ego, OrgDb = org.Hs.eg.db)
GO <- ego3[, c(1,2,8,6)]
#GO <- ego3[c(2,3,6,7,8,9), c(1,2,8,6)]
GO$geneID <- str_replace_all(GO$geneID, "/", ",")
names(GO) <- c("ID", "Term", "Genes", "adj_pval")
GO$Category = "BP"

load(paste("./MultiOmicsData/", cancers.v[c], "_MultiOmicsProcess.Rdata", sep=""))
tmp.idx <- match(hubtarget.v, rownames(DEG.m))
tmp.m <- DEG.m[tmp.idx,]
hubtargetName.v <-  unlist(mget(hubtarget.v, org.Hs.egSYMBOL, ifnotfound=NA))
genedata <- data.frame(ID= hubtargetName.v, logFC=tmp.m[,1])

circ <- circle_dat(GO, genedata)
#GOCircle(circ,table.legend = F,label.size=5,nsub=nrow(GO))
chord <- chord_dat(data=circ, genes= genedata, process=GO$Term)

pdf("HubTarget_GO_Chord_ENSG00000241684_BLCA.pdf", width=19, height=20)
GOChord(chord, space=0.01, gene.order='logFC', gene.space=0.25, gene.size=12, ribbon.col=brewer.pal(length(GO$Term),'Dark2'), process.label=16)
dev.off()


pdf("HubTarget_GO_Chord_ENSG00000241684_COAD.pdf", width=19, height=20)
GOChord(chord, space=0.01, gene.order='logFC', gene.space=0.25, gene.size=8, ribbon.col=c(brewer.pal(length(GO$Term),'Dark2'), brewer.pal(4,'Set2')), process.label=14)
dev.off()




