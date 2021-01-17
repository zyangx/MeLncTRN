options(stringsAsFactors=FALSE)

load("../CancerRdata/exp_FPKM/BLCA_LNC_expFPKM.Rdata")
load("../CancerRdata/Exp_gene/BLCA_expGenelevel.Rdata")


cancers.v  <-  c("BLCA", "BRCA(ER+)", "BRCA(ER-)", "CESC", "CHOL", "COAD", "ESCA", "HNSC",
				 "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PAAD"  , "PRAD", "READ",  "THCA", "UCEC")

LncPCGInterAct.v <- c()
for(i in 1:length(rownames(LNCexpFPKM.m)) ){
	tmp.v <- paste(rownames(LNCexpFPKM.m)[i] , rownames(PCGexpFPKM.m), sep="_")
	LncPCGInterAct.v <- c(LncPCGInterAct.v , tmp.v)
}



LncPCGInterAct.v <- c()
for(c in 1:length(cancers.v)) {
	print(cancers.v[c])
	#load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
	load(paste("./InteractNetworkStat1/", cancers.v[c], "_InteractNetStat.Rdata", sep="") )
	tmp.v <- paste(validLncTargetStat.df[,1], validLncTargetStat.df[,2], sep="_")
        LncPCGInterAct.v <- c(LncPCGInterAct.v, tmp.v)
}
LncPCGInterAct.v <- unique(LncPCGInterAct.v)


LncCancerCount.m <- matrix(0, nrow=nrow(LNCexpFPKM.m), ncol=18)
rownames(LncCancerCount.m) <- rownames(LNCexpFPKM.m)
colnames(LncCancerCount.m) <- cancers.v

PCGCancerCount.m <- matrix(0, nrow=nrow(PCGexpFPKM.m), ncol=18)
rownames(PCGCancerCount.m) <- rownames(PCGexpFPKM.m)
colnames(PCGCancerCount.m) <- cancers.v

LncPCGActCancerCount.m <- matrix(0, nrow=length(LncPCGInterAct.v), ncol=18)
rownames(LncPCGActCancerCount.m) <- LncPCGInterAct.v
colnames(LncPCGActCancerCount.m) <- cancers.v


for(c in 1:length(cancers.v)) {

	print(cancers.v[c])
	load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
	load(paste("./InteractNetworkStat1/", cancers.v[c], "_InteractNetStat.Rdata", sep="") )

	tmp.idx <- match(rownames(InteractLG.m), rownames(LncCancerCount.m))
	LncCancerCount.m[tmp.idx,c] <- 1

	tmp.idx <- match(colnames(InteractLG.m), rownames(PCGCancerCount.m))
	PCGCancerCount.m[tmp.idx,c] <- 1

	tmp.idx <- match(paste(validLncTargetStat.df[,1], validLncTargetStat.df[,2], sep="_"), rownames(LncPCGActCancerCount.m))
	LncPCGActCancerCount.m[tmp.idx,c] <- 1

}

#################
LncCancerCount.v <- rowSums(LncCancerCount.m)
tmp.idx <- which(LncCancerCount.v > 0)
LncCancerCount.v <- LncCancerCount.v[tmp.idx]
LncCancerCount.m <- LncCancerCount.m[tmp.idx,]

popp <- data.frame(table(LncCancerCount.v))
popp[,2] <- popp[,2] / sum(popp[,2] )
colnames(popp)[1] <- "Num"
######
PCGCancerCount.v <- rowSums(PCGCancerCount.m)
tmp.idx <- which(PCGCancerCount.v > 0)
PCGCancerCount.v <- PCGCancerCount.v[tmp.idx]
PCGCancerCount.m <- PCGCancerCount.m[tmp.idx,]

popp1 <- data.frame(table(PCGCancerCount.v))
popp1[,2] <- popp1[,2] / sum(popp1[,2] )
colnames(popp1)[1] <- "Num"
######
LncPCGActCancerCount.v <- rowSums(LncPCGActCancerCount.m)
tmp.idx <- which(LncPCGActCancerCount.v > 0)
LncPCGActCancerCount.v <- LncPCGActCancerCount.v[tmp.idx]
LncPCGActCancerCount.m <- LncPCGActCancerCount.m[tmp.idx,]

popp2 <- data.frame(table(LncPCGActCancerCount.v))
popp2[,2] <- popp2[,2] / sum(popp2[,2] )
colnames(popp2)[1] <- "Num"
#################

prop.m <- matrix(NA, nrow=18, ncol=4)
colnames(prop.m) <- c("Num", "Lnc","Gene","Pair")
prop.m[,1:2] <-  as.matrix(popp)
prop.m[1:15,3] <-  popp1[,2]
prop.m[1:12,4] <-  popp2[,2]
prop.df <- as.data.frame(prop.m)
##
prop.df <- rbind(rbind(popp, popp1), popp2)
prop.df[,2] <- prop.df[,2] *100
prop.df$type <- c(rep("Lnc",18), rep("Gene",15), rep("Pair",12))
colnames(prop.df)  <- c("Number_of_Cancer_Type", "Percent", "Type")

#plot(as.numeric(levels(popp[,1])[popp[,1]] ), popp[,2], main = cancers.v[c], xlab="Degree of Genes",ylab="Frequency", type = "p", col = "blue", pch=18,  cex=1.2, cex.main=1.2, cex.lab=1.1)


#plot(as.numeric(prop.df[,1]),prop.df[,2],xlab="Num of Cancer Type", ylab="Frequency",col = "black", lwd=2,main = "")

library(ggplot2)
library(reshape2)

ggplot(prop.df, aes(x=Number_of_Cancer_Type, y=Percent, group=Type)) + geom_line(aes(colour=Type)) + geom_point(size=1, shape=20) + theme_bw() + scale_x_discrete(breaks=c(0,3,6,9,12,15,18), labels=c(0,3,6,9,12,15,18)) +  theme(legend.position = c(0.86,0.79))
ggsave("Cancer_lnc_gene_inter_percent.pdf", width=3.5, height=3.2)


library(RColorBrewer)
library(extrafont)
library(ggsci)
library(scales)
mycolors.v <- pal_d3("category20")(18)

###### Lnc
tmp.idx <- which(rowSums(LncCancerCount.m) == 1)
uniqCount.m <- LncCancerCount.m[tmp.idx,]
LncStat.df <- as.data.frame(colSums(uniqCount.m))
colnames(LncStat.df) <- "Num"

slicesNum <- as.numeric(LncStat.df[,1])
names(slicesNum) <- cancers.v;
lbls <- paste(names(slicesNum), "\n(", slicesNum, ")", sep= "");

pdf("Cancer_Specific_Lnc_Count.pdf", width=3, height=3)
par(family="Arial")
par(mar=c(0.28,0.2,0.1,0.52))
pie(slicesNum, labels=lbls, col=mycolors.v, clockwise=TRUE, radius = 0.8, cex= 0.38, font = 2  )
text(-0.85, 0.8, labels = "A", xpd = TRUE, cex = 0.95, font=2)
#text(0,1, labels = "Exp", xpd = TRUE, cex = 1.15, font=2)
dev.off()


###### Gene
tmp.idx <- which(rowSums(PCGCancerCount.m) == 1)
uniqCount.m <- PCGCancerCount.m[tmp.idx,]
GeneStat.df <- as.data.frame(colSums(uniqCount.m))
colnames(GeneStat.df) <- "Num"

slicesNum <- as.numeric(GeneStat.df[,1])
names(slicesNum) <- cancers.v;
lbls <- paste(names(slicesNum), "\n(", slicesNum, ")", sep= "");

pdf("Cancer_Specific_Gene_Count.pdf", width=3, height=3)
par(family="Arial")
par(mar=c(0.28,0.2,0.1,0.52))
pie(slicesNum, labels=lbls, col=mycolors.v, clockwise=TRUE, radius = 0.8, cex= 0.38, font = 2  )
text(-0.85, 0.8, labels = "B", xpd = TRUE, cex = 0.95, font=2)
#text(0,1, labels = "Exp", xpd = TRUE, cex = 1.15, font=2)
dev.off()


###### Lnc_Gene_interact
tmp.idx <- which(rowSums(LncPCGActCancerCount.m) == 1)
uniqCount.m <- LncPCGActCancerCount.m[tmp.idx,]
PairStat.df <- as.data.frame(colSums(uniqCount.m))
colnames(PairStat.df) <- "Num"

slicesNum <- as.numeric(PairStat.df[,1])
names(slicesNum) <- cancers.v;
lbls <- paste(names(slicesNum), "\n(", slicesNum, ")", sep= "");

pdf("Cancer_Specific_LncGenePair_Count.pdf", width=3, height=3)
par(family="Arial")
par(mar=c(0.28,0.2,0.1,0.52))
pie(slicesNum, labels=lbls, col=mycolors.v, clockwise=TRUE, radius = 0.8, cex= 0.38, font = 2  )
text(-0.85, 0.8, labels = "C", xpd = TRUE, cex = 0.95, font=2)
#text(0,1, labels = "Exp", xpd = TRUE, cex = 1.15, font=2)
dev.off()


#################
#################
LncPanCancerCount.df <- as.data.frame(LncCancerCount.v)
LncPanCancerCount.df$type <- "Moderate"
LncPanCancerCount.df$type[which(LncPanCancerCount.df[,1] > 14)] <- "Pan-cancer"
LncPanCancerCount.df$type[which(LncPanCancerCount.df[,1]  ==  1)] <- "Cancer-specific"
LncPanCancerCount.df$class <- "Others"
LncPanCancerCount.df$class[which(LncPanCancerCount.df[,1] > 14)] <- "Pan-cancer"

table(LncPanCancerCount.df$class)

library(ggpubr)

pcom.lv <- list()
for(c in 1:length(cancers.v)) {

print(cancers.v[c])
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )

LncDegree.v <- rowSums(InteractLG.m)
tmp.idx <- match(names(LncDegree.v),  rownames(LncPanCancerCount.df))

data.df <- as.data.frame(cbind(LncDegree.v, LncPanCancerCount.df[tmp.idx,2]))
colnames(data.df) <- c("Degree", "Type")
data.df$Degree <- as.numeric(data.df$Degree)
data.df <- within(data.df, Type <- factor(Type, levels = c("Cancer-specific", "Moderate", "Pan-cancer")))

p <- ggboxplot(data.df, x = "Type", y = "Degree",  fill = "Type", width = 0.8, palette = brewer.pal(3,"Blues"),  add = "none", ylim=c(0, max(data.df$Degree)*1.4), title=cancers.v[c], ylab="Gene Number",xlab="", ggtheme = theme_bw()) +theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_blank())
p <- p+stat_compare_means(comparisons = list(c("Pan-cancer", "Moderate"), c("Moderate", "Cancer-specific"), c("Pan-cancer", "Cancer-specific")),label= "p.format")

pcom.lv[[c]] <- p

}

ggarrange(pcom.lv[[1]], pcom.lv[[2]], pcom.lv[[3]], pcom.lv[[4]], pcom.lv[[5]], pcom.lv[[6]], pcom.lv[[7]], pcom.lv[[8]], pcom.lv[[9]], pcom.lv[[10]], pcom.lv[[11]], pcom.lv[[12]], pcom.lv[[13]], pcom.lv[[14]], pcom.lv[[15]], pcom.lv[[16]], pcom.lv[[17]], pcom.lv[[18]], nrow=5, ncol=4, legend="bottom", common.legend=TRUE)

ggsave("LncCategory_Compare_by_Cancer_Distribution.pdf", width=10,height=13);


dataC.df <- data.frame(Degree=0, Type=0, Cancer=0)
dataC.df <- dataC.df[-1,]

for(c in 1:length(cancers.v)) {

print(cancers.v[c])
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )

LncDegree.v <- rowSums(InteractLG.m)
tmp.idx <- match(names(LncDegree.v),  rownames(LncPanCancerCount.df))

data.df <- as.data.frame(cbind(LncDegree.v, LncPanCancerCount.df[tmp.idx,3]))
colnames(data.df) <- c("Degree", "Type")
data.df$Degree <- as.numeric(data.df$Degree)
data.df <- within(data.df, Type <- factor(Type, levels = c("Pan-cancer", "Others")))
data.df$Cancer <- cancers.v[c]

dataC.df <- rbind(dataC.df, data.df)

}


p <- ggboxplot(dataC.df, x="Cancer", y="Degree", fill = "Cancer", color="Type", palette = mycolors.v, add = "none", size=0.2, ylab="Number of genes regulated", xlab="Cancer Types", ggtheme = theme_bw()) + theme(plot.title = element_text(hjust = 0.5),  axis.text.x=element_blank())
p <- p+ stat_compare_means(aes(group=Type),label = "p.signif")
p <- p+ scale_colour_manual(values=c("darkred", "black"))
p <- p+ annotate("text", x=3, y=600, label="****  P < 0.0001") 
p
ggsave("LncCategory_Compare_by_PanCancer_Distribution.pdf", width=12,height=6)


#dataC.df <- within(dataC.df, Type <- factor(Type, levels = c("Others", "Pan-cancer")))
dataC.df <- within(dataC.df, Type <- factor(Type, levels = c("Pan-cancer", "Others")))
p <- ggboxplot(dataC.df, x="Cancer", y="Degree", fill = "Cancer", color="Type", palette = mycolors.v, add = "none", size=0.2, ylab="Number of genes regulated", xlab="Cancer Types", ggtheme = theme_bw()) + theme(plot.title = element_text(hjust = 0.5),  axis.text.y=element_blank() )
p <- p+ stat_compare_means(aes(group=Type),label = "p.signif")
p <- p+ scale_colour_manual(values=c("darkred", "black"))
p <- p+ coord_flip()
p <- p+ annotate("text", x=15, y=550, label="****  P < 0.0001") 
p <- p+ scale_x_discrete(limits=rev(cancers.v))
p
ggsave("LncCategory_Compare_by_PanCancer_Distribution1.pdf", width=6,height=10)



#######################
#######################
load("./GeneConversation/GenecodePhastCons_100way.Rdata")
tmp.idx <- match(rownames(LncPanCancerCount.df), rownames(GenecodePhastCons.df))
LncPanCancerCount.df$Conservation <- GenecodePhastCons.df$Conservation[tmp.idx]

LncPanCancerCount.df <- within(LncPanCancerCount.df, type <- factor(type, levels = c("Cancer-specific", "Moderate","Pan-cancer")))
with(LncPanCancerCount.df, levels(type))

p <- ggboxplot(LncPanCancerCount.df, x = "type", y = "Conservation",  fill = "type", width = 0.8, palette = brewer.pal(3,"Blues"),  add = "none", title="", ylab="Conservation",xlab="", ggtheme = theme_bw()) +theme(plot.title = element_text(hjust = 0.5))

p <- p+stat_compare_means(comparisons = list(c("Pan-cancer", "Moderate"), c("Moderate", "Cancer-specific"), c("Pan-cancer", "Cancer-specific")),label= "p.format")
p <- p + guides(fill=FALSE)
p
ggsave("LncCategory_Conservation.pdf", width=4,height=5)

####
load("./GeneConversation/GenecodePhastCons_100wayPromoter.Rdata")
tmp.idx <- match(rownames(LncPanCancerCount.df), rownames(GenecodePhastCons.df))
LncPanCancerCount.df$Conservation <- NA
LncPanCancerCount.df$Conservation <- GenecodePhastCons.df$Conservation[tmp.idx]

p <- ggboxplot(LncPanCancerCount.df, x = "type", y = "Conservation",  fill = "type", width = 0.8, palette = brewer.pal(3,"Blues"),  add = "none", title="", ylab="Conservation",xlab="", ggtheme = theme_bw()) +theme(plot.title = element_text(hjust = 0.5))

p <- p+stat_compare_means(comparisons = list(c("Pan-cancer", "Moderate"), c("Moderate", "Cancer-specific"), c("Pan-cancer", "Cancer-specific")),label= "p.format")
p <- p + guides(fill=FALSE)
p
ggsave("LncCategory_Conservation_promoter.pdf", width=4,height=5)


pdf("LncCategory_Conservation_scatter1.pdf", width=5,height=5)
par(mar=c(4,4,3,1));
#smoothScatter(as.numeric(as.vector(LncPanCancerCount.df[,1])), LncPanCancerCount.df[,4], colramp = Lab.palette, nrpoints=0, main="", xlab="Number of Cancer", ylab="Conservation", cex.lab=1.2, las=1)
smoothScatter(as.numeric(as.vector(LncPanCancerCount.df[,1])), LncPanCancerCount.df[,4],  nrpoints=0, main="", xlab="Number of Cancer", ylab="Conservation", cex.lab=1.2, las=1)
lm.o <- lm(as.numeric(LncPanCancerCount.df[,4]) ~ as.numeric(as.vector(LncPanCancerCount.df[,1])));
rsquare=parse(text=paste('R^2==',round(summary(lm.o)$r.sq,2),sep=''));
cor.o <- cor.test(as.numeric(LncPanCancerCount.df[,4]) , as.numeric(as.vector(LncPanCancerCount.df[,1])))
#text(x=max(data.m[,1])*0.65, y=max(data.m[,2])*0.95, adj= 0, cex=0.9, font=2, col="#FF0000", rsquare);
text(x=13, y=0.9, adj= 0, cex=0.9, font=2, col="#FF0000", label= paste('PCC = ',round(cor.o$estimate,digits=2),sep=''));
text(x=13, y=0.8, adj= 0, cex=0.9, font=2, col="#FF0000", label= paste('P = ',format(cor.o$p.value,digits=3, scientific=T),sep=''));

dev.off()


#######################
#######################
load("./TissueSpecifity/TissueSpecifityTau.Rdata")
tmp.idx <- match(rownames(LncPanCancerCount.df), names(TissueSpecifity.v))
LncPanCancerCount.df$TissueSpecifity <- as.numeric(TissueSpecifity.v[tmp.idx])

p <- ggboxplot(LncPanCancerCount.df, x = "type", y = "TissueSpecifity",  fill = "type", width = 0.8, palette = brewer.pal(3,"Blues"),  add = "none", title="", ylab="Tau", xlab="", ggtheme = theme_bw()) +theme(plot.title = element_text(hjust = 0.5))
p <- p+stat_compare_means(comparisons = list(c("Pan-cancer", "Moderate"), c("Moderate", "Cancer-specific"), c("Pan-cancer", "Cancer-specific")),label= "p.format")
p <- p + guides(fill=FALSE)
p
ggsave("LncCategory_TissueSpecifity.pdf", width=4,height=5)

####
load("./TissueSpecifity/TissueSpecifityTauGTEx.Rdata")
#tmp.idx <- match( strsplit(rownames(LncPanCancerCount.df),"\.",fixed= T) ,  names(TissueSpecifity.v)  )
tmp.idx <- match( unlist(lapply( strsplit(rownames(LncPanCancerCount.df),".",fixed= T) , function(x) x[1])), unlist(lapply( strsplit(names(TissueSpecifity.v),".",fixed= T) , function(x) x[1])) )

LncPanCancerCount.df$TissueSpecifity <- NA
LncPanCancerCount.df$TissueSpecifity[which(!is.na(tmp.idx))] <- TissueSpecifity.v[setdiff(tmp.idx, NA)]

p <- ggboxplot(LncPanCancerCount.df, x = "type", y = "TissueSpecifity",  fill = "type", width = 0.8, palette = brewer.pal(3,"Blues"),  add = "none", title="", ylab="Tau", xlab="", ggtheme = theme_bw()) +theme(plot.title = element_text(hjust = 0.5))
p <- p+stat_compare_means(comparisons = list(c("Pan-cancer", "Moderate"), c("Moderate", "Cancer-specific"), c("Pan-cancer", "Cancer-specific")),label= "p.format")
p <- p + guides(fill=FALSE)
p
ggsave("LncCategory_TissueSpecifity_by_GTEx.pdf", width=4,height=5)



