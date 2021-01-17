options(stringsAsFactors=FALSE)

cancers.v  <-  c("BLCA", "BRCA(ER+)", "BRCA(ER-)", "CESC", "CHOL", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC",
				"LUAD", "LUSC", "PAAD" , "PRAD", "READ", "THCA", "UCEC")


DELstat.m <- matrix(NA, nrow=18, ncol=6)
rownames(DELstat.m) <- cancers.v
colnames(DELstat.m) <- c("UpPercent", "DownPercent", "OnPercent", "OffPercent", "UpOnPercent", "DownOffPercent")

DELPvalue.m <- matrix(NA, nrow=18, ncol=3)
rownames(DELPvalue.m) <- cancers.v
colnames(DELPvalue.m) <- c("Up/On", "Down/Off", "All")

load("~/backup/TCGA-XenaBrowser/gencode/gencode.Rdata")


for(c in 1:length(cancers.v)) {
#rm(list=ls())
print(cancers.v[c])
load(paste("./MultiOmicsData/", cancers.v[c], "_MultiOmicsProcess.Rdata", sep="") )
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
load(paste("./LncDEstatus/", cancers.v[c], "_LncDEStatus.Rdata", sep="") )

tmp.idx <-  intersect(which(DEL.m[,1] > 0),  which(DEL.m[,5] < 0.05))
upLnc.v <-  rownames(DEL.m)[tmp.idx]
DELstat.m[c,1] <-  length(intersect(upLnc.v, rownames(InteractLG.m) ) ) / nrow(InteractLG.m)

tmp.idx <-  intersect(which(DEL.m[,1] < 0),  which(DEL.m[,5] < 0.05))
downLnc.v <-  rownames(DEL.m)[tmp.idx]
DELstat.m[c,2] <-  length(intersect(downLnc.v, rownames(InteractLG.m))) / nrow(InteractLG.m)

tmp.idx <-  intersect(which(LncDEStatus.m[,9] > 2), which(LncDEStatus.m[,11] < 0.05))
onLnc.v <- rownames(LncDEStatus.m)[tmp.idx]
DELstat.m[c,3] <-  length(intersect(onLnc.v, rownames(InteractLG.m))) / nrow(InteractLG.m)

tmp.idx <-  intersect(which(LncDEStatus.m[,9] < 0.5), which(LncDEStatus.m[,11] < 0.05))
offLnc.v <- rownames(LncDEStatus.m)[tmp.idx]
DELstat.m[c,4] <-  length(intersect(offLnc.v, rownames(InteractLG.m))) / nrow(InteractLG.m)

uponLnc.v <- union(upLnc.v, onLnc.v)
DELstat.m[c,5] <-  length(intersect(uponLnc.v, rownames(InteractLG.m))) / nrow(InteractLG.m)

downoffLnc.v <- union(downLnc.v, offLnc.v)
DELstat.m[c,6] <-  length(intersect(downoffLnc.v, rownames(InteractLG.m))) / nrow(InteractLG.m)


#########
backgroundInter.v <- intersect(uponLnc.v, names(statLMG.lm))
DELinter.v <- intersect(uponLnc.v, rownames(InteractLG.m))
data.m <- matrix(NA, nrow=2, ncol=2)
data.m[1,1] <- length(names(statLMG.lm))
data.m[2,1] <- length(backgroundInter.v)
data.m[1,2] <- nrow(InteractLG.m)
data.m[2,2] <- length(DELinter.v)
print(fisher.test(data.m)$p.value)
DELPvalue.m[c, 1] <- fisher.test(data.m)$p.value

#########
backgroundInter.v <- intersect(downoffLnc.v, names(statLMG.lm))
DELinter.v <- intersect(downoffLnc.v, rownames(InteractLG.m))
data.m <- matrix(NA, nrow=2, ncol=2)
data.m[1,1] <- length(names(statLMG.lm))
data.m[2,1] <- length(backgroundInter.v)
data.m[1,2] <- nrow(InteractLG.m)
data.m[2,2] <- length(DELinter.v)
print(fisher.test(data.m)$p.value)
DELPvalue.m[c, 2] <- fisher.test(data.m)$p.value

#########
allDiffLnc.v <- c(uponLnc.v, downoffLnc.v)
backgroundInter.v <- intersect(allDiffLnc.v, names(statLMG.lm))
DELinter.v <- intersect(allDiffLnc.v, rownames(InteractLG.m))
data.m <- matrix(NA, nrow=2, ncol=2)
data.m[1,1] <- length(names(statLMG.lm))
data.m[2,1] <- length(backgroundInter.v)
data.m[1,2] <- nrow(InteractLG.m)
data.m[2,2] <- length(DELinter.v)
print(fisher.test(data.m)$p.value)
DELPvalue.m[c, 3] <- fisher.test(data.m)$p.value

}

############################

library("reshape2")
library("ggpubr")
library(RColorBrewer)
library(extrafont)
font_import()

tmp.df <- as.data.frame(melt(DELstat.m[,c(5,6)]))
colnames(tmp.df) <- c("Cancer", "Type", "Value")
p <- ggbarplot(tmp.df, x = "Cancer", y = "Value", fill= "Type", add = "none", title="", palette = c("#FFED6F", "#377EB8"), color="white", ylab="Proportion", xlab="", ggtheme = theme_bw()) 
p <- p + scale_fill_manual( values=c("#EFC000FF", "#0073C2FF"), breaks=c("UpOnPercent","DownOffPercent"), labels=c("Up/On", "Down/Off")) + theme(legend.position=c(0.89,0.91),legend.text=element_text(family="Arial"), axis.text.x = element_text(angle = 45, vjust = 1,  size = 8, hjust = 1)) + guides(fill=guide_legend(title=NULL))
p
ggsave("Lnc_Differential_Exp_Stat.pdf",width=7.5,height=5.5)


p <- ggbarplot(tmp.df, x = "Cancer", y = "Value", fill= "Type", add = "none", title="", palette = c("#FFED6F", "#377EB8"), color="white", ylab="Proportion", xlab="", ggtheme = theme_bw()) 
p <- p + scale_fill_manual( values = c("#EFC000FF", "#0073C2FF"), breaks=c("UpOnPercent","DownOffPercent"), labels=c("Up/On", "Down/Off")) + theme(legend.position="top", legend.text=element_text(family="Arial", size = 6), axis.text.x = element_text( vjust = 1,  size = 6, hjust = 1), axis.text.y = element_text( vjust = 1,  size = 6, hjust = 1)) + guides(fill=guide_legend(title=NULL))
p <- p+ coord_flip()  + scale_x_discrete(limits=rev(cancers.v))
p
ggsave("Lnc_Differential_Exp_Stat1.pdf",width=3.5,height=5)


DELstat.m <- round(DELstat.m, digits=2)

pcom.lv <- list()
for(c in 1:nrow(DELstat.m)) {
tmp.m <- matrix( as.vector(DELstat.m[c,c(1,2,3,4)]), nrow=2) 
rownames(tmp.m) <- c("Up/On","Down/Off")
colnames(tmp.m) <- c("Up/Down","On/Off")
tmp.df <- as.data.frame(tmp.m)
tmp.df <- cbind(rownames(tmp.df), tmp.df)
tmp.df <- melt(tmp.df)
colnames(tmp.df) <- c("Status", "Group", "Value")
p <- ggbarplot(tmp.df, x = "Group", y = "Value", fill= "Status", cex.main=0.1 , cex.axis=0.1,  label = TRUE,   lab.col = "white",lab.pos = "in", lab.size = 2.5, palette = c( "#0073C2FF", "#EFC000FF"), color="white", title=cancers.v[c], ylab="", xlab="", ggtheme = theme_bw()) + theme(legend.position="left", plot.title = element_text(hjust = 0.5))
pcom.lv[[c]] <- p
}
ggarrange(pcom.lv[[1]], pcom.lv[[2]], pcom.lv[[3]], pcom.lv[[4]], pcom.lv[[5]], pcom.lv[[6]], pcom.lv[[7]], pcom.lv[[8]], pcom.lv[[9]], pcom.lv[[10]], pcom.lv[[11]], pcom.lv[[12]], pcom.lv[[13]], pcom.lv[[14]], pcom.lv[[15]], pcom.lv[[16]], pcom.lv[[17]], pcom.lv[[18]], nrow=2, ncol=9, legend="right", common.legend = TRUE)
ggsave("Lnc_Differential_Exp_Statatus_by_Cancer1.pdf", width=16,height=6);


###########
pcom.lv <- list()
for(c in 1:nrow(DELstat.m)) {
tmp.m <- matrix( as.vector(DELstat.m[c,c(1,2,3,4,5,6)]), nrow=2) 
rownames(tmp.m) <- c("Up/On","Down/Off")
colnames(tmp.m) <- c("Up/Down","On/Off", "UPON/DOWNOFF")

tmp.df <- as.data.frame(tmp.m)
tmp.df <- cbind(rownames(tmp.df), tmp.df)
tmp.df <- melt(tmp.df)
colnames(tmp.df) <- c("Status", "Group", "Value")
p <- ggbarplot(tmp.df, x = "Group", y = "Value", fill= "Status",   label = TRUE, ylim=c(0,1), palette = c( "#0073C2FF", "#EFC000FF"), color="white", title=cancers.v[c], ylab="", xlab="", ggtheme = theme_bw()) + theme(legend.position="none", plot.title = element_text(hjust = 0.5))
pcom.lv[[c]] <- p
}
ggarrange(pcom.lv[[1]], pcom.lv[[2]], pcom.lv[[3]], pcom.lv[[4]], pcom.lv[[5]], pcom.lv[[6]], pcom.lv[[7]], pcom.lv[[8]], pcom.lv[[9]], pcom.lv[[10]], pcom.lv[[11]], pcom.lv[[12]], pcom.lv[[13]], pcom.lv[[14]], pcom.lv[[15]], pcom.lv[[16]], pcom.lv[[17]], pcom.lv[[18]], nrow=3, ncol=6, legend=NULL)
ggsave("Lnc_Differential_Exp_Statatus_by_Cancer1.pdf", width=10,height=6.5);


####################
DEL_logP.m <- -log10(DELPvalue.m)
sigP.df <- melt(DEL_logP.m)
colnames(sigP.df) <- c("Cancer", "Status", "Value")

p <- ggplot(sigP.df, aes(x=Cancer, y=Value, group=Status)) + geom_line(aes(colour=Status)) + geom_point(size=1.5, shape=20) + theme_light() 
p <- p + theme(legend.position=c(0.89,0.85),legend.text=element_text(family="Arial"), axis.text.x = element_text(angle = 45, vjust = 1,  size = 8, hjust = 1)) + guides(fill=guide_legend(title=NULL))
p <- p + scale_color_manual( values = c("Up/On"="#EFC000FF", "Down/Off"="#0073C2FF", "All"="red")) 
p <- p + geom_hline(yintercept= -log10(0.05), linetype="dashed", colour="gray")
p <- p + xlab("")+ ylab("-log10(P_value)")
p 
ggsave("Lnc_Differential_Exp_NetworkEnrichP_scatter.pdf", width=7, height=5)

p <- ggplot(sigP.df, aes(x=Cancer, y=Value, group=Status)) + geom_line(aes(colour=Status)) + geom_point(size=1.5, shape=20) + theme_light() 
p <- p + theme(legend.position=c(0.8,0.25),legend.text=element_text(family="Arial"), axis.text.x = element_text(vjust = 1,  size = 8, hjust = 1)) + guides(fill=guide_legend(title=NULL))
p <- p + scale_color_manual( values = c("Up/On"="#EFC000FF", "Down/Off"="#0073C2FF", "All"="red")) 
p <- p + geom_hline(yintercept= -log10(0.05), linetype="dashed", colour="gray")
p <- p + xlab("")+ ylab("-log10(P_value)")
p <- p + coord_flip()  + scale_x_discrete(limits=rev(cancers.v))  + guides(fill=guide_legend(title=NULL))
p
ggsave("Lnc_Differential_Exp_NetworkEnrichP_scatter1.pdf", width=4, height=5)

#################
#################
options(stringsAsFactors=FALSE)
cancers.v  <-  c("BLCA", "BRCA(ER+)", "BRCA(ER-)", "CESC", "CHOL", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC",
				"LUAD", "LUSC", "PAAD" , "PRAD", "READ", "THCA", "UCEC")

load("../CancerRdata/exp_FPKM/BLCA_LNC_expFPKM.Rdata")

LncDECount.m <- matrix(0, nrow=nrow(LNCexpFPKM.m), ncol=18)
rownames(LncDECount.m) <- rownames(LNCexpFPKM.m)
colnames(LncDECount.m) <- cancers.v

LncSwitchCount.m <- matrix(0, nrow=nrow(LNCexpFPKM.m), ncol=18)
rownames(LncSwitchCount.m) <- rownames(LNCexpFPKM.m)
colnames(LncSwitchCount.m) <- cancers.v

LncDESwitchCount.m <- matrix(0, nrow=nrow(LNCexpFPKM.m), ncol=18)
rownames(LncDESwitchCount.m) <- rownames(LNCexpFPKM.m)
colnames(LncDESwitchCount.m) <- cancers.v

##########---------
for(c in 1:length(cancers.v)) {

	print(cancers.v[c])
	load(paste("./MultiOmicsData/", cancers.v[c], "_MultiOmicsProcess.Rdata", sep="") )
	#load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
	#load(paste("./InteractNetworkStat1/", cancers.v[c], "_InteractNetStat.Rdata", sep="") )
	load(paste("./LncDEstatus/", cancers.v[c], "_LncDEStatus.Rdata", sep="") )

	tmp.idx <- match(rownames(DEL.m), rownames(LncDECount.m))
	LncDECount.m[tmp.idx,c] <- 1

	tmp.idx <-  intersect( sort(c(which(LncDEStatus.m[,9] > 2), which(LncDEStatus.m[,9] < 0.5))), which(LncDEStatus.m[,11] < 0.05) )
	switchLnc.v <- rownames(LncDEStatus.m)[tmp.idx]
	tmp.idx <- match(switchLnc.v, rownames(LncSwitchCount.m))
	LncSwitchCount.m[tmp.idx,c] <- 1

	allLnc.v <- union(rownames(DEL.m), switchLnc.v)
	tmp.idx <- match(allLnc.v, rownames(LncDESwitchCount.m))
	LncDESwitchCount.m[tmp.idx,c] <- 1

}


############+++++++++++
for(c in 1:length(cancers.v)) {

	print(cancers.v[c])
	load(paste("./MultiOmicsData/", cancers.v[c], "_MultiOmicsProcess.Rdata", sep="") )
	load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
	#load(paste("./InteractNetworkStat1/", cancers.v[c], "_InteractNetStat.Rdata", sep="") )
	load(paste("./LncDEstatus/", cancers.v[c], "_LncDEStatus.Rdata", sep="") )

	DELnc.v <- intersect(rownames(DEL.m), rownames(InteractLG.m))
	tmp.idx <- match(DELnc.v, rownames(LncDECount.m))
	LncDECount.m[tmp.idx,c] <- 1

	tmp.idx <-  intersect( sort(c(which(LncDEStatus.m[,9] > 2), which(LncDEStatus.m[,9] < 0.5))), which(LncDEStatus.m[,11] < 0.05) )
	switchLnc.v <- rownames(LncDEStatus.m)[tmp.idx]
	switchLnc.v <- intersect(switchLnc.v, rownames(InteractLG.m))
	tmp.idx <- match(switchLnc.v, rownames(LncSwitchCount.m))
	LncSwitchCount.m[tmp.idx,c] <- 1

	allLnc.v <- union(DELnc.v, switchLnc.v)
	allLnc.v <- intersect(allLnc.v, rownames(InteractLG.m))
	tmp.idx <- match(allLnc.v, rownames(LncDESwitchCount.m))
	LncDESwitchCount.m[tmp.idx,c] <- 1

}


rowSums(LncDECount.m)
hist(rowSums(LncDECount.m))
tmp.idx <- which(rowSums(LncDECount.m) > 0)
LncDECount.m <- LncDECount.m[tmp.idx,]
hist(rowSums(LncDECount.m))

rowSums(LncSwitchCount.m)
hist(rowSums(LncSwitchCount.m))
tmp.idx <- which(rowSums(LncSwitchCount.m) > 0)
LncSwitchCount.m <- LncSwitchCount.m[tmp.idx,]
hist(rowSums(LncSwitchCount.m))

rowSums(LncDESwitchCount.m)
hist(rowSums(LncDESwitchCount.m))
tmp.idx <- which(rowSums(LncDESwitchCount.m) > 0)
LncDESwitchCount.m <- LncDESwitchCount.m[tmp.idx,]
hist(rowSums(LncDESwitchCount.m))

require("RColorBrewer")
display.brewer.all()
brewer.pal(9, "Purples")
#colors.v <- colorRampPalette(brewer.pal(9, "Greens"))(16)
#colors.v <- colorRampPalette(brewer.pal(9, "Blues"))(16)
#colors.v <- colorRampPalette(brewer.pal(9, "BuPu"))(16)
colors.v <- colorRampPalette(brewer.pal(9, "BuGn"))(16)

pdf("Lnc_DExpSwitch_PancerStat.pdf", width=5, height=4)
par(mar=c(4,4,3,1))
barplot(table(rowSums(LncDESwitchCount.m)), col=colors.v, ylab="Number of Genes", xlab="Number of Cancer Types")
dev.off()

###
prop.df <-as.data.frame(table(rowSums(LncDECount.m)), stringsAsFactors=F)
colnames(prop.df) <- c("Cancers", "Freq")
#ggbarplot(prop.df, x="Cancers", y="Freq", fill=colors.v, color="gray", width=0.95, xlab="Number of Cancer Types",ylab="Number of Genes", ggtheme = theme_bw())
p <- ggbarplot(prop.df, x="Cancers", y="Freq", fill="Cancers", palette=colors.v, color="gray",  width=0.95, xlab="Number of Cancer Types", ylab="Number of Genes", label = TRUE,  lab.size = 2, ggtheme = theme_bw()) + theme(legend.position="none")
p
ggsave("Lnc_DExp_PancerStat.pdf", width=4, height=3.5)

###
prop.df <-as.data.frame(table(rowSums(LncSwitchCount.m)), stringsAsFactors=F)
colnames(prop.df) <- c("Cancers", "Freq")
p <- ggbarplot(prop.df, x="Cancers", y="Freq", fill="Cancers", palette=colorRampPalette(brewer.pal(9, "BuGn"))(11), color="gray",  width=0.95, xlab="Number of Cancer Types", ylab="Number of Genes", label = TRUE,  lab.size = 2, ggtheme = theme_bw()) + theme(legend.position="none")
p
ggsave("Lnc_Switch_PancerStat.pdf", width=4, height=3.5)

###
prop.df <-as.data.frame(table(rowSums(LncDESwitchCount.m)), stringsAsFactors=F)
colnames(prop.df) <- c("Cancers", "Freq")
p <- ggbarplot(prop.df, x="Cancers", y="Freq", fill="Cancers", palette=colors.v, color="gray",  width=0.95, xlab="Number of Cancer Types",ylab="Number of Genes", label = TRUE,  lab.size = 2, ggtheme = theme_bw()) + theme(legend.position="none")
p
ggsave("Lnc_DExpSwitch_PancerStat.pdf", width=4, height=3.5)

###########
###########
table(rowSums(LncDESwitchCount.m))

#colors.v <- colorRampPalette(brewer.pal(9, "BuGn"))(7)[1:7]
#colors.v <- brewer.pal(7, "BuGn")

colors.v <- colorRampPalette(brewer.pal(9, "BuGn"))(10)[2:8]

prop.df <- data.frame(category=as.character(c("1","2","3","4","5","6~10",">=11")), Num=NA)
prop.df[1,2] <- length(which(rowSums(LncDECount.m) == 1))
prop.df[2,2] <- length(which(rowSums(LncDECount.m) == 2))
prop.df[3,2] <- length(which(rowSums(LncDECount.m) == 3))
prop.df[4,2] <- length(which(rowSums(LncDECount.m) == 4))
prop.df[5,2] <- length(which(rowSums(LncDECount.m) == 5))
prop.df[6,2] <- length( intersect(which(rowSums(LncDECount.m) > 5), which(rowSums(LncDECount.m)  <= 10) ) )
prop.df[7,2] <- length( intersect(which(rowSums(LncDECount.m) > 10), which(rowSums(LncDECount.m)  <= 20) ) )
#prop.df[8,2] <- length(which(rowSums(LncDECount.m) >= 15))
prop.df <- within(prop.df, category <- factor(category, levels = c("1","2","3","4","5","6~10",">=11")))

prop.v <- as.numeric(prop.df[,2])
names(prop.v) <- prop.df[,1];
lbls <- paste(names(prop.v), "\n(", prop.v, ")", sep= "");

pdf("Lnc_DExp_PancerStat_pie.pdf", width=3, height=3)
par(family="Arial")
par(mar=c(0.28,0.2,0.1,0.52))
pie(prop.v, labels=lbls, col=colors.v, clockwise=TRUE, radius = 0.6, cex= 0.8, font = 2  )
dev.off()

library(ggpubr)

ggpie(prop.df, "Num", label = lbls, fill = "category",  color = "white", palette = colors.v) + theme(legend.position="none")
ggdonutchart(prop.df, "Num", label = rev(lbls), fill = "category",  color = "white", palette = colors.v ) + theme(legend.position="none")
ggsave("Lnc_DExp_PancerStat_donutchart.pdf", width=4, height=3.5)

p1 <- ggdonutchart(prop.df, "Num",  fill = "category",  color = "white", palette = colors.v ) 


##########
prop.df <- data.frame(category=c("1","2","3","4","5","6~10",">=11"), Num=NA)
prop.df[1,2] <- length(which(rowSums(LncSwitchCount.m) == 1))
prop.df[2,2] <- length(which(rowSums(LncSwitchCount.m) == 2))
prop.df[3,2] <- length(which(rowSums(LncSwitchCount.m) == 3))
prop.df[4,2] <- length(which(rowSums(LncSwitchCount.m) == 4))
prop.df[5,2] <- length(which(rowSums(LncSwitchCount.m) == 5))
prop.df[6,2] <- length( intersect(which(rowSums(LncSwitchCount.m) >= 6), which(rowSums(LncSwitchCount.m)  <= 10) ) )
prop.df[7,2] <- length( intersect(which(rowSums(LncSwitchCount.m) > 10), which(rowSums(LncSwitchCount.m)  <= 20) ) )
#prop.df[8,2] <- length(which(rowSums(LncDECount.m) >= 15))
prop.df <- within(prop.df, category <- factor(category, levels = c("1","2","3","4","5","6~10",">=11")))

prop.v <- as.numeric(prop.df[,2])
names(prop.v) <- prop.df[,1];
lbls <- paste(names(prop.v), "\n(", prop.v, ")", sep= "");

pdf("Lnc_Switch_PancerStat_pie.pdf", width=5, height=5)
par(family="Arial")
par(mar=c(0.28,0.2,0.1,0.52))
pie(prop.v, labels=lbls, col=colors.v, clockwise=TRUE, radius = 0.6, cex= 0.8, font = 2  )
dev.off()

#prop.df[7,2] <- 20
ggpie(prop.df, "Num", label = lbls, fill = "category",  color = "white", palette = colors.v) + theme(legend.position="none")
ggdonutchart(prop.df, "Num", label = rev(lbls), fill = "category",  color = "white", palette = colors.v ) + theme(legend.position="none")
ggsave("Lnc_Switch_PancerStat_donutchart.pdf", width=4, height=3.5)

p2 <- ggdonutchart(prop.df, "Num",  fill = "category",  color = "white", palette = colors.v ) 

ggarrange(p1, p2, nrow=2, ncol=1, legend="right", common.legend = TRUE)

ggsave("Lnc_DiffExp_Switch_PancerStat_donutchart.pdf", width=4, height=6);

#############
prop.df <- data.frame(category=c("1","2","3","4","5","6~10",">=11"), Num=NA)
prop.df[1,2] <- length(which(rowSums(LncDESwitchCount.m) == 1))
prop.df[2,2] <- length(which(rowSums(LncDESwitchCount.m) == 2))
prop.df[3,2] <- length(which(rowSums(LncDESwitchCount.m) == 3))
prop.df[4,2] <- length(which(rowSums(LncDESwitchCount.m) == 4))
prop.df[5,2] <- length(which(rowSums(LncDESwitchCount.m) == 5))
prop.df[6,2] <- length( intersect(which(rowSums(LncDESwitchCount.m) > 5), which(rowSums(LncDESwitchCount.m)  <= 10) ) )
prop.df[7,2] <- length( intersect(which(rowSums(LncDESwitchCount.m) > 10), which(rowSums(LncDESwitchCount.m)  <= 20) ) )
#prop.df[8,2] <- length(which(rowSums(LncDECount.m) >= 15))
prop.df <- within(prop.df, category <- factor(category, levels = c("1","2","3","4","5","6~10",">=11")))

prop.v <- as.numeric(prop.df[,2])
names(prop.v) <- prop.df[,1];
lbls <- paste(names(prop.v), "\n(", prop.v, ")", sep= "");

pdf("Lnc_DExpSwitch_PancerStat_pie.pdf", width=3, height=3)
par(family="Arial")
par(mar=c(0.28,0.2,0.1,0.52))
pie(prop.v, labels=lbls, col=colors.v, clockwise=TRUE, radius = 0.6, cex= 0.8, font = 2  )
dev.off()

prop.df[7,2] <- 20
ggpie(prop.df, "Num", label = lbls, fill = "category",  color = "white", palette = colors.v) + theme(legend.position="none")
ggdonutchart(prop.df, "Num", label = rev(lbls), fill = "category",  color = "white", palette = colors.v ) + theme(legend.position="none")
ggsave("Lnc_DExpSwitch_PancerStat_donutchart.pdf", width=4, height=3.5)



################
################
LncDESwitchPanCancerCount.df <- as.data.frame(rowSums(LncDESwitchCount.m))
LncDESwitchPanCancerCount.df$type <- "Moderate"
LncDESwitchPanCancerCount.df$type[which(LncDESwitchPanCancerCount.df[,1] >= 12)] <- "Pan-cancer"
LncDESwitchPanCancerCount.df$type[which(LncDESwitchPanCancerCount.df[,1]  ==  1)] <- "Cancer-specific"
LncDESwitchPanCancerCount.df$class <- "Others"
LncDESwitchPanCancerCount.df$class[which(LncDESwitchPanCancerCount.df[,1] >= 12)] <- "Pan-cancer"
 


tmp.idx <- which(LncDESwitchPanCancerCount.df[,1] >= 12)
RecurrDELnc.v <- rownames(LncDESwitchPanCancerCount.df)[tmp.idx]

load(file="LncPanCancerCount.Rdata")
tmp.idx <- which(LncPanCancerCount.df[,1] >= 15)
PanCancerLnc.v <- rownames(LncPanCancerCount.df)[tmp.idx]

intersect(RecurrDELnc.v , PanCancerLnc.v ) -> PanCancerDELnc.v

data.m <- matrix(0, nrow=nrow(LncPanCancerCount.df),ncol=2)
rownames(data.m) <-  rownames(LncPanCancerCount.df)
colnames(data.m) <- c("RecurrDELnc", "PanCancerLnc")
tmp.idx <- match(RecurrDELnc.v , rownames(data.m))
data.m[tmp.idx,1] <- 1
tmp.idx <- match(PanCancerLnc.v , rownames(data.m))
data.m[tmp.idx,2] <- 1

tmp.m <- table(data.m[,1],data.m[,2])
p.v <- fisher.test(tmp.m)$p.v


library("VennDiagram")
vennDiagram <-venn.diagram(x=list(PanCancer_Lnc_Modulator=PanCancerLnc.v, Recurrent_DE_Lnc=RecurrDELnc.v), filename =NULL, lwd = 3, 
       fill = c("orange", "darkgreen"),  col = "transparent", alpha = 0.4, label.col = "black",  cex = 1.5,
  fontfamily = "Arial",  fontface = "bold",  cat.col = c("black", "black"),  cat.cex = 1,
  cat.fontfamily = "Arial",  cat.fontface = "bold",  margin = 0.05,  cat.dist = c(0.06, 0.07),  cat.pos = c(-10, 15))


pdf("PanCancerLnc_RecurrDELnc.v_vennDiagram12.pdf",width=5,height=3)
plot(1:10, 1:10, col="white",ylab="", xlab="", ann = F, bty = "n", xaxt = "n", yaxt ="n")
grid.draw(vennDiagram)
text(x=5, y=10, labels=paste("Fisher P =", format(p.v, digits=3)))
dev.off()

intersect(PanCancerLnc.v, RecurrDELnc.v) -> PanCancerDELnc.v

c <- 1
load(paste("./LncDEstatus/", cancers.v[c], "_LncDEStatus.Rdata", sep="") )
tmp.idx <- match(PanCancerDELnc.v,  rownames(LncDEStatus.m))
datatmp.m <- LncDEStatus.m[tmp.idx, c(3,5,6)]
colnames(datatmp.m) <- paste0(c("logFC", "P_value", "FDR"), " (", cancers.v[c], ")")

PanCancerDELncStat.df <- as.data.frame(datatmp.m)

for(c in 2:length(cancers.v)) {
	print(cancers.v[c])
	#load(paste("./InteractNetworkStat1/", cancers.v[c], "_InteractNetStat.Rdata", sep="") )
	load(paste("./LncDEstatus/", cancers.v[c], "_LncDEStatus.Rdata", sep="") )
	tmp.idx <- match(PanCancerDELnc.v, rownames(LncDEStatus.m))
	datatmp.m <- LncDEStatus.m[tmp.idx, c(3,5,6)]
	colnames(datatmp.m) <- paste0(c("logFC", "P_value", "FDR"), " (", cancers.v[c], ")")
	PanCancerDELncStat.df <- cbind(PanCancerDELncStat.df, as.data.frame(datatmp.m) )

}

PanCancerDELncStat.df <- format(PanCancerDELncStat.df, digits=3)
tmp.idx <- match(rownames(PanCancerDELncStat.df), LNC.df$gene_id)
anno.df <- LNC.df[tmp.idx, c(11,13)]

PanCancerDELncStat.df <- cbind(anno.df, PanCancerDELncStat.df)
rownames(PanCancerDELncStat.df) <- PanCancerDELnc.v

write.table(PanCancerDELncStat.df, file="PanCancerDELncStat.txt", sep="\t", quote=F, col.names=NA)


####################### valid 85 lncRNAs


library(RColorBrewer)
library(extrafont)
library(ggsci)
library(scales)
mycolors.v <- pal_d3("category20")(18)


########### ENSG00000227036

dataE.df <- data.frame(Exp=0, Type=0, Cancer=0)
dataE.df <- dataE.df[-1,]

for(c in 1:length(cancers.v)) {

print(cancers.v[c])
#load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
load(paste("../CancerRdata/LncValExp0/", cancers.v[c], "_valFPKM.Rdata", sep=""))
row.idx <- which(rownames(LncValFPKM.m) == "ENSG00000227036.5")
PhenoTypesLnc.lv <- PhenoTypes.lv

LncExp.v <- as.vector(LncValFPKM.m[row.idx,])

data.df <- as.data.frame(cbind(LncExp.v, as.vector(PhenoTypesLnc.lv[[1]]) ))
colnames(data.df) <- c("Exp", "Type")
data.df$Exp <- as.numeric(data.df$Exp)
data.df$Type <- ifelse(data.df$Type == "0", "Normal", "Cancer")
data.df$Cancer <- cancers.v[c]

dataE.df <- rbind(dataE.df, data.df)

}


p <- ggboxplot(dataE.df, x="Cancer", y="Exp", fill = "Cancer", color="Type", palette = mycolors.v, add = "none", size=0.2, ylab="Exp", xlab="Cancer Types", ggtheme = theme_bw()) + theme(plot.title = element_text(hjust = 0.5),  axis.text.x=element_blank())
p <- p+ stat_compare_means(aes(group=Type),label = "p.signif")
p <- p+ scale_colour_manual(values=c("darkred", "black"))
#p <- p+ annotate("text", x=3, y=600, label="****  P < 0.0001") 
p
ggsave("LncExp_Compare_by_PanCancer_Distribution_ENSG00000227036.pdf", width=12,height=6)


#dataC.df <- within(dataC.df, Type <- factor(Type, levels = c("Others", "Pan-cancer")))
#dataE.df <- within(dataE.df, Type <- factor(Type, levels = c("Normal", "Cancer")))
p <- ggboxplot(dataE.df, x="Cancer", y="Exp", fill = "Cancer", color="Type", palette = mycolors.v, add = "none", size=0.2, ylab="Exp", xlab="Cancer Types", ggtheme = theme_bw()) + theme(plot.title = element_text(hjust = 0.5),  axis.text.y=element_blank() , legend.position = "left")
p <- p+ stat_compare_means(aes(group=Type),label = "p.signif")
p <- p+ scale_colour_manual(values=c("darkred", "black"))
p <- p+ coord_flip()
#p <- p+ annotate("text", x=15, y=550, label="****  P < 0.0001") 
p <- p+ scale_x_discrete(limits=rev(cancers.v))
p
ggsave("LncExp_Compare_by_PanCancer_Distribution_ENSG00000227036_2.pdf", width=6,height=10)



cox_results.m <- matrix(NA, nrow=nrow(diffgeneFPKM.m), ncol=5)
rownames(cox_results.m) <- rownames(diffgeneFPKM.m)
colnames(cox_results.m) <- c("Zscore",  "HR",  "HR.95L",  "HR.95H",  "P")
######## Cox Regression ##########
library(survival)
library(forestplot)
library(survminer)

ComputeHR <- function(cox.o) {tmp.v <- c("Zscore" = summary(cox.o)$coeff[1,4],
					"HR"=summary(cox.o)$coeff[1,2],
					"lower95"=summary(cox.o)$conf[1,3],
					"upper95"=summary(cox.o)$conf[1,4],
					summary(cox.o)$sctest[3]) }

cox_results.m <- matrix(NA, nrow=18, ncol=5)
rownames(cox_results.m) <- cancers.v
colnames(cox_results.m) <- c("Zscore",  "HR",  "HR.95L",  "HR.95H",  "P")

pcom.lv <- list()
for(c in 1:length(cancers.v)){

  print(cancers.v[c])
  #load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
  load(paste("../CancerRdata/LncValExp0/", cancers.v[c], "_valFPKM.Rdata", sep=""))
  row.idx <- which(rownames(LncValFPKM.m) == "ENSG00000227036.5")
  PhenoTypesLnc.lv <- PhenoTypes.lv

  survival.df <- data.frame("Time" = as.numeric(PhenoTypes.lv$Overall_Survival) / 365.25, "Event" = as.numeric(as.vector(PhenoTypes.lv$Vital_status))  )
  survival.df$Time[which(survival.df$Time > 10)] <-  10
  rownames(survival.df) <- colnames(LncValFPKM.m)

  gene= as.vector(LncValFPKM.m[row.idx,])

  datatmp <- cbind(survival.df[,1:2], gene)
  datatmp <- datatmp[which(PhenoTypesLnc.lv[[1]] == 1), ] ### cancer samples only
  datatmp$group <- ifelse(datatmp$gene > median(datatmp$gene),'high','low') 


  sfit.o <- survfit(Surv( Time, Event ) ~ group, data=datatmp)
  sdiff.o <- survdiff(Surv(Time, Event) ~ group, data = datatmp)
  #cox.o <- coxph(Surv( Time, Event ) ~ group, data=datatmp)
  cox.o <- coxph(Surv( Time, Event ) ~ gene, data=datatmp)

  coeff.v <- ComputeHR(cox.o)
  cox_results.m[c,] <- round(coeff.v, digits=3)
  p.val = 1 - pchisq(sdiff.o$chisq, length(sdiff.o$n) - 1)
  #cox_results.m[i,6] <- p.val

  labels.v <- ifelse(names(sfit.o$strata) == "group=high",  "high",  "low") ### note the order of strata
  labelCol.v <- c("#E7B800", "#2E9FDF")

  p <- ggsurvplot(sfit.o,  data = datatmp,  pval = TRUE,  surv.median.line = "hv", pval.size = 3.5, pval.coord = c(2.1, 0.4), pval.method=TRUE, pval.method.size = 3.5, pval.method.coord = c(0, 0.4),   legend = c(0.20, 0.17),  legend.title = "Group",  legend.labs = c(paste0(labels.v[1]," (n= ", sfit.o$n[1], ")"), paste0(labels.v[2]," (n= ", sfit.o$n[2], ")") ),  conf.int = TRUE,  conf.int.style = "ribbon",  conf.int.alpha= 0.1,  palette =labelCol.v,  risk.table = FALSE,  tables.height = 0.2,  tables.theme = theme_cleantable(),  cumevents=FALSE,  cumcensor=FALSE,  ncensor.plot=FALSE,  ggtheme = theme_bw() )
  
  pcom.lv[[c]] <- p

}

#ggarrange(pcom.lv[[1]], pcom.lv[[2]], pcom.lv[[3]], pcom.lv[[4]], pcom.lv[[5]], pcom.lv[[6]], pcom.lv[[7]], pcom.lv[[8]], pcom.lv[[9]], pcom.lv[[10]], pcom.lv[[11]], pcom.lv[[12]], pcom.lv[[13]], pcom.lv[[14]], pcom.lv[[15]], pcom.lv[[16]], pcom.lv[[17]], pcom.lv[[18]], nrow=5, ncol=4,  common.legend=TRUE)
#p <- pcom.lv[[9]]
#c <- 9
ggsave(paste0("Lnc_ENSG00000227036", "_KMplot_KIRC.pdf"), width=4, height=3.5)


c <- 9  ####### KIRC
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
row.idx <- which(rownames(InteractLG.m) == "ENSG00000227036.5")
col.idx <- which(InteractLG.m[row.idx,] == 1)
cancerTarget.v <- colnames(InteractLG.m)[col.idx]
#[1] "11187"  "11237"  "25984"  "3866"   "5329"   "730755" "84958"

##
target.v <- c()
for(c in 1:length(cancers.v)){

  print(cancers.v[c])
  load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
  row.idx <- which(rownames(InteractLG.m) == "ENSG00000227036.5")
  col.idx <- which(InteractLG.m[row.idx,] == 1)
  target.v <- c(target.v , colnames(InteractLG.m)[col.idx])

}


c <- 9 
####### gene level correlation #######
load(paste("../CancerRdata/LncValExp0/", cancers.v[c], "_valFPKM.Rdata", sep=""))
row.idx <- which(rownames(LncValFPKM.m) == "ENSG00000227036.5")
LncExp.v <- as.vector(LncValFPKM.m[row.idx,])
PhenoTypesLnc.lv <- PhenoTypes.lv

load(paste("../CancerRdata/DNAm_gene/", cancers.v[c], "_DNAmGenelevel.Rdata", sep=""))
row.idx <- which(rownames(avbeta.m) == "3866")
targetDNAm.v <- as.vector(avbeta.m[row.idx,])
PhenoTypesDNAm.lv <- PhenoTypes.lv
pheno.v <- as.numeric(as.vector(PhenoTypesDNAm.lv[[1]]))
pheno.v <- ifelse(pheno.v  == "0", "Normal", "Cancer")

interSample.v <- intersect(colnames(LncValFPKM.m), colnames(avbeta.m))
tmp1.idx <- match(interSample.v, colnames(LncValFPKM.m))
tmp2.idx <- match(interSample.v, colnames(avbeta.m))

LncExp.v <- LncExp.v[tmp1.idx]
targetDNAm.v <- targetDNAm.v[tmp2.idx]
pheno.v <- pheno.v[tmp2.idx]


mydata <- data.frame("DNAm" = as.vector(targetDNAm.v), "Exp" = as.vector(LncExp.v), "Sample" = pheno.v)
ggscatter(data=mydata, x = "Exp", y = "DNAm", title = "", color = "#426671", size =1.2,  add = "reg.line", add.params = list(color = "#764C29", fill = "#E7E1D7"), conf.int = TRUE,   cor.coef = TRUE,  cor.coeff.args = list(method = "pearson", label.x = 0.15,  label.y = 0.25, label.sep = "\n")) +  theme_bw() +  theme( plot.title = element_text(hjust = 0.5))
ggsave("Lnc_ENSG00000227036_exp_Gene3866_DNAm_Cor_KIRC.pdf", width=5, height=5)

ggscatterhist(data=mydata, x = "Exp", y = "DNAm", color = "Sample" , margin.plot = "boxplot", conf.int = TRUE,   cor.coef = TRUE,  cor.coeff.args = list(method = "pearson", label.x = 0.15,  label.y = 0.25, label.sep = "\n"), margin.params = list(fill = "Sample", color = "black", size = 0.2), margin.ggtheme = theme_void(), ggtheme = theme_bw())

library("ggstatsplot")
ggscatterstats(data=mydata, x = "Exp", y = "DNAm", marginal.type="boxplot")

####### probe level correlation #######
load("../CancerRdata/DNAm_gene/probe450kfemanno.rda")
cpgs.idx <- which(probe450kfemanno$eid == "3866")
probeID.v <- probe450kfemanno$probeID[cpgs.idx]
probeID.v <- probeID.v[-5]

load(paste("../CancerRdata/LncValExp0/", cancers.v[c], "_valFPKM.Rdata", sep=""))
row.idx <- which(rownames(LncValFPKM.m) == "ENSG00000227036.5")
LncExp.v <- as.vector(LncValFPKM.m[row.idx,])
PhenoTypesLnc.lv <- PhenoTypes.lv

load(paste("../CancerRdata/DNAm_450k/", cancers.v[c], "_bmiqbeta.Rdata", sep=""))
row.idx <- match(probeID.v, rownames(bmiqbeta.m))
row.idx <- row.idx[-2]  ## missing value
probeID.v <- probeID.v[-2]  ## missing value

genebeta.m <- bmiqbeta.m[row.idx,]

PhenoTypesDNAm.lv <- PhenoTypes.lv
pheno.v <- as.numeric(as.vector(PhenoTypesDNAm.lv[[1]]))
pheno.v <- ifelse(pheno.v  == "0", "Normal", "Cancer")

interSample.v <- intersect(colnames(LncValFPKM.m), colnames(genebeta.m))
tmp1.idx <- match(interSample.v, colnames(LncValFPKM.m))
tmp2.idx <- match(interSample.v, colnames(genebeta.m))

pheno.v <- pheno.v[tmp2.idx]
LncExp.v <- LncExp.v[tmp1.idx]

targetDNAm.v <- genebeta.m[1,tmp2.idx]
mydata <- data.frame("DNAm" = as.vector(targetDNAm.v), "Exp" = as.vector(LncExp.v), "Sample" = pheno.v)


pcom.lv <- list()
for(i in 1:nrow(genebeta.m)){

targetDNAm.v <- genebeta.m[i,tmp2.idx]
mydata <- data.frame("DNAm" = as.vector(targetDNAm.v), "Exp" = as.vector(LncExp.v), "Sample" = pheno.v)

#p <- ggscatterhist(data=mydata, x = "Exp", y = "DNAm", color = "Sample", title=rownames(genebeta.m)[i], margin.plot = "boxplot", conf.int = TRUE,   cor.coef = TRUE,  cor.coeff.args = list(method = "pearson", label.x = 0.15,  label.y = 0.25, label.sep = "\n"), margin.params = list(fill = "Sample", color = "black", size = 0.2), margin.ggtheme = theme_void(), ggtheme = theme_bw())

p <- ggscatter(data=mydata, x = "Exp", y = "DNAm", title = rownames(genebeta.m)[i],  color = "Sample", size =1.2, ylim=c(0,1) , add = "reg.line", add.params = list(color = "#764C29", fill = "#E7E1D7"),  conf.int = TRUE,   cor.coef = TRUE,  cor.coeff.args = list(method = "pearson", label.x = 0.15,  label.y = 0.15, label.sep = "\n")) +  theme_bw() +  theme( plot.title = element_text(hjust = 0.5))
#p <- ggscatter(data=mydata, x = "Exp", y = "DNAm")
pcom.lv[[i]] <- p

}

ggarrange(pcom.lv[[1]], pcom.lv[[2]], pcom.lv[[3]], pcom.lv[[4]], pcom.lv[[5]], pcom.lv[[6]], nrow=2, ncol=3,  legend="bottom", common.legend=TRUE)
ggsave("Lnc_ENSG00000227036_exp_Gene3866_cpgs_DNAm_Cor_KIRC.pdf", width=12,height=8);



####### gene level correlation  in pan-cancer#######
pcom.lv <- list()
for(c in  1:length(cancers.v)){
#if (c != 14)
print(cancers.v[c])
load(paste("../CancerRdata/LncValExp0/", cancers.v[c], "_valFPKM.Rdata", sep=""))
row.idx <- which(rownames(LncValFPKM.m) == "ENSG00000227036.5")
LncExp.v <- as.vector(LncValFPKM.m[row.idx,])
PhenoTypesLnc.lv <- PhenoTypes.lv

load(paste("../CancerRdata/DNAm_gene/", cancers.v[c], "_DNAmGenelevel.Rdata", sep=""))
row.idx <- which(rownames(avbeta.m) == "3866")
targetDNAm.v <- as.vector(avbeta.m[row.idx,])
PhenoTypesDNAm.lv <- PhenoTypes.lv
pheno.v <- as.numeric(as.vector(PhenoTypesDNAm.lv[[1]]))
pheno.v <- ifelse(pheno.v  == "0", "Normal", "Cancer")

interSample.v <- intersect(colnames(LncValFPKM.m), colnames(avbeta.m))
tmp1.idx <- match(interSample.v, colnames(LncValFPKM.m))
tmp2.idx <- match(interSample.v, colnames(avbeta.m))

LncExp.v <- LncExp.v[tmp1.idx]
targetDNAm.v <- targetDNAm.v[tmp2.idx]
pheno.v <- pheno.v[tmp2.idx]


mydata <- data.frame("DNAm" = as.vector(targetDNAm.v), "Exp" = as.vector(LncExp.v), "Sample" = pheno.v)
p <- ggscatter(data=mydata, x = "Exp", y = "DNAm", title = cancers.v[c], color = "#426671", size =1.2, ylim=c(0,1), add = "reg.line", add.params = list(color = "#764C29", fill = "#E7E1D7"), conf.int = TRUE,   cor.coef = TRUE,  cor.coeff.args = list(method = "pearson", label.x = 0.15,  label.y = 0.15, label.sep = "\n")) +  theme_bw() +  theme( plot.title = element_text(hjust = 0.5))

pcom.lv[[c]] <- p
#}
}

#ggarrange(pcom.lv[[1]], pcom.lv[[2]], pcom.lv[[3]], pcom.lv[[4]], pcom.lv[[5]], pcom.lv[[6]], pcom.lv[[7]], pcom.lv[[8]], pcom.lv[[9]], pcom.lv[[10]], pcom.lv[[11]], pcom.lv[[12]], pcom.lv[[13]], pcom.lv[[14]], pcom.lv[[15]], pcom.lv[[16]], pcom.lv[[17]], pcom.lv[[18]], nrow=3, ncol=6, legend=NULL)

ggarrange(pcom.lv[[1]], pcom.lv[[2]], pcom.lv[[3]], pcom.lv[[4]], pcom.lv[[5]], pcom.lv[[6]], pcom.lv[[7]], pcom.lv[[8]], pcom.lv[[10]],    pcom.lv[[12]], pcom.lv[[13]], pcom.lv[[14]], pcom.lv[[15]], pcom.lv[[17]], pcom.lv[[18]], nrow=5, ncol=3, legend=NULL)

ggsave("Lnc_ENSG00000227036_exp_Gene3866_DNAm_Cor_pan_cancer.pdf",  width=10,height=13)


###################### 
########### ENSG00000203499.9

dataE.df <- data.frame(Exp=0, Type=0, Cancer=0)
dataE.df <- dataE.df[-1,]

for(c in 1:length(cancers.v)) {

print(cancers.v[c])
#load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
load(paste("../CancerRdata/LncValExp0/", cancers.v[c], "_valFPKM.Rdata", sep=""))
row.idx <- which(rownames(LncValFPKM.m) == "ENSG00000203499.9")
PhenoTypesLnc.lv <- PhenoTypes.lv

LncExp.v <- as.vector(LncValFPKM.m[row.idx,])

data.df <- as.data.frame(cbind(LncExp.v, as.vector(PhenoTypesLnc.lv[[1]]) ))
colnames(data.df) <- c("Exp", "Type")
data.df$Exp <- as.numeric(data.df$Exp)
data.df$Type <- ifelse(data.df$Type == "0", "Normal", "Cancer")
data.df$Cancer <- cancers.v[c]

dataE.df <- rbind(dataE.df, data.df)

}


p <- ggboxplot(dataE.df, x="Cancer", y="Exp", fill = "Cancer", color="Type", palette = mycolors.v, add = "none", size=0.2, ylab="Exp", xlab="Cancer Types", ggtheme = theme_bw()) + theme(plot.title = element_text(hjust = 0.5),  axis.text.x=element_blank())
p <- p+ stat_compare_means(aes(group=Type),label = "p.signif")
p <- p+ scale_colour_manual(values=c("darkred", "black"))
#p <- p+ annotate("text", x=3, y=600, label="****  P < 0.0001") 
p
ggsave("LncExp_Compare_by_PanCancer_Distribution_ENSG00000203499.pdf", width=12,height=6)




#dataC.df <- within(dataC.df, Type <- factor(Type, levels = c("Others", "Pan-cancer")))
#dataE.df <- within(dataE.df, Type <- factor(Type, levels = c("Normal", "Cancer")))
p <- ggboxplot(dataE.df, x="Cancer", y="Exp", fill = "Cancer", color="Type", palette = mycolors.v, add = "none", size=0.2, ylab="Exp", xlab="Cancer Types", ggtheme = theme_bw()) + theme(plot.title = element_text(hjust = 0.5),  axis.text.y=element_blank() , legend.position = "left")
p <- p+ stat_compare_means(aes(group=Type),label = "p.signif")
p <- p+ scale_colour_manual(values=c("darkred", "black"))
p <- p+ coord_flip()
#p <- p+ annotate("text", x=15, y=550, label="****  P < 0.0001") 
p <- p+ scale_x_discrete(limits=rev(cancers.v))
p
ggsave("LncExp_Compare_by_PanCancer_Distribution_ENSG00000203499_2.pdf", width=6,height=10)



cox_results.m <- matrix(NA, nrow=nrow(diffgeneFPKM.m), ncol=5)
rownames(cox_results.m) <- rownames(diffgeneFPKM.m)
colnames(cox_results.m) <- c("Zscore",  "HR",  "HR.95L",  "HR.95H",  "P")
######## Cox Regression ##########
library(survival)
library(forestplot)
library(survminer)

ComputeHR <- function(cox.o) {tmp.v <- c("Zscore" = summary(cox.o)$coeff[1,4],
					"HR"=summary(cox.o)$coeff[1,2],
					"lower95"=summary(cox.o)$conf[1,3],
					"upper95"=summary(cox.o)$conf[1,4],
					summary(cox.o)$sctest[3]) }

cox_results.m <- matrix(NA, nrow=18, ncol=5)
rownames(cox_results.m) <- cancers.v
colnames(cox_results.m) <- c("Zscore",  "HR",  "HR.95L",  "HR.95H",  "P")

pcom.lv <- list()
for(c in 1:length(cancers.v)){

  print(cancers.v[c])
  #load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
  load(paste("../CancerRdata/LncValExp0/", cancers.v[c], "_valFPKM.Rdata", sep=""))
  row.idx <- which(rownames(LncValFPKM.m) == "ENSG00000203499.9")
  PhenoTypesLnc.lv <- PhenoTypes.lv

  survival.df <- data.frame("Time" = as.numeric(PhenoTypes.lv$Overall_Survival) / 365.25, "Event" = as.numeric(as.vector(PhenoTypes.lv$Vital_status))  )
  survival.df$Time[which(survival.df$Time > 10)] <-  10
  rownames(survival.df) <- colnames(LncValFPKM.m)

  gene= as.vector(LncValFPKM.m[row.idx,])

  datatmp <- cbind(survival.df[,1:2], gene)
  datatmp <- datatmp[which(PhenoTypesLnc.lv[[1]] == 1), ] ### cancer samples only
  datatmp$group <- ifelse(datatmp$gene > median(datatmp$gene),'high','low') 


  sfit.o <- survfit(Surv( Time, Event ) ~ group, data=datatmp)
  sdiff.o <- survdiff(Surv(Time, Event) ~ group, data = datatmp)
  #cox.o <- coxph(Surv( Time, Event ) ~ group, data=datatmp)
  cox.o <- coxph(Surv( Time, Event ) ~ gene, data=datatmp)

  coeff.v <- ComputeHR(cox.o)
  cox_results.m[c,] <- round(coeff.v, digits=3)
  p.val = 1 - pchisq(sdiff.o$chisq, length(sdiff.o$n) - 1)
  #cox_results.m[i,6] <- p.val

  labels.v <- ifelse(names(sfit.o$strata) == "group=high",  "high",  "low") ### note the order of strata
  labelCol.v <- c("#E7B800", "#2E9FDF")

  p <- ggsurvplot(sfit.o,  data = datatmp,  pval = TRUE,  surv.median.line = "hv", pval.size = 3.5, pval.coord = c(6.1, 0.75), pval.method=TRUE, pval.method.size = 3.5, pval.method.coord = c(4.4, 0.75),   legend = c(0.20, 0.17),  legend.title = "Group",  legend.labs = c(paste0(labels.v[1]," (n= ", sfit.o$n[1], ")"), paste0(labels.v[2]," (n= ", sfit.o$n[2], ")") ),  conf.int = TRUE,  conf.int.style = "ribbon",  conf.int.alpha= 0.1,  palette =labelCol.v,  risk.table = FALSE,  tables.height = 0.2,  tables.theme = theme_cleantable(),  cumevents=FALSE,  cumcensor=FALSE,  ncensor.plot=FALSE,  ggtheme = theme_bw() )
  
  pcom.lv[[c]] <- p

}

#ggarrange(pcom.lv[[1]], pcom.lv[[2]], pcom.lv[[3]], pcom.lv[[4]], pcom.lv[[5]], pcom.lv[[6]], pcom.lv[[7]], pcom.lv[[8]], pcom.lv[[9]], pcom.lv[[10]], pcom.lv[[11]], pcom.lv[[12]], pcom.lv[[13]], pcom.lv[[14]], pcom.lv[[15]], pcom.lv[[16]], pcom.lv[[17]], pcom.lv[[18]], nrow=5, ncol=4,  common.legend=TRUE)
#p <- pcom.lv[[14]]
#c <- 14
ggsave(paste0("Lnc_ENSG00000203499", "_KMplot_PAAD.pdf"), width=4, height=3.5)
#c <- 2
ggsave(paste0("Lnc_ENSG00000203499", "_KMplot_BRCA_ERpos.pdf"), width=4, height=3.5)


c <- 14 #### PAAD

load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
row.idx <- which(rownames(InteractLG.m) == "ENSG00000203499.9")
col.idx <- which(InteractLG.m[row.idx,] == 1)
cancerTarget.v <- colnames(InteractLG.m)[col.idx]



##
target.v <- c()
for(c in 1:length(cancers.v)){

  print(cancers.v[c])
  load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
  row.idx <- which(rownames(InteractLG.m) == "ENSG00000203499.9")
  col.idx <- which(InteractLG.m[row.idx,] == 1)
  target.v <- c(target.v , colnames(InteractLG.m)[col.idx])

}



c <- 14 #### PAAD
####### gene level correlation #######
load(paste("../CancerRdata/LncValExp0/", cancers.v[c], "_valFPKM.Rdata", sep=""))
row.idx <- which(rownames(LncValFPKM.m) == "ENSG00000203499.9")
LncExp.v <- as.vector(LncValFPKM.m[row.idx,])
PhenoTypesLnc.lv <- PhenoTypes.lv

load(paste("../CancerRdata/DNAm_gene/", cancers.v[c], "_DNAmGenelevel.Rdata", sep=""))
row.idx <- which(rownames(avbeta.m) == "25900")
targetDNAm.v <- as.vector(avbeta.m[row.idx,])
PhenoTypesDNAm.lv <- PhenoTypes.lv
pheno.v <- as.numeric(as.vector(PhenoTypesDNAm.lv[[1]]))
pheno.v <- ifelse(pheno.v  == "0", "Normal", "Cancer")

interSample.v <- intersect(colnames(LncValFPKM.m), colnames(avbeta.m))
tmp1.idx <- match(interSample.v, colnames(LncValFPKM.m))
tmp2.idx <- match(interSample.v, colnames(avbeta.m))

LncExp.v <- LncExp.v[tmp1.idx]
targetDNAm.v <- targetDNAm.v[tmp2.idx]
pheno.v <- pheno.v[tmp2.idx]


mydata <- data.frame("DNAm" = as.vector(targetDNAm.v), "Exp" = as.vector(LncExp.v), "Sample" = pheno.v)
ggscatter(data=mydata, x = "Exp", y = "DNAm", title = "", color = "#426671", size =1.2,  add = "reg.line", add.params = list(color = "#764C29", fill = "#E7E1D7"), conf.int = TRUE,   cor.coef = TRUE,  cor.coeff.args = list(method = "pearson", label.x = 4,  label.y = 0.25, label.sep = "\n")) +  theme_bw() +  theme( plot.title = element_text(hjust = 0.5))
ggsave("Lnc_ENSG00000203499_exp_Gene25900_DNAm_Cor_PAAD.pdf", width=5, height=5)

ggscatterhist(data=mydata, x = "Exp", y = "DNAm", color = "Sample" , add = "reg.line", margin.plot = "boxplot", conf.int = TRUE,   cor.coef = TRUE,  cor.coeff.args = list(method = "pearson", label.x = 0.15,  label.y = 0.25, label.sep = "\n"), margin.params = list(fill = "Sample", color = "black", size = 0.2), margin.ggtheme = theme_void(), ggtheme = theme_bw())

library("ggstatsplot")
ggscatterstats(data=mydata, x = "Exp", y = "DNAm", marginal.type="boxplot")



####### probe level correlation #######
load("../CancerRdata/DNAm_gene/probe450kfemanno.rda")
cpgs.idx <- which(probe450kfemanno$eid == "25900")
probeID.v <- probe450kfemanno$probeID[cpgs.idx]
probe450kfemanno$GeneGroup[cpgs.idx]
valid.v <- c(1,2,3,4,7, 14, 23, 30, 33)
probeID.v <- probeID.v[valid.v]

load(paste("../CancerRdata/LncValExp0/", cancers.v[c], "_valFPKM.Rdata", sep=""))
row.idx <- which(rownames(LncValFPKM.m) == "ENSG00000203499.9")
LncExp.v <- as.vector(LncValFPKM.m[row.idx,])
PhenoTypesLnc.lv <- PhenoTypes.lv

load(paste("../CancerRdata/DNAm_450k/", cancers.v[c], "_bmiqbeta.Rdata", sep=""))
row.idx <- match(probeID.v, rownames(bmiqbeta.m))
#row.idx <- row.idx[-2]
#probeID.v <- probeID.v[-2]

genebeta.m <- bmiqbeta.m[row.idx,]

PhenoTypesDNAm.lv <- PhenoTypes.lv
pheno.v <- as.numeric(as.vector(PhenoTypesDNAm.lv[[1]]))
pheno.v <- ifelse(pheno.v  == "0", "Normal", "Cancer")

interSample.v <- intersect(colnames(LncValFPKM.m), colnames(genebeta.m))
tmp1.idx <- match(interSample.v, colnames(LncValFPKM.m))
tmp2.idx <- match(interSample.v, colnames(genebeta.m))

pheno.v <- pheno.v[tmp2.idx]
LncExp.v <- LncExp.v[tmp1.idx]

targetDNAm.v <- genebeta.m[1,tmp2.idx]
mydata <- data.frame("DNAm" = as.vector(targetDNAm.v), "Exp" = as.vector(LncExp.v), "Sample" = pheno.v)


pcom.lv <- list()
for(i in 1:nrow(genebeta.m)){

targetDNAm.v <- genebeta.m[i,tmp2.idx]
mydata <- data.frame("DNAm" = as.vector(targetDNAm.v), "Exp" = as.vector(LncExp.v), "Sample" = pheno.v)

#p <- ggscatterhist(data=mydata, x = "Exp", y = "DNAm", color = "Sample", title=rownames(genebeta.m)[i], margin.plot = "boxplot", conf.int = TRUE,   cor.coef = TRUE,  cor.coeff.args = list(method = "pearson", label.x = 0.15,  label.y = 0.25, label.sep = "\n"), margin.params = list(fill = "Sample", color = "black", size = 0.2), margin.ggtheme = theme_void(), ggtheme = theme_bw())

p <- ggscatter(data=mydata, x = "Exp", y = "DNAm", title = rownames(genebeta.m)[i],  color = "Sample", size =1.2, ylim=c(0,1) , add = "reg.line", add.params = list(color = "#764C29", fill = "#E7E1D7"),  conf.int = TRUE,   cor.coef = TRUE,  cor.coeff.args = list(method = "pearson", label.x = 4,  label.y = 0.15, label.sep = "\n")) +  theme_bw() +  theme( plot.title = element_text(hjust = 0.5))
#p <- ggscatter(data=mydata, x = "Exp", y = "DNAm")
pcom.lv[[i]] <- p

}

ggarrange(pcom.lv[[1]], pcom.lv[[2]], pcom.lv[[3]], pcom.lv[[4]], pcom.lv[[5]], pcom.lv[[6]], pcom.lv[[7]], pcom.lv[[8]], pcom.lv[[9]], nrow=3, ncol=3,  legend="bottom", common.legend=TRUE)
ggsave("Lnc_ENSG00000203499_exp_Gene25900_cpgs_DNAm_Cor_PAAD.pdf", width=12,height=12);




####### gene level correlation  in pan-cancer#######
pcom.lv <- list()
for(c in  1:length(cancers.v)){
#if (c != 14)
print(cancers.v[c])
load(paste("../CancerRdata/LncValExp0/", cancers.v[c], "_valFPKM.Rdata", sep=""))
row.idx <- which(rownames(LncValFPKM.m) == "ENSG00000203499.9")
LncExp.v <- as.vector(LncValFPKM.m[row.idx,])
PhenoTypesLnc.lv <- PhenoTypes.lv

load(paste("../CancerRdata/DNAm_gene/", cancers.v[c], "_DNAmGenelevel.Rdata", sep=""))
row.idx <- which(rownames(avbeta.m) == "25900")
targetDNAm.v <- as.vector(avbeta.m[row.idx,])
PhenoTypesDNAm.lv <- PhenoTypes.lv
pheno.v <- as.numeric(as.vector(PhenoTypesDNAm.lv[[1]]))
pheno.v <- ifelse(pheno.v  == "0", "Normal", "Cancer")

interSample.v <- intersect(colnames(LncValFPKM.m), colnames(avbeta.m))
tmp1.idx <- match(interSample.v, colnames(LncValFPKM.m))
tmp2.idx <- match(interSample.v, colnames(avbeta.m))

LncExp.v <- LncExp.v[tmp1.idx]
targetDNAm.v <- targetDNAm.v[tmp2.idx]
pheno.v <- pheno.v[tmp2.idx]


mydata <- data.frame("DNAm" = as.vector(targetDNAm.v), "Exp" = as.vector(LncExp.v), "Sample" = pheno.v)
p <- ggscatter(data=mydata, x = "Exp", y = "DNAm", title = cancers.v[c], color = "#426671", size =1.2,  ylim=c(0,1), add = "reg.line", add.params = list(color = "#764C29", fill = "#E7E1D7"), conf.int = TRUE,   cor.coef = TRUE,  cor.coeff.args = list(method = "pearson", label.x = 3.5,  label.y = 0.2, label.sep = "\n")) +  theme_bw() +  theme( plot.title = element_text(hjust = 0.5))

pcom.lv[[c]] <- p
#}
}

#ggarrange(pcom.lv[[1]], pcom.lv[[2]], pcom.lv[[3]], pcom.lv[[4]], pcom.lv[[5]], pcom.lv[[6]], pcom.lv[[7]], pcom.lv[[8]], pcom.lv[[9]], pcom.lv[[10]], pcom.lv[[11]], pcom.lv[[12]], pcom.lv[[13]], pcom.lv[[14]], pcom.lv[[15]], pcom.lv[[16]], pcom.lv[[17]], pcom.lv[[18]], nrow=3, ncol=6, legend=NULL)

ggarrange(pcom.lv[[1]], pcom.lv[[2]], pcom.lv[[3]], pcom.lv[[4]], pcom.lv[[5]], pcom.lv[[6]], pcom.lv[[7]], pcom.lv[[8]],  pcom.lv[[11]], pcom.lv[[12]], pcom.lv[[13]], pcom.lv[[15]], pcom.lv[[17]], pcom.lv[[18]], nrow=5, ncol=3, legend=NULL)

ggsave("Lnc_ENSG00000203499_exp_Gene25900_DNAm_Cor_pan_cancer.pdf",  width=10,height=13)



##################
##################
##################
##################
intersect(PanCancerLnc.v, RecurrDELnc.v) -> validID.v
load("~/backup/TCGA-XenaBrowser/gencode/gencode.Rdata")
tmp.idx <- match(validID.v, LNC.df$gene_id)
LNC.df[tmp.idx,]


load("../CancerRdata/Exp_gene/BLCA_expGenelevel.Rdata")
backgroundGenes.v  <- rownames(PCGexpFPKM.m)

Cosmic.df <- read.csv("./CancerHallMark/cancer_gene_census.csv", header=T,check.names=F)
backgroundInter.v <- intersect(Cosmic.df[,3], backgroundGenes.v)


FisherP.m <- matrix(NA, nrow=60, ncol=18)
#rownames(FisherP.m) <- cancers.v
#colnames(FisherP.m) <- "COSMIC"


for(c in 1:length(cancers.v)) {
#rm(list=ls())
print(cancers.v[c])
#load(paste("./MultiOmicsData/", cancers.v[c], "_MultiOmicsProcess.Rdata", sep="") )
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
#load(paste("./LncDEstatus/", cancers.v[c], "_LncDEStatus.Rdata", sep="") )
    for(g in 1:length(validID.v)){
    	tmp.idx <- which( rownames(InteractLG.m) == validID.v[g] )
   	CancerGenes.idx <- which(InteractLG.m[tmp.idx,] == 1)
   	CancerGenes.v <- colnames(InteractLG.m)[CancerGenes.idx]
    	CancerGenesInter.v <- intersect( CancerGenes.v ,Cosmic.df[,3])
	
	if(length(CancerGenesInter.v) >= 1){
		data.m <- matrix(NA, nrow=2, ncol=2)
		data.m[1,1] <- length(backgroundGenes.v)
		data.m[2,1] <- length(backgroundInter.v)
		data.m[1,2] <- length(CancerGenes.v )
		data.m[2,2] <- length(CancerGenesInter.v)
		print(fisher.test(data.m)$p.value)
        	FisherP.m[g, c] <- fisher.test(data.m)$p.value
	}
    }
}


library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(topGO)
library(VennDiagram)

intersect(PanCancerLnc.v, RecurrDELnc.v) -> validLnc.v
#####
enrichedTerm.v <-  c()
enrichedTrem.lv <- list()
for(c in 1:length(cancers.v)) {
#rm(list=ls())
print(cancers.v[c])

load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
tmp.idx <- match(validLnc.v, rownames(InteractLG.m))
validInterAct.m <- InteractLG.m[setdiff(tmp.idx, NA), ]
validInterAct.v <- colSums(validInterAct.m)
validTarget.v <- names(validInterAct.v)[which(validInterAct.v > 0)]
#######GO 
ego <- enrichGO(validTarget.v, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.01, 
                 qvalueCutoff = 0.05, keyType = 'ENTREZID')
ego2 <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
#dim(ego2)
#ego2 <- ego
ego3 <- setReadable(ego2, OrgDb = org.Hs.eg.db)
dim(ego3)
#head(summary(ego3))
enrichedResult.df <- as.data.frame(ego3)
enrichedResult.df$FDR <- -log10(enrichedResult.df$p.adjust)

enrichedTerm.v <- c(enrichedTerm.v , enrichedResult.df$Description)
enrichedTrem.lv[[c]] <-  enrichedResult.df
}

rev(sort(table(enrichedTerm.v)))


load("../CancerRdata/Exp_gene/BLCA_expGenelevel.Rdata")
allGenes.v <- rownames(PCGexpFPKM.m)

enrichedTerm.v <-  c()
enrichedTrem.lv <- list()

for(c in 1:length(cancers.v)) {
#rm(list=ls())
print(cancers.v[c])

load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
tmp.idx <- match(validLnc.v, rownames(InteractLG.m))
validInterAct.m <- InteractLG.m[setdiff(tmp.idx, NA), ]
validInterAct.v <- colSums(validInterAct.m)
validTarget.v <- names(validInterAct.v)[which(validInterAct.v > 0)]
#######GO 
geneList <- factor(as.integer(allGenes.v %in% validTarget.v))
names(geneList) <- allGenes.v 
str(geneList)

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, nodeSize = 10, annot = annFUN.org, mapping = "org.Hs.eg.db", ID="entrez")
result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
sig.tab <- GenTable(GOdata, Fis = result,  topNodes = 100)


sig.tab$p.adjust  <- p.adjust(sig.tab$Fis, method = "BH")
sig.tab$FDR <-  -log10(as.numeric(sig.tab$Fis))


enrichedTerm.v <- c(enrichedTerm.v , sig.tab$Term)
enrichedTrem.lv[[c]] <-  sig.tab

}

rev(sort(table(enrichedTerm.v)))


colors.v <- colorRampPalette(brewer.pal(9, "BuGn"))(10)[2:8]

goterms.v <- c("cell adhesion", "inflammatory response", "cell proliferation", "cell migration", "cell motility", "cell communication", "signal transduction", "epithelium development", "cell differentiation", "lymphocyte activation")
enrichP.m <- matrix(0, nrow=18, ncol=10)
rownames(enrichP.m)  <- cancers.v
colnames(enrichP.m)  <- goterms.v


for(c in 1:length(cancers.v)) {
#rm(list=ls())
print(cancers.v[c])
tmp.idx <- match(goterms.v , enrichedTrem.lv[[c]]$Term)
col.idx <- which(!is.na(tmp.idx))
fdr.v <- as.numeric(enrichedTrem.lv[[c]][setdiff(tmp.idx, NA),8])
enrichP.m[c, col.idx] <- fdr.v

}

data.m <- melt(enrichP.m)
colnames(data.m)  <- c("Cancer", "GOterm", "FDR")

p <- ggplot(data.m, aes(x=GOterm,y=Cancer)) + geom_point(aes( size = FDR), shape=21, colour="white", fill=brewer.pal(9, "BuGn")[6])
p <- p + theme_bw() + theme(axis.text.x = element_text(angle=30, hjust=1,vjust=1), text=element_text(family="Arial"), panel.grid.major=element_blank())
p <- p + xlab("") + ylab("")
p <- p + labs(size="-log10(FDR)")
p <- p + scale_y_discrete(limits=rev(cancers.v))
p
ggsave("Lnc_85_enriched_GO_term_bubble.pdf", width=5, height=6)



######## kegg
ego <- enrichKEGG(gene = validTarget.v, keyType = "kegg", organism = 'hsa', pvalueCutoff = 0.01,  pAdjustMethod= "BH", qvalueCutoff  = 0.05)
enrichedResult.df <- as.data.frame(ego)
enrichedResult.df$FDR <- -log10(enrichedResult.df$p.adjust)

######## topGO
load("../CancerRdata/Exp_gene/BLCA_expGenelevel.Rdata")
allGenes.v <- rownames(PCGexpFPKM.m)

geneList <- factor(as.integer(allGenes.v %in% validTarget.v))
names(geneList) <- allGenes.v 
str(geneList)

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, nodeSize = 10, annot = annFUN.org, mapping = "org.Hs.eg.db", ID="entrez")
result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
sig.tab <- GenTable(GOdata, Fis = result,  topNodes = 100)


sig.tab$p.adjust  <- p.adjust(sig.tab$Fis, method = "BH")*10^2
sig.tab$FDR <-  -log10(as.numeric(sig.tab$Fis))
row.idx <- c(6,11,21,35,38,47,48)
sig.tab  <- sig.tab[row.idx,]




###############
tmp.idx <- which(LncDESwitchPanCancerCount.df[,1] >= 15)
RecurrDELnc.v <- rownames(LncDESwitchPanCancerCount.df)[tmp.idx]

load(file="LncPanCancerCount.Rdata")
tmp.idx <- which(LncPanCancerCount.df[,1] >= 15)
PanCancerLnc.v <- rownames(LncPanCancerCount.df)[tmp.idx]

intersect(RecurrDELnc.v , PanCancerLnc.v )

