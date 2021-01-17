options(stringsAsFactors=FALSE)

library(reshape2)
library(ggpubr)
library(RColorBrewer)
library(ggsci)
library(scales)
mycolors.v <- pal_d3("category20")(18)

cancers.v  <-  c("BLCA", "BRCA(ER+)", "BRCA(ER-)", "CESC", "CHOL", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC",
				"LUAD", "LUSC", "PAAD" , "PRAD", "READ", "THCA", "UCEC")


dataC.df <- data.frame(Cancer=0, Type=0, Train=0, Test=0)
dataC.df <- dataC.df[-1,]

for(c in 1:length(cancers.v)) {
#rm(list=ls())
print(cancers.v[c])
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
load(paste("./Survival_Lnc_Gene/", cancers.v[c], "_Survival.Rdata", sep="") )

tmp.idx <- match(rownames(InteractLG.m), rownames(Lnc_cox_results_Train.m))
marker.m <- cbind(rep(cancers.v[c], nrow(InteractLG.m)),  rep("Lnc", nrow(InteractLG.m)) )
pvalues.m <- cbind(Lnc_cox_results_Train.m[tmp.idx,6], Lnc_cox_results_Test.m[tmp.idx,6])
LncPvalue.m <- cbind(marker.m, pvalues.m)
#plot(-log10(Lnc_cox_results_Train.m[tmp.idx,6]),  -log10(Lnc_cox_results_Test.m[tmp.idx,6]) )

tmp.idx <- match(colnames(InteractLG.m), rownames(PCG_cox_results_Train.m))
marker.m <- cbind(rep(cancers.v[c], ncol(InteractLG.m)),  rep("Gene", ncol(InteractLG.m)) )
pvalues.m <- cbind(PCG_cox_results_Train.m[tmp.idx,6], PCG_cox_results_Test.m[tmp.idx,6])
PCGPvalue.m <- cbind(marker.m, pvalues.m)

allPvalue.m <- rbind(LncPvalue.m, PCGPvalue.m )

dataC.df <- rbind(dataC.df, allPvalue.m)

}

colnames(dataC.df) <- c("Cancer", "Type", "Train", "Test")
dataC.df[,3] <- -log10(as.numeric(dataC.df[,3]))
dataC.df[,4] <- -log10(as.numeric(dataC.df[,4]))

sig.idx <- intersect(which(dataC.df[,3] > -log10(0.05)), which(dataC.df[,4] > -log10(0.05)))
length(sig.idx) / nrow(dataC.df)

#plot(-log10(as.numeric(dataC.df[,3])),  -log10(as.numeric(dataC.df[,4])) , ylim=c(1,50), xlim=c(1,50) )

#data.df <- melt(dataC.df, id.vars=c("Cancer", "Type"), variable.name="Class", value.name="Value")

p <- ggscatter(dataC.df, x="Train", y="Test", color="Type", ylim=c(0,20), xlim=c(0,20), ylab="-log10(P)_Validation", xlab="-log10(P)_Discovery", palette = brewer.pal(2,"Set1")[1:2], cor.coef = TRUE,  cor.coeff.args = list(method = "pearson", label.x = 10,  label.y = 19, label.sep = ", ", size=5) )
#p <- p +stat_cor(method = "pearson", label.x = 20, label.y = 35, size=7) 
p <- p + theme_bw() +theme(legend.position="top", legend.text=element_text(size=15))
p <- p +  geom_line(aes(x=-log10(0.05)), linetype=2, color="black")  + geom_line(aes(y=-log10(0.05)), linetype=2, color="black")
p <- p + annotate("text", x=11, y=13, label="5.32%",size=7)
p

ggsave("Lnc_Gene_Survival_P_corScatter.pdf", width=6, height=6)
ggsave("Lnc_Gene_Survival_P_corScatter.png", width=6, height=6)



#################################
dataP.df <- as.data.frame(matrix(NA, nrow=18, ncol=5))
colnames(dataP.df) <- c("Cancer","TrainL", "TrainP",  "TestL", "TestP")
rownames(dataP.df) <- cancers.v
dataP.df[,1] <- cancers.v

for(c in 1:length(cancers.v)) {
#rm(list=ls())
print(cancers.v[c])

load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
load(paste("./Survival_Lnc_Gene/", cancers.v[c], "_Survival.Rdata", sep="") )

tmp.idx <- match(rownames(InteractLG.m), rownames(Lnc_cox_results_Train.m))
dataP.df[c,2]  <- length(which(Lnc_cox_results_Train.m[tmp.idx,6] < 0.05)) / nrow(InteractLG.m)

tmp.idx <- match(colnames(InteractLG.m), rownames(PCG_cox_results_Train.m))
dataP.df[c,3]  <- length(which(PCG_cox_results_Train.m[tmp.idx,6] < 0.05)) / ncol(InteractLG.m)

tmp.idx <- match(rownames(InteractLG.m), rownames(Lnc_cox_results_Test.m))
dataP.df[c,4]  <- length(which(Lnc_cox_results_Test.m[tmp.idx,6] < 0.05)) / nrow(InteractLG.m)

tmp.idx <- match(colnames(InteractLG.m), rownames(PCG_cox_results_Test.m))
dataP.df[c,5]  <- length(which(PCG_cox_results_Test.m[tmp.idx,6] < 0.05)) / ncol(InteractLG.m)

}

dataTrain.df <- dataP.df[,c(1,2,3)]
colnames(dataTrain.df) <- c("Cancer", "Lnc",  "PCG")
dataTrain.df <- melt(dataTrain.df)
names(dataTrain.df) <- c("Cancer", "Type",  "Percent")
p1 <- ggplot(dataTrain.df, aes(x=Type, y=Percent, group=Cancer)) + scale_colour_manual(values =brewer.pal(2,"Set1")[2:1])
p1 <- p1 + geom_line( )  + geom_point(size=3, shape=20, aes(colour=Type)) + theme_bw() + scale_fill_manual() +theme(legend.position="none")
p1
ggsave("Lnc_Gene_Survival_PerecentScatter_Train.pdf", width=4, height=5)


dataTest.df <- dataP.df[,c(1,4,5)]
colnames(dataTest.df) <- c("Cancer", "Lnc",  "PCG")
dataTest.df <- melt(dataTest.df)
names(dataTest.df) <- c("Cancer", "Type",  "Percent")
p1 <- ggplot(dataTest.df, aes(x=Type, y=Percent, group=Cancer)) + scale_colour_manual(values =brewer.pal(2,"Set1")[2:1])
p1 <- p1 + geom_line( )  + geom_point(size=3, shape=20, aes(colour=Type)) + theme_bw() + scale_fill_manual() +theme(legend.position="none")
p1
ggsave("Lnc_Gene_Survival_PerecentScatter_Test.pdf", width=4, height=5)

#################################

load("LncPanCancerCount.Rdata")


dataH.df <- data.frame(HazardRatio=0, Type=0, Cancer=0)
dataH.df <- dataH.df[-1,]

for(c in 1:length(cancers.v)) {
print(cancers.v[c])
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
load(paste("./Survival_Lnc_Gene/", cancers.v[c], "_Survival.Rdata", sep="") )

tmp.idx <- match(rownames(InteractLG.m), rownames(Lnc_cox_results_Train.m))

HazardRatio.v <- Lnc_cox_results_Train.m[tmp.idx,3]
tmp.idx <- match(names(HazardRatio.v),  rownames(LncPanCancerCount.df))

data.df <- as.data.frame(cbind(HazardRatio.v, LncPanCancerCount.df[tmp.idx,3]))
colnames(data.df) <- c("HazardRatio", "Type")
data.df$HazardRatio <- as.numeric(data.df$HazardRatio)
data.df <- within(data.df, Type <- factor(Type, levels = c("Pan-cancer", "Others")))
data.df$Cancer <- cancers.v[c]

dataH.df <- rbind(dataH.df, data.df)

}

p <- ggboxplot(dataH.df, x="Cancer", y="HazardRatio", fill = "Cancer", color="Type", palette = mycolors.v, add = "none", size=0.2, ylab="HR", xlab="Cancer Types", ggtheme = theme_bw()) + theme(plot.title = element_text(hjust = 0.5),  axis.text.x=element_blank())
p <- p+ stat_compare_means(aes(group=Type),label = "p.signif")
p <- p+ scale_colour_manual(values=c("darkred", "black"))
p <- p+ annotate("text", x=3, y=9, label="****  P < 0.0001", size=6) 
#p <- p   + geom_line(aes(x=-log2(0.05)))  + geom_line(aes(y=-log2(0.05)))
p  
ggsave("LncCategory_Compare_by_PanCancer_HazardRatio.pdf", width=12,height=6)


#dataH.df <- within(dataH.df, Type <- factor(Type, levels = c("Others", "Pan-cancer")))
dataH.df <- within(dataH.df, Type <- factor(Type, levels = c("Pan-cancer", "Others")))
p <- ggboxplot(dataH.df, x="Cancer", y="HazardRatio", fill = "Cancer", color="Type", palette = mycolors.v, add = "none", size=0.2, ylab="HR", xlab="Cancer Types", ggtheme = theme_bw()) + theme(plot.title = element_text(hjust = 0.5),  axis.text.y=element_blank() )
p <- p+ stat_compare_means(aes(group=Type),label = "p.signif")
p <- p+ scale_colour_manual(values=c("darkred", "black"))
p <- p+ coord_flip()
p <- p+ annotate("text", x=15, y=9, label="****  P < 0.0001") 
p <- p+ scale_x_discrete(limits=rev(cancers.v))
p
ggsave("LncCategory_Compare_by_PanCancer_HazardRatio1.pdf", width=6,height=10)


#################################
dataH.df <- as.data.frame(matrix(NA, nrow=18, ncol=4))
colnames(dataH.df) <- c("TrainL", "TrainP",  "TestL", "TestP")
rownames(dataH.df) <- cancers.v
dataH.df[,1] <- cancers.v

for(c in 1:length(cancers.v)) {
#rm(list=ls())
print(cancers.v[c])

load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
load(paste("./Survival_Lnc_Gene/", cancers.v[c], "_Survival.Rdata", sep="") )

tmp.idx <- match(rownames(InteractLG.m), rownames(Lnc_cox_results_Train.m))
dataH.df[c,1]  <- mean(Lnc_cox_results_Train.m[tmp.idx,3])

tmp.idx <- match(colnames(InteractLG.m), rownames(PCG_cox_results_Train.m))
dataH.df[c,2]  <- mean(PCG_cox_results_Train.m[tmp.idx,3])

tmp.idx <- match(rownames(InteractLG.m), rownames(Lnc_cox_results_Test.m))
dataH.df[c,3]  <- mean(Lnc_cox_results_Test.m[tmp.idx,3])

tmp.idx <- match(colnames(InteractLG.m), rownames(PCG_cox_results_Test.m))
dataH.df[c,4]  <- mean(PCG_cox_results_Test.m[tmp.idx,3])

}



#################################
dataH.df <- as.data.frame(matrix(NA, nrow=18, ncol=2))
colnames(dataH.df) <- c("Cancer", "Num")
rownames(dataH.df) <- cancers.v
dataH.df[,1] <- cancers.v

#for(c in 1:length(cancers.v)) {
for(c in c(1,2,3,4,6,7,8,9,10,11,13,14,15,16,17,18)) {
#rm(list=ls())
print(cancers.v[c])

load(paste("./Biclique/biclique_MBEA_R/Survival_Results/", cancers.v[c], "_module_cox.Rdata", sep="") )

#dataH.df[c,2]  <- length(which(Module_Cox_Results.m[,6] < 0.05)) / nrow(Module_Cox_Results.m)

#dataH.df[c,3]  <-  length(which(Module_Cox_Results.m[,13] < 0.05)) / nrow(Module_Cox_Results.m)

print(length(which(Module_Cox_Results.m[,13] < 0.05)) )
dataH.df[c,2] <- length(which(Module_Cox_Results.m[,13] < 0.05))
}

dataH.df[c(5,12),2] <- 0

dataH.df <- dataH.df[which(dataH.df[,2] > 0),]

dataH.df <- data.frame("Cancer" =c("THCA", "KIRP", "Others"), Num=c(1584, 245, 232))
dataH.df$percent = round(dataH.df$Num/sum(dataH.df$Num) * 100)
lbls <- paste(dataH.df$Cancer, "\n(", dataH.df$Num, ")", sep= "");
ggpie(dataH.df, "Num", label  = lbls, fill = "Cancer", lab.pos = "in",  lab.font = "white", color = "white", palette = c( "#FC4E07","#E7B800","#00AFBB"), orientation="reverse") + theme(legend.position="none")
ggsave("Module_Survival_percentPie1.pdf", width=5, height=5)

