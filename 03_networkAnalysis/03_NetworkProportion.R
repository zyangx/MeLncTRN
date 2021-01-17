options(stringsAsFactors=FALSE)
library(igraph)
library(ggsci)
library(scales)
library(ggpubr)
library(stringr) 


cancers.v  <-  c("BLCA", "BRCA(ER+)", "BRCA(ER-)", "CESC", "CHOL", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC",
				"LUAD", "LUSC", "PAAD" , "PRAD", "READ", "THCA", "UCEC")


LncPCGInterAct.v <- c()
for(c in 1:length(cancers.v)) {
	print(cancers.v[c])
	#load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
	load(paste("./InteractNetworkStat1/", cancers.v[c], "_InteractNetStat.Rdata", sep="") )
	tmp.v <- paste(validLncTargetStat.df[,1], validLncTargetStat.df[,2], sep="_")
        LncPCGInterAct.v <- c(LncPCGInterAct.v, tmp.v)
}
LncPCGInterAct.v <- unique(LncPCGInterAct.v)

LncInterAct.v <-  unique(str_split(LncPCGInterAct.v,'[_]',simplify = T)[,1])
PCGInterAct.v <-  unique(str_split(LncPCGInterAct.v,'[_]',simplify = T)[,2])

length(LncInterAct.v)
length(PCGInterAct.v)
length(LncPCGInterAct.v)


###############
ValidGeneProp.m <- matrix(NA, nrow=length(cancers.v), ncol=1)
rownames(ValidGeneProp.m)  <-  cancers.v
colnames(ValidGeneProp.m)[1]  <-  c("Num")

ValidLncProp.m <- matrix(NA, nrow=length(cancers.v), ncol=1)
rownames(ValidLncProp.m)  <-  cancers.v
colnames(ValidLncProp.m)[1]  <-  c("Num")

CancerSpecificStat.m <- matrix(NA, nrow=18, ncol=3)
rownames(CancerSpecificStat.m) <- cancers.v
colnames(CancerSpecificStat.m) <- c("Lnc", "Gene", "Interact")


for(c in 1:length(cancers.v)) {
print(cancers.v[c])

load(paste("./MultiOmicsData/", cancers.v[c], "_MultiOmicsProcess.Rdata", sep=""))
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
load(paste("./InteractNetworkStat1/", cancers.v[c], "_InteractNetStat.Rdata", sep="") )

ValidGeneProp.m[c,1] <- round( ncol(InteractLG.m)/nrow(statLMG.lm[[1]]) , digits=2 )
ValidLncProp.m[c,1] <- round( nrow(InteractLG.m)/length(statLMG.lm) , digits=2 )


CancerSpecificStat.m[c,1] <- round( nrow(InteractLG.m)/length(LncInterAct.v), digits=2 )
CancerSpecificStat.m[c,2] <- round( ncol(InteractLG.m)/length(PCGInterAct.v), digits=2 )
CancerSpecificStat.m[c,3] <- round( nrow(validLncTargetStat.df)/length(LncPCGInterAct.v), digits=2 )

}



mycolors.v <- pal_d3("category20")(18)


library(ggplot2)
library(ggthemes)
ValidGeneProp.df <- as.data.frame(ValidGeneProp.m)
ValidGeneProp.df$Cancer <- rownames(ValidGeneProp.df )

p <- ggbarplot(ValidGeneProp.df, x = "Cancer", y = "Num", fill=mycolors.v, color=mycolors.v,width=0.85, xlab="")+ theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1), text=element_text(family="Arial"), legend.position='none')
p <- p + scale_fill_manual(values =mycolors.v ) +xlab(NULL)
p 
ggsave("DNAmDrivenGenes_Prop_InNetwork1.pdf", device=cairo_pdf, width=6,height=4)


p <- ggplot(ValidGeneProp.df, aes(x = Cancer, y = Num, fill = Cancer)) + geom_bar(stat ="identity") + theme_bw() + theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1), text=element_text(family="Arial"), legend.position='none')
p <- p + scale_fill_manual(values =mycolors.v ) +xlab(NULL)
p 
ggsave("DNAmDrivenGenes_Prop_InNetwork1.pdf", device=cairo_pdf, width=6,height=4)

########
ValidLncProp.df <- as.data.frame(ValidLncProp.m)
ValidLncProp.df$Cancer <- rownames(ValidLncProp.df )

p <- ggplot(ValidLncProp.df, aes(x = Cancer, y = Num, fill = Cancer)) + geom_bar(stat ="identity") + theme_bw() + theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1), text=element_text(family="Arial"), legend.position='none')
p <- p + scale_fill_manual(values =mycolors.v ) +xlab(NULL)
p 
ggsave("Lnc_Prop_ExpressedLnc1.pdf", device=cairo_pdf, width=6,height=4)


###############
###############
Prop.df <- as.data.frame(CancerSpecificStat.m[,1])
Prop.df$Cancer <- rownames(Prop.df )
colnames(Prop.df) <- c("Num", "Cancer");

p <- ggplot(Prop.df, aes(x = Cancer, y = Num, fill = Cancer)) + geom_bar(stat ="identity") + theme_bw() + theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1), text=element_text(family="Arial"), legend.position='none')
p <- p + scale_fill_manual(values =mycolors.v ) +xlab(NULL)
p 
ggsave("Lnc_Prop_InNetwork.pdf", device=cairo_pdf, width=6,height=4)




Prop.df <- as.data.frame(CancerSpecificStat.m[,2])
Prop.df$Cancer <- rownames(Prop.df )
colnames(Prop.df) <- c("Num", "Cancer");

p <- ggplot(Prop.df, aes(x = Cancer, y = Num, fill = Cancer)) + geom_bar(stat ="identity") + theme_bw() + theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1), text=element_text(family="Arial"), legend.position='none')
p <- p + scale_fill_manual(values =mycolors.v ) +xlab(NULL)
p 
ggsave("Gene_Prop_InNetwork.pdf", device=cairo_pdf, width=6,height=4)



Prop.df <- as.data.frame(CancerSpecificStat.m[,3])
Prop.df$Cancer <- rownames(Prop.df )
colnames(Prop.df) <- c("Num", "Cancer");

p <- ggplot(Prop.df, aes(x = Cancer, y = Num, fill = Cancer)) + geom_bar(stat ="identity") + theme_bw() + theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1), text=element_text(family="Arial"), legend.position='none')
p <- p + scale_fill_manual(values =mycolors.v ) +xlab(NULL)
p 
ggsave("LncGeneAct_Prop_InNetwork.pdf", device=cairo_pdf, width=6,height=4)











