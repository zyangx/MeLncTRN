library(igraph)
library(ggsci)
library(scales)
library(ggpubr)

cancers.v  <-  c("BLCA", "BRCA(ER+)", "BRCA(ER-)", "CESC", "CHOL", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC",
				"LUAD", "LUSC", "PAAD" , "PRAD", "READ", "THCA", "UCEC")
mycolors.v <- pal_d3("category20")(18)

pdf("Degree_distribution_of_networks.pdf", width=10, height=13)
par(mfrow=c(c(5,4)))
par(mar=c(4,4,3,1))
for(c in 1:length(cancers.v)) {
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )

igraph.o <- graph_from_incidence_matrix(InteractLG.m)

d <- degree(igraph.o, mode="in")
dd <- degree.distribution(igraph.o, mode="in", cumulative=TRUE)
alpha <- power.law.fit(d, xmin=20)$alpha
#plot(dd, log="xy", xlab="degree", ylab="cumulative frequency",  col=1, main="Nonlinear preferential attachment")

popp<-data.frame(table(d))
#plot(popp[,1],popp[,2],xlab="In-degree",ylab="Frequency",type = "p", col = "black", lwd=2,main = "")

powerfit <- lm(log10(popp[,2])~log10( as.numeric(levels(popp[,1])[popp[,1]] )) )
#summary(powerfit)
rsquare=parse(text=paste('R^2==',round(summary(powerfit)$r.sq,2),sep=''));
fitfun =parse(text=paste('Y==X^-',format(alpha,digits=3),sep=''));
#plot( log( as.numeric(levels(popp[,1])[popp[,1]] )), log(popp[,2]), main = "Cancer", xlab="Degree of Genes",ylab="Frequency", type = "p", col = "skyblue",pch=18, xaxt="n", yaxt="n", cex=0.4, cex.main=0.8, cex.lab=0.9)

plot( log10( as.numeric(levels(popp[,1])[popp[,1]] )), log10(popp[,2]), main = cancers.v[c], xlab="Degree of Genes",ylab="Frequency", type = "p", col = mycolors.v[c], pch=18, xaxt="n", yaxt="n", cex=1.2, cex.main=1.2, cex.lab=1.1)
axis(side=1,at = c(0,1,2,3,4), labels = c("1", "10", "100", "1k", "10k"), tick = TRUE, line = NA, cex=0.5)
axis(side=2,at = c(0,1,2,3,4), labels = c("1", "10", "100", "1k", "10k"), tick = TRUE, line = NA, cex=0.5, las=2)

abline(powerfit, col = "grey",lwd=1.5)
text(x=max(log10( as.numeric(levels(popp[,1])[popp[,1]] )))*0.8, y=max(log10(popp[,2]))*0.9,rsquare);
text(x=max(log10( as.numeric(levels(popp[,1])[popp[,1]] )))*0.8, y=max(log10(popp[,2]))*0.8, fitfun)

}
dev.off()

########################
########################

pcom.lv <- list()
for(c in 1:length(cancers.v)) {
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )

igraph.o <- graph_from_incidence_matrix(InteractLG.m)

d <- degree(igraph.o, mode="in")
dd <- degree.distribution(igraph.o, mode="in", cumulative=TRUE)
alpha <- power.law.fit(d, xmin=20)$alpha
#plot(dd, log="xy", xlab="degree", ylab="cumulative frequency",  col=1, main="Nonlinear preferential attachment")

popp<- as.data.frame(table(d))
#plot(popp[,1],popp[,2],xlab="In-degree",ylab="Frequency",type = "p", col = "black", lwd=2,main = "")
popp[,1] <- as.numeric(as.vector(popp[,1]))
popp[,2] <- as.numeric(as.vector(popp[,2]))

#powerfit <- lm(log10(popp[,2])~log10( as.numeric(levels(popp[,1])[popp[,1]] ) ) )
powerfit <- lm(log10(popp[,2])~log10(popp[,1]) )
intercept <- summary(powerfit)$coeff[1,1]
slope <- summary(powerfit)$coeff[2,1]
#summary(powerfit)
rsquare=parse(text=paste('R^2==',round(summary(powerfit)$r.sq,2),sep=''));
fitfun =parse(text=paste('Y==X^-',format(alpha,digits=3),sep=''));
### + geom_abline(intercept=intercept, slope=slope)  ## + stat_smooth(method=lm, se=FALSE, colour=mycolors.v[c])  ### ylim=c(1, max(popp[,2])), add = "reg.line", 
p <- ggscatter(popp, x="d", y="Freq", main = cancers.v[c], xlab="Degree of Genes",ylab="Frequency",col = mycolors.v[c], ylim=c(1, max(popp[,2])), add = "reg.line", size = 0.8, ggtheme = theme_bw()) + theme(plot.title = element_text(hjust = 0.5) ) + scale_x_log10() + scale_y_log10()   
p <- p + annotate("text", x=max(popp[,1])*0.35, y=max(popp[,2])*0.7, label=rsquare) 
p <- p + annotate("text", x=max(popp[,1])*0.35, y=max(popp[,2])*0.3, label=fitfun) 
p
pcom.lv[[c]] <- p
}

ggarrange(pcom.lv[[1]], pcom.lv[[2]], pcom.lv[[3]], pcom.lv[[4]], pcom.lv[[5]], pcom.lv[[6]], pcom.lv[[7]], pcom.lv[[8]], pcom.lv[[9]], pcom.lv[[10]], pcom.lv[[11]], pcom.lv[[12]], pcom.lv[[13]], pcom.lv[[14]], pcom.lv[[15]], pcom.lv[[16]], pcom.lv[[17]], pcom.lv[[18]], nrow=5, ncol=4, legend="bottom", common.legend=TRUE)
ggsave("Degree_distribution_of_networks2.pdf", width=10,height=13);

c <- 2 
#p <- ggscatter(popp, x="d", y="Freq", xlab="Degree of Genes",ylab="Frequency",col = mycolors.v[c], ylim=c(1, max(popp[,2])), add = "reg.line", size = 0.8, ggtheme = theme_bw()) + theme(plot.title = element_text(hjust = 0.5) ) + scale_x_log10() + scale_y_log10()   
p <- ggscatter(popp, x="d", y="Freq", xlab="Degree of Genes",ylab="Frequency",col = mycolors.v[c], ylim=c(1, max(popp[,2])),  size = 0.8, ggtheme = theme_bw()) + theme(plot.title = element_text(hjust = 0.5) ) + scale_x_log10() + scale_y_log10()   + stat_smooth(method=lm, se=FALSE, colour=mycolors.v[c])
p <- p + annotate("text", x=max(popp[,1])*0.35, y=max(popp[,2])*0.7, label=rsquare) 
p <- p + annotate("text", x=max(popp[,1])*0.35, y=max(popp[,2])*0.3, label=fitfun) 
p

ggsave("Degree_distribution_of_networks_BRCA.pdf", width=3.5,height=3.5);


########################
########################
pdf("LncDegree_distribution_of_networks.pdf", width=10, height=13)
par(mfrow=c(c(5,4)))
par(mar=c(4,4,3,1))
for(c in 1:length(cancers.v)) {
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
LncDegrees.v <- rowSums(InteractLG.m)
popp <- data.frame(table(LncDegrees.v))

powerfit <- lm(log10(popp[,2])~log10( as.numeric(levels(popp[,1])[popp[,1]] )) )
#summary(powerfit)
rsquare <- parse(text=paste('R^2==',round(summary(powerfit)$r.sq,2), sep=''));
pvalue <- paste('P = ', format(summary(powerfit)$coeff[2,4], digits=3), sep='');

hist(LncDegrees.v, breaks=30, las=1, main = cancers.v[c], xlab="", ylab="", col =  mycolors.v[c], cex=1.2, cex.main=1.2, cex.lab=1.1)
text(x=max(as.numeric(levels(popp[,1])[popp[,1]] ))*0.7, y=900,rsquare);

text(x=max(as.numeric(levels(popp[,1])[popp[,1]] ))*0.7, y=750, pvalue)


}
dev.off()


########################
pdf("GeneDegree_distribution_of_networks.pdf", width=10, height=13)
par(mfrow=c(c(5,4)))
par(mar=c(4,4,3,1))
for(c in 1:length(cancers.v)) {
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )
GeneDegrees.v <- colSums(InteractLG.m)
popp <- data.frame(table(GeneDegrees.v))

powerfit <- lm(log10(popp[,2])~log10( as.numeric(levels(popp[,1])[popp[,1]] )) )
#summary(powerfit)
rsquare <- parse(text=paste('R^2==',round(summary(powerfit)$r.sq,2), sep=''));
pvalue <- paste('P = ', format(summary(powerfit)$coeff[2,4], digits=3), sep='');

hist(GeneDegrees.v, breaks=15, las=1, main = cancers.v[c], xlab="", ylab="", col = mycolors.v[c], cex=1.2, cex.main=1.2, cex.lab=1.1)
text(x=max(as.numeric(levels(popp[,1])[popp[,1]] ))*0.8, y=max(popp[,2])*2.5,rsquare);

text(x=max(as.numeric(levels(popp[,1])[popp[,1]] ))*0.8, y=max(popp[,2])*2, pvalue)


}
dev.off()

