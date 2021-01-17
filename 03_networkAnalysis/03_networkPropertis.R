library(igraph)
library(ggsci)
library(scales)

cancers.v  <-  c("BLCA", "BRCA(ER+)", "BRCA(ER-)", "CESC", "CHOL", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC",
				"LUAD", "LUSC", "PAAD" , "PRAD", "READ", "THCA", "UCEC")
mycolors.v <- pal_d3("category20")(18)

networkStat.m <- matrix(NA, nrow=18, ncol=7)
rownames(networkStat.m) <- cancers.v;
colnames(networkStat.m) <- c("LncNum", "GeneNum", "LncDegree","GeneDegree", "EdgeNum", "AveDegree", "Rsquare")


for(c in 1:length(cancers.v)) {
print(cancers.v[c])
load(paste("./InteractNetwork1/", cancers.v[c], "_InteractNet.Rdata", sep="") )

networkStat.m[c, 1] <- nrow(InteractLG.m)
networkStat.m[c, 2] <- ncol(InteractLG.m)
networkStat.m[c, 3] <- round(mean(rowSums(InteractLG.m)),1)
networkStat.m[c, 4] <- round(mean(colSums(InteractLG.m)),1)

igraph.o <- graph_from_incidence_matrix(InteractLG.m)
num.edges <- length(E(igraph.o))
num.vertices <- length(V(igraph.o))

connectance <- edge_density(igraph.o,loops=FALSE)

average.degree <- mean(igraph::degree(igraph.o))
average.degree

average.path.length <-  average.path.length(igraph.o) # equals to: mean_distance(igraph)
average.path.length

diameter <- diameter(igraph.o, directed = FALSE, unconnected = TRUE, weights = NULL)
diameter

edge.connectivity <- edge_connectivity(igraph.o)
edge.connectivity

clustering.coefficient <- transitivity(igraph.o) 
clustering.coefficient

transitivity(igraph.o, type = "average")

no.clusters <- no.clusters(igraph.o)
no.clusters

centralization.betweenness <-  centralization.betweenness(igraph.o)$centralization 
centralization.betweenness

centralization.degree <-  centralization.degree(igraph.o)$centralization
centralization.degree

networkStat.m[c, 5] <- num.edges
networkStat.m[c, 6] <- round(average.degree,1)

d <- degree(igraph.o, mode="in")
dd <- degree.distribution(igraph.o, mode="in", cumulative=TRUE)
#networkStat.m[c, 7] <- round(power.law.fit(d, xmin=20)$alpha,2)

popp<-data.frame(table(d))
#plot(popp[,1],popp[,2],xlab="In-degree",ylab="Frequency",type = "p", col = "black", lwd=2,main = "")

powerfit <- lm(log10(popp[,2])~log10( as.numeric(levels(popp[,1])[popp[,1]] )) )
#summary(powerfit)
networkStat.m[c, 7] <- round(summary(powerfit)$r.sq,2) 

}
write.table(networkStat.m, file="networkStat.txt", quote=F, sep="\t", col.names=NA)


##############
tmp.idx <- 18:1
networkStat.m <- networkStat.m[tmp.idx,]


library(marray)
library(plotrix);
library(RColorBrewer)


draw_image <- function(data, colorgroup=1, axis1 = FALSE, axis2 = FALSE, label = FALSE) {
  # set parameters
  breaks.frequency <- seq(from=min(data), to=max(data), length.out=10)
  if(colorgroup == 1){
    myColors <- colorRampPalette(c("orange","red"))
  }else if(colorgroup == 2){
    myColors <- colorRampPalette(c("yellowgreen", "yellow"))
  }else if(colorgroup == 3){
    myColors <- colorRampPalette(c("Violet","purple", "purple"))
  }else if(colorgroup == 4){
    myColors <- colorRampPalette(c("skyblue","blue"))
  }else if(colorgroup == 5){
    myColors <- colorRampPalette(c("bisque","chocolate"))
  }
  # generate figure
  image(1:nrow(data), 1:ncol(data), as.matrix(data),  breaks=breaks.frequency, col=myColors(length(breaks.frequency)-1), axes = FALSE, cex = 0.8, xlab = "", ylab = "")

  # draw axis
  if (axis1) {
    axis(2,at=1:nrow(tmp.m),labels=rownames(tmp.m),las=2);
    #axis(2, 1:ncol(data), colnames(data), cex.axis=2.5)
    #axis(1, 1:nrow(data), rownames(mxdata), cex.axis=2.5)
  }

  if (axis2) {
    #text(x=1:ncol(tmp.m), y =-0.1, labels=colnames(tmp.m), srt = 90,xpd=T);
    axis(1,at=1:ncol(tmp.m),labels=colnames(tmp.m),las=2);
    #axis(2, 1:ncol(data), colnames(data), cex.axis=2.5)
    #axis(1, 1:nrow(data), rownames(mxdata), cex.axis=2.5)
  }

 #  draw numbers
  if (label) {
    for (x in 1:nrow(data)) {
      for (y in 1:ncol(data))  {
        text(x, y, data[x, y], cex = 0.9)
        }
    }
  }

}


pdf("Network_Stat_1adj.pdf", width=3.5, height=7)
layout(matrix(c(1,1,1,1,1,1,1,1,1,  2, 2, 2,   3, 3, 3,   4, 4, 4,  5, 5, 5,  6, 6, 6,  7, 7, 7,7), nr=1, byrow=TRUE))
#layout(matrix(c(1,1,1,  2,   3,    4,  5,  6 ,7), nr=1, byrow=TRUE))
#par(mfrow=c(1,7))

par(mar=c(6, 5.8, 1.5, 0.07))
tmp.m <- as.matrix(networkStat.m[,1])
colnames(tmp.m) <- colnames(networkStat.m)[1]
draw_image(t(tmp.m),1, TRUE,TRUE,TRUE)

par(mar=c(6, 0.07, 1.5, 0.07))
tmp.m <- as.matrix(networkStat.m[,2])
colnames(tmp.m) <- colnames(networkStat.m)[2]
draw_image(t(tmp.m),1, FALSE,TRUE,TRUE)

tmp.m <- as.matrix(networkStat.m[,3])
colnames(tmp.m) <- colnames(networkStat.m)[3]
draw_image(t(tmp.m),2, FALSE,TRUE,TRUE)

tmp.m <- as.matrix(networkStat.m[,4])
colnames(tmp.m) <- colnames(networkStat.m)[4]
draw_image(t(tmp.m),2, FALSE,TRUE,TRUE)

tmp.m <- as.matrix(networkStat.m[,5])
colnames(tmp.m) <- colnames(networkStat.m)[5]
draw_image(t(tmp.m),3, FALSE,TRUE,TRUE)

tmp.m <- as.matrix(networkStat.m[,6])
colnames(tmp.m) <- colnames(networkStat.m)[6]
draw_image(t(tmp.m),4, FALSE,TRUE,TRUE)

par(mar=c(6, 0.05, 1.5, 0.8))
tmp.m <- as.matrix(networkStat.m[,7])
colnames(tmp.m) <- colnames(networkStat.m)[7]
draw_image(t(tmp.m),5, FALSE,TRUE,TRUE)

dev.off()

