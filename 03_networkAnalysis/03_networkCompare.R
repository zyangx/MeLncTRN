Jaccard  <- function(A,B){
	Dist <- length(intersect(A,B)) /length(union(A,B))
}

get_lower_tri <- function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
 }

cancers.v  <-  c("BLCA", "BRCA(ER+)", "BRCA(ER-)", "CESC", "CHOL", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC",
				"LUAD", "LUSC", "PAAD" , "PRAD", "READ", "THCA", "UCEC")

Dist.m <- matrix(NA, nrow=18, ncol=18)
rownames(Dist.m) <- cancers.v
colnames(Dist.m) <- cancers.v


for(i in 1:18){
    load(paste("./InteractNetwork1/", cancers.v[i], "_InteractNet.Rdata", sep="") )
    A.v <- rownames(InteractLG.m)
    for(j in 1:18){
        load(paste("./InteractNetwork1/", cancers.v[j], "_InteractNet.Rdata", sep="") )
	B.v <- rownames(InteractLG.m)
	dist <- Jaccard(A.v, B.v)
	Dist.m[i,j] <- dist

    }
}

Dist.m <- round(Dist.m,2)

col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white", "cyan", "#007FFF", "blue","#00007F"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))
col3 <- colorRampPalette(c("red", "white", "blue"))
col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F", "cyan", "#007FFF", "blue", "#00007F"))
col5 <- colorRampPalette(c( "green","green", "white", "red","red"))

col = col1(20)

library(corrplot)
corrplot(Dist.m, method=c("square"), is.corr=T, type = "upper", col = col5(200), cl.lim=c(0,1))
corrplot.mixed(Dist.m, tl.pos = "lt", diag = "l",  is.corr=F,cl.lim=c(0,1))


library(corrgram)
corrgram(Dist.m, order=NULL, lower.panel=panel.shade, upper.panel=panel.pie, text.panel=panel.txt)
corrgram(Dist.m, order=NULL, lower.panel = panel.cor, upper.panel = panel.pie, col.regions = colorRampPalette(c( "blue","blue", "white", "red")))

library(ggplot2)
library(reshape2)
upper_tri <- get_upper_tri(Dist.m)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

lower_tri <- get_lower_tri(Dist.m)
melted_cormat <- melt(lower_tri, na.rm = TRUE)

####
melted_cormat  -> lncRNA.m

ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "#1B9E77", high = "#D95F02", mid = "white", 
   midpoint = 0.5, limit = c(0,1), space = "Lab", 
    name="Jaccard") +
  theme_minimal()+ # minimal theme
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 8, hjust = 1), axis.text.y = element_text( size = 8))+
 coord_fixed() +  scale_y_discrete(position = "right")

# Print the heatmap
print(ggheatmap)
ggheatmap <- ggheatmap +
geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  legend.justification = c(1, 0),
  legend.position = c(0.6, 0.7),
  legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                title.position = "top", title.hjust = 0.5)) 
ggheatmap <- ggheatmap +  scale_y_discrete(position = "right")

print(ggheatmap)

ggsave("Lnc_Jaccard_heatplot.pdf", width=5,height=5)



##################33


Dist.m <- matrix(NA, nrow=18, ncol=18)
rownames(Dist.m) <- cancers.v
colnames(Dist.m) <- cancers.v


for(i in 1:18){
    load(paste("./InteractNetwork1/", cancers.v[i], "_InteractNet.Rdata", sep="") )
    A.v <- colnames(InteractLG.m)
    for(j in 1:18){
        load(paste("./InteractNetwork1/", cancers.v[j], "_InteractNet.Rdata", sep="") )
	B.v <- colnames(InteractLG.m)
	dist <- Jaccard(A.v, B.v)
	Dist.m[i,j] <- dist

    }
}

Dist.m <- round(Dist.m,2)

upper_tri <- get_upper_tri(Dist.m)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

####
melted_cormat  -> mRNA.m

ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "#1B9E77", high = "#D95F02", mid = "white", 
   midpoint = 0.5, limit = c(0,1), space = "Lab", 
    name="Jaccard") +
  theme_minimal()+ # minimal theme
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 8, hjust = 1), axis.text.y = element_text( size = 8))+
 coord_fixed() +  scale_y_discrete(position = "right")

# Print the heatmap
print(ggheatmap)
ggheatmap <- ggheatmap +
geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  legend.justification = c(1, 0),
  legend.position = c(0.6, 0.7),
  legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                title.position = "top", title.hjust = 0.5)) 
ggheatmap <- ggheatmap +  scale_y_discrete(position = "right")

print(ggheatmap)

ggsave("Gene_Jaccard_heatplot.pdf", width=5,height=5)


##################33


Dist.m <- matrix(NA, nrow=18, ncol=18)
rownames(Dist.m) <- cancers.v
colnames(Dist.m) <- cancers.v


for(i in 1:18){
    load(paste("./InteractNetwork1/", cancers.v[i], "_InteractNet.Rdata", sep="") )
    tmp.df <- melt(InteractLG.m)
    tmp.idx <- which(tmp.df[,3] == 1)
    tmp.df <- tmp.df[tmp.idx,]
    A.v <- paste0(tmp.df[,1], "_", tmp.df[,2])

    for(j in 1:18){
        load(paste("./InteractNetwork1/", cancers.v[j], "_InteractNet.Rdata", sep="") )
        tmp.df <- melt(InteractLG.m)
        tmp.idx <- which(tmp.df[,3] == 1)
        tmp.df <- tmp.df[tmp.idx,]
        B.v <- paste0(tmp.df[,1], "_", tmp.df[,2])
	dist <- Jaccard(A.v, B.v)
	Dist.m[i,j] <- dist

    }
}

Dist.m <- round(Dist.m,2)


upper_tri <- get_upper_tri(Dist.m)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

####
melted_cormat  -> pair.m

ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "#1B9E77", high = "#D95F02", mid = "white", 
   midpoint = 0.5, limit = c(0,1), space = "Lab", 
    name="Jaccard") +
  theme_minimal()+ # minimal theme
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 8, hjust = 1), axis.text.y = element_text( size = 8))+
 coord_fixed() +  scale_y_discrete(position = "right")

# Print the heatmap
print(ggheatmap)
ggheatmap <- ggheatmap +
geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  legend.justification = c(1, 0),
  legend.position = c(0.6, 0.7),
  legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                title.position = "top", title.hjust = 0.5)) 
ggheatmap <- ggheatmap +  scale_y_discrete(position = "right")

print(ggheatmap)

ggsave("GeneLncInteract_Jaccard_heatplot.pdf", width=5,height=5)

######################

data.m <- matrix(NA, nrow=171, ncol=3)
colnames(data.m) <- c("LncRNA" ,"PCG","Pair")
data.m[,1] <- as.numeric(lncRNA.m[,3])
data.m[,2] <- as.numeric(mRNA.m[,3])
data.m[,3] <- as.numeric(pair.m[,3])
tmp.idx <- which(data.m[,1] == 1)
data.m <- data.m[-tmp.idx,]


data.df <- as.data.frame(data.m)
require(ggpubr)
library(reshape2)
library(ggsci)
mydata <- melt(data.df)
colnames(mydata) <- c("Class", "Jaccard_Coefficient")
p <- ggboxplot(mydata, x = "Class", y = "Jaccard_Coefficient", color = "Class", add = "jitter",   line.color = "gray")
p <- p+stat_compare_means(comparisons = list(c("LncRNA", "PCG"), c("LncRNA", "Pair")),label= "p.format")
p

ggsave("Group_Jaccard_boxplot.pdf", width=5,height=5)


