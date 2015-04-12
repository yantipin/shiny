# Code started 04/12/2015 
## If running the code for the 1st time install all packages:
# Un-comment lines 5 - 7 and run the code 

# install.packages("RColorBrewer", "gplots", "devtools")
# library(devtools)
# install_github("ririzarr/rafalib")
library(gplots); library(RColorBrewer); library(rafalib)

# reading in data 
df <- read.table("../data/COAD.HiSeq.16genes.matrix.txt", sep="\t", header = TRUE)

# making df a numeric matrix 
mat <- as.matrix(df[,2:dim(df)[2]])
rownames(mat) <- df$Hybridization.REF

# taking only 50 random samples of samples(columns) (without replacement)
# indeces of columns to be kept (will transpose df)
set.seed(1)
idx <- t(df)[sample(1:nrow(t(df)), 50, replace=FALSE),]
matS <- as.matrix(df[, which(colnames(df) %in% row.names(idx))])
# assigning rownames of Hybridization.REF to keep numeric matrix
rownames(matS) <- df$Hybridization.REF

# Hierarchical Clustering 
d <- dist(t(matS))
hc <- hclust(d)
hc
plot(hc, cex=0.5)
abline(h=5000)
# cutting the tree
hclusters <- cutree(hc, h=5000)
length(table(cluster=hclusters))
# adding colors 
# http://rstudio-pubs-static.s3.amazonaws.com/1876_df0bf890dd54461f98719b461d987c3d.html
labelColors <- colorRampPalette(brewer.pal(9, "RdYlGn"))(length(table(cluster=hclusters)))
# function to get color labels
colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- labelColors[hclusters[which(names(hclusters) == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
  }
  n
}
# using dendrapply
hcd <- as.dendrogram(hc)
clusDendro <- dendrapply(hcd, colLab)
# make plot
par(cex=0.6)
plot(clusDendro)
abline(h=5000, col="red")


# Heatmap 
# image(matS)
cols <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(ncol(matS))
colsR <- colorRampPalette(rev(brewer.pal(11,"Set3")))(nrow(matS))
heatmap.2(matS, col=cols,
          trace="none", scale="row",
          #RowSideColors=colsR,
          key=TRUE)


