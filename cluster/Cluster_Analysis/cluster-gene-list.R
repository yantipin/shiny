
# Code started 04/12/2015 
## If running the code for the 1st time install all packages:
# Un-comment lines 5 - 7 and run the code 

# install.packages("RColorBrewer", "gplots", "devtools")
# library(devtools)
# install_github("ririzarr/rafalib")
library(gplots); library(RColorBrewer); library(rafalib)

# reading sample clinical data
path <- "~/Dropbox/Yevgeniy-to-Maria/clustering-tool/TCGA-COAD-data/COAD_x_QN.txt"
cli <- read.table(path, sep="\t", header = TRUE)
cli_samples <- as.vector(cli$tcga_barcode)  # sample names present in clinical data
cliQN <- as.vector(cli$label)               # sample QN status
cliMut <- as.vector(cli$mutations)          # sample's number of mutation


# reading in data 
path <- "../data/COAD.HiSeq.16genes.matrix.txt"
path <- "~/Dropbox/Yevgeniy-to-Maria/clustering-tool/TCGA-COAD-data/COAD.HiSeq.matrix.txt"
df <- read.table(path, sep="\t", header = TRUE)

# making df a numeric matrix 
mat <- as.matrix(df[ 2:dim(df)[1], 2:dim(df)[2] ])

# split gene name and GeneID
a <- df$Hybridization.REF[ 2:dim(df)[1] ]
b <- as.character(a)
genenames <- unlist(lapply(strsplit(b, "|", fixed="TRUE"), `[[`, 1))
GeneIDs <- unlist(lapply(strsplit(b, "|", fixed="TRUE"), `[[`, 2))

# read gene list from a file
path <- "~/Dropbox/Yevgeniy-to-Maria/clustering-tool/TCGA-COAD-data/diff-expressed-genes.txt"
df2 <- read.table(path, sep="\t", header = TRUE)
mat <- as.matrix(df[ 2:dim(df)[1], 2:dim(df)[2] ])
GeneIDlist <- df2[1:dim(df2)[1],2]

## b <- noquote(strsplit(names[1], "|", fixed=TRUE)[[1]])
## class(mat) <- "numeric"

# --- fix samples names ------------------
colnames(mat) <- lapply(colnames(mat), function(x) substr(gsub(".","-",x,fixed=TRUE), 1, 16))

# --- assign GeneIDs as row names --------
rownames(mat) <- GeneIDs

# --- now lets load external list of samples from the table -------
path <- "~/Dropbox/Yevgeniy-to-Maria/clustering-tool/TCGA-COAD-data/COAD_x_QN.txt"
qdf <- read.table(path, sep="\t", header=TRUE)

# --- leave in the matrix only sample present in qdf
smat <- mat[, which(colnames(mat) %in% qdf$tcga_barcode)]

# --- remove genes with variace < 4 ---------
rows <- apply(smat, 1, function(row) var(row) > 4 )
fmat <- smat[ rows, ]

# --- leave only gene list loaded earlier ---
fmat2 <- fmat[ which(rownames(fmat) %in% GeneIDlist), ]

class(fmat2) <- "numeric"

# hierarchical clustering
#
# eucledian distance -------
d <- dist(fmat2)
hc <- hclust(d)
hc
plot(hc, cex=0.5, xlim=1000, ylim=500)


# spearman correlations ----
# sp <- cor(fmat, use="all.obs", method="spearman")
# dist <- as.dist(sp)
# hc <- hclust(dist)
# hc
# plot(hc, cex=0.5, xlim=1000, ylim=500)
# abline(h=5000)

# cutting the tree
hclusters <- cutree(hc, h=200000)
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
abline(h=200000, col="red")


# Heatmap 
# image(fmat2)
cols <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(ncol(fmat2))
colsR <- colorRampPalette(rev(brewer.pal(11,"Set3")))(nrow(fmat2))

# convert gene IDs to gene names on the heatmap
rownames_GeneID <- rownames(fmat2)
rownames(fmat2) <- genenames[which(GeneIDs %in% rownames(fmat2))]

# convert samples ID to QN status
colnames_sample <- colnames(fmat2)
colnames(fmat2) <- cliQN[ which(cli_samples %in% colnames(fmat2)) ]

# creates a own color palette from red to green
betterPal <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
# (optional) defines the color breaks manually for a "skewed" color transition
betterBreaks = c(
  seq( 0, 100, length=100),   # for red
  seq( 100, 400, length=100),   # for yellow
  seq( 400, 1000, length=100) )  # for green

heatmap.2(fmat2, trace="none", col=cols, scale="row", key=TRUE)

#          col=betterPal, 
#          breaks=betterBreaks






