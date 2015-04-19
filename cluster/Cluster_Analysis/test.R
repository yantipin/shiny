# Code started 04/12/2015 
## If running the code for the 1st time install all packages:
# Un-comment lines 5 - 7 and run the code 

# install.packages("RColorBrewer", "gplots", "devtools")
# library(devtools)
# install_github("ririzarr/rafalib")
library(gplots); library(RColorBrewer); library(rafalib)

# reading in data 
path <- "../data/COAD.HiSeq.16genes.matrix.txt"
path <- "~/Dropbox/Yevgeniy-to-Maria/clustering-tool/TCGA-COAD-data/COAD.HiSeq.matrix.txt"
df <- read.table(path, sep="\t", header = TRUE)

# making df a numeric matrix 
mat <- as.matrix(df[ 2:dim(df)[1], 2:dim(df)[2] ])

# leave gene names only, remove GeneID
genenames <- df$Hybridization.REF[ 2:dim(df)[1] ]
aa <- lapply(genenames, function(x) { 
   a <- strsplit(x, "|", fixed=TRUE) 
   return(a[1]) 
} )


b <- noquote(strsplit(names[1], "|", fixed=TRUE)[[1]])

b[1]

b <- "123"
class(b) <- "numeric"
b

str <- "?|100130426"
a <- strsplit(str, "|", fixed=TRUE)

q <- 5
q <- 





rownames(mat) <- 

class(mat) <- "numeric"

# -------------- removing genes with low variance -------------
# first find the desired quantile breaks for the entire matrix
# TODO: what we really want is to find variance within each gene individually
# qt <- quantile( mat , probs = c(0.1,0.9) )
# qt
# 20%  80% 
# 5.17 6.62 
# Next get a logical vector of the rows that have any values outside these breaks
# rows <- apply( mat, 1, function(x) any( x < qt[1] | x > qt[2] ) )
#  Subset on this vector
# mat <- mat[ rows , ]

# --- fix samples names ------------------
colnames(mat) <- lapply(colnames(mat), function(x) substr(gsub(".","-",x,fixed=TRUE), 1, 16))

# --- now lets load external list of samples from the table -------
path <- "~/Dropbox/Yevgeniy-to-Maria/clustering-tool/TCGA-COAD-data/COAD_x_QN.txt"
qdf <- read.table(path, sep="\t", header=TRUE)

# --- leave in the matrix only sample present in qdf
smat <- mat[, which(colnames(mat) %in% qdf$tcga_barcode)]

dim(mat)
dim(smat)

# --- remove genes with variace < 4 ---------
rows <- apply(smat, 1, function(row) var(row) > 4 )
fmat <- smat[ rows, ]
fmat <- fmat[1:10, ]

dim(fmat)
tfmat <- t(fmat)
dim(tfmat)

# hierarchical clustering
# cor(tfmat, use="all.obs", method="spearman")

d <- dist(tfmat)
hc <- hclust(d)
hc
plot(hc, cex=0.5, xlim=1000, ylim=500)
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
heatmap.2(fmat, col=cols,
          trace="none", scale="row",
          #RowSideColors=colsR,
          key=TRUE)






