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
rownames(mat) <- df$Hybridization.REF[ 2:dim(df)[1] ]

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
dim(fmat)


hist(fmat, breaks=40, col="blue")


d <- density(fmat[1,])
plot(d)
polygon(d, col="red")

var(mat[1,]) < 1

qt <- quantile( mat[1,], probs = c(0.1, 0.9) )
qt
any(mat[1, ] > qt[2])




