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

# this is not working
rows <- apply(mat, 1, function(row) {
  qt <- quantile( row, probs = c(0.1, 0.9) )
  q <- any( row < qt[1] | row > qt[2] )
  return(q)
} )
matF <- mat[ rows, ]

# --- this works fabulously ---
rows <- apply(mat, 1, function(row) var(row) < 1 )
matF <- mat[ rows, ]
dim(matF)


hist(mat[1,], breaks=40, col="blue")

d <- density(mat[1,])
plot(d)
polygon(d, col="red")

var(mat[1,]) < 1

qt <- quantile( mat[1,], probs = c(0.1, 0.9) )
qt
any(mat[1, ] > qt[2])




