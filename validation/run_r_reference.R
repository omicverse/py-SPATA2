args <- commandArgs(trailingOnly=TRUE)
coords <- as.matrix(read.csv(args[[1]]))
eps <- as.numeric(args[[2]])
flags <- dbscan::dbscan(coords, eps=eps, minPts=4)$cluster == 0
write.csv(data.frame(outlier=as.integer(flags)),args[[3]],row.names=FALSE)
