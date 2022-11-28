#(1) Cluster patients using the hierarchical clustering algorithm (complete linkage), and cluster genes using the k-means algorithm (try K= 20, 30, 40, use 100 random starts for each K).
names = as.character(bulk[,1])
rownames(bulk) = names
bulk = as.matrix(bulk[-c(1)])

plot(apply(bulk,2,mean), main="mean of genes")
plot(apply(bulk,2,median), main="median of genes")
plot(apply(bulk,2,sd), main="SD of genes")
plot(apply(bulk,2,mad), main="MAD of genes")
heatmap(bulk) #using package

# Prepare Data
mydata <- na.omit(bulk) # listwise deletion of missing

# Clustering the Observations
data.dist=dist(t(mydata), , method = "euclidean") # distance matrix

# Hierarchical Clustering (Complete Linkage) for patients
hc.c <- hclust(data.dist, method="complete")

plot(hc.c, main="Complete Linkage", hang = -1, xlab="", sub="",ylab="")

set.seed(0)

# fit kmeans models with different number of clusters for genes

#k=20 (try k=20,30,40)
set.seed(0)
nclust <- 20
kms = kmeans(mydata, nclust, nstart=100) # transpose

km.clusters <- kms$cluster
table(km.clusters)

# make plots of 2-d clustered data
library(gplots) # for the color function "redgreen"
tdata = t(mydata)
image(1:nrow(tdata), 1:ncol(tdata), tdata[hc.c$order,order(kms$cluster)], col=redgreen(1024), main="k: 20")

abline(col=0, h=cumsum(kms$size)-0.5)

#######################################################
#k=30 (try k=20,30,40) k=the number of clusters for genes
set.seed(0)
nclust <- 30
kms = kmeans(mydata, nclust, nstart=100) # transpose

km.clusters <- kms$cluster
table(km.clusters)

image(1:nrow(tdata), 1:ncol(tdata), tdata[hc.c$order,order(kms$cluster)], col=redgreen(1024), main="k: 30")

abline(col=0, h=cumsum(kms$size)-0.5)

#######################################################
#k=40 (try k=20,30,40) k=the number of clusters for genes
set.seed(0)
nclust <- 40
kms = kmeans(mydata, nclust, nstart=100) # transpose

km.clusters <- kms$cluster
table(km.clusters)

# make plots of 2-d clustered data
image(1:nrow(tdata), 1:ncol(tdata), tdata[hc.c$order,order(kms$cluster)], col=redgreen(1024), main="k: 40")

abline(col=0, h=cumsum(kms$size)-0.5)

# some other packages
heatmap.2(mydata, col=redgreen(1024), scale="row", 
          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)

######################################################################
#(5) Repeat for log-transferred gene expression levels.

#(1) Cluster patients using the hierarchical clustering algorithm (complete linkage), and cluster genes using the k-means algorithm (try K= 20, 30, 40, use 100 random starts for each K).
# Log transfer data
mydata = (log(mydata+1)) 

# Clustering the Observations
data.dist=dist(t(mydata), , method = "euclidean") # distance matrix

# Hierarchical Clustering (Complete Linkage) for patients
set.seed(0)
hc.c <- hclust(data.dist, method="complete")

plot(hc.c, main="Complete Linkage", hang = -1, xlab="", sub="",ylab="")

set.seed(0)

# fit kmeans models with different number of cclusters for genes

#k=20 (try k=20,30,40) k=the number of clusters for genes
set.seed(0)
nclust <- 20
kms = kmeans(mydata, nclust, nstart=100) # transpose

km.clusters <- kms$cluster
table(km.clusters)

# make plots of 2-d clustered data
tdata = t(mydata)
image(1:nrow(tdata), 1:ncol(tdata), tdata[hc.c$order,order(kms$cluster)], col=redgreen(1024), main="k: 20")

abline(col=0, h=cumsum(kms$size)-0.5)

#######################################################
#k=30 (try k=20,30,40) k=the number of clusters for genes?
set.seed(0)
nclust <- 30
kms = kmeans(mydata, nclust, nstart=100) # transpose

km.clusters <- kms$cluster
table(km.clusters)

# make plots of 2-d clustered data
image(1:nrow(tdata), 1:ncol(tdata), tdata[hc.c$order,order(kms$cluster)], col=redgreen(1024), main="k: 30")

abline(col=0, h=cumsum(kms$size)-0.5)

#######################################################
# k=40 (try k=20,30,40) k=the number of clusters for genes
set.seed(0)
nclust <- 40
kms = kmeans(mydata, nclust, nstart=100) # transpose

km.clusters <- kms$cluster
table(km.clusters)

# make plots of 2-d clustered data
image(1:nrow(tdata), 1:ncol(tdata), tdata[hc.c$order,order(kms$cluster)], col=redgreen(1024), main="k: 40")

abline(col=0, h=cumsum(kms$size)-0.5)

# some other packages
heatmap.2(mydata, col=redgreen(1024), scale="row", 
          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
