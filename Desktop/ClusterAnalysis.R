#download dataset
#install.packages("downloader")
#library(downloader) 
#url <- "https://github.com/genomicsclass/tissuesGeneExpression/tree/master/data/tissuesGeneExpression.rda"
#filename <- "tissuesGeneExpression.rda"
#download(url, destfile=filename)

#download dataset
install.packages("devtools")
library(devtools)
install_github("genomicsclass/tissuesGeneExpression")

library(tissuesGeneExpression) 
data(tissuesGeneExpression)
dim(e) ##e contains the expression data

table(tissue) ##tissue[i] tells us what tissue is represented by e[,i]

#To illustrate one application of clustering, let’s pretend that we don’t know
#these are different tissues and are interested in clustering. 
#The first step is to compute the distance between each sample:
d <- dist( t(e) )

#Hierarchical clustering 
install.packages("rafalib")
library(rafalib)
mypar()
hc <- hclust(d)
hc
plot(hc,labels=tissue,cex=0.5)#plots dendogram 

#add color to visualize different tissues
myplclust(hc, labels=tissue, lab.col=as.fumeric(tissue), cex=0.5)

#add cuttoff value
myplclust(hc, labels=tissue, lab.col=as.numeric(tissue),cex=0.5)

#how do our clusters align with actual tissues
hclusters <- cutree(hc, h=120)
table(true=tissue, cluster=hclusters)

hclusters <- cutree(hc, k=8)
table(true=tissue, cluster=hclusters)


##K means

#As example, select first two genes to run 
set.seed(1)
km <- kmeans(t(e[1:2,]), centers=7)
names(km)

mypar(1,2)
plot(e[1,], e[2,], col=as.fumeric(tissue), pch=16) #color represents actual tissues
plot(e[1,], e[2,], col=km$cluster, pch=16) #color reps clusters defined by kmeans


table(true=tissue,cluster=km$cluster) #compare actual vs. expected results


#K means clustering using all genes
km <- kmeans(t(e), centers=7)
mds <- cmdscale(d)
mypar(1,2)
plot(mds[,1], mds[,2]) #MDSplot
plot(mds[,1], mds[,2], col=km$cluster, pch=16)


table(true=tissue,cluster=km$cluster) #how does this change results
abline(h=120)

#Heatmaps

#define color pallete
install.packages("RColorBrewer")
library(RColorBrewer)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

#install bioconductor

#install_bioc("genefilter")
library(genefilter)
rv <- rowVars(e)
idx <- order(-rv)[1:40]

#install.packages("gplots")
library(gplots) ##Available from CRAN
cols <- palette(brewer.pal(9, "GnBu"))[as.fumeric(tissue)] 
head(cbind(colnames(e),cols))

heatmap.2(e[idx,], labCol=tissue,
          trace="none",
          ColSideColors=cols,
          col=hmcol)



