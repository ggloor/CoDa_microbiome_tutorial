---
title: "explore schurch RNA-seq set"
author: "gg"
date: "March 19,2017"
output: pdf_document
---

```{r lib, include=FALSE, message=FALSE, warning=FALSE, echo=TRUE}
# CRAN
# install.packages("mypackage") propr, zCompositions, car

#bioconductor
# source("https://bioconductor.org/biocLite.R")
# biocLite("ALDEx2") and DESeq and edgeR

# github
# install.packages('devtools')
# devtools::install_github('ggloor/CoDaSeq/CoDaSeq')

library(zCompositions)
library(ALDEx2)
library(CoDaSeq)

plot=FALSE
```

```{r explore1, message=FALSE, warning=FALSE, echo=TRUE, fig.cap='Initial exploratory PCA plot of the SNF1 KO dataset. .', fig.height=4, fig.width=14}

d <- read.table("data/countfinal2.tsv", header=T, row.names=1, sep="\t", stringsAsFactors=F, comment.char="")

m <- read.table("data/ERP004763_sample_mapping.tsv", , header=T, row.names=1, sep="\t", stringsAsFactors=F, comment.char="")

library(CoDaSeq)
library(zCompositions)

min(apply(d,1,sum))
# at least one gene has 0 in every sample

# need to remove genes with 0 reads in all samples
d.n0 <- codaSeq.filter(d, min.reads=0, min.prop=0, min.count=0, samples.by.row=FALSE)

# replace 0 values
d.czm <- cmultRepl(t(d.n0), method="CZM", label=0)

# convert to clr
d.clr <- codaSeq.clr(d.czm)

# do the SVD
d.pcx <- prcomp(d.clr)

# make the biplot
biplot(d.pcx, var.axes=F, scale=0, cex=c(0.5,0.1))

```

```{r agg, message=FALSE, warning=FALSE, echo=TRUE, fig.cap='aggregate the dataset by techrep.', fig.height=4, fig.width=14}

# make a vector of names concatenated with biological replicate number
nms.agg <- paste(m[,2], m[,3], sep=".")

#sum by the set of unique names
d.agg <- aggregate(t(d), by=list(nms.agg), FUN=sum)
rownames(d.agg) <- d.agg$Group.1
d.agg$Group.1 <- NULL

# remove rows with 0 counts
# returns data with samples by columns!!!!
d.agg.gt0 <- codaSeq.filter(d.agg, min.reads=0, min.prop=0, min.count=0, samples.by.row=TRUE)

# estimate 0 values (zCompositions)
d.agg.n0 <- cmultRepl(t(d.agg.gt0), method="CZM", label=0)

# clr transform
d.agg.n0.clr <- codaSeq.clr(d.agg.n0)

# do the SVD and get the total variance
pcx.agg  <- prcomp(d.agg.n0.clr)
mvar.agg.clr <- sum(pcx.agg$sdev^2)
PC1.var <- pcx.agg$sdev[1]^2/mvar.agg.clr
PC2.var <- pcx.agg$sdev[2]^2/mvar.agg.clr

# scree plot
plot(pcx.agg$sdev[1:10]^2/mvar.agg.clr)

# check it
par(mfrow=c(1,1))
biplot(pcx.agg, var.axes=F, scale=0, cex=c(1,.5))

```

```{r outlier, message=FALSE, warning=FALSE, echo=TRUE, fig.cap='outliers', fig.height=4, fig.width=14}
# get the outliers from each group. See codaSeq.outlier
# get WT indices
WT <- grep("WT", rownames(d.agg))
# subset
WT.agg <- d.agg[WT,]

# filter
wt.gt0 <- codaSeq.filter(WT.agg, min.reads=0, min.prop=0, min.count=0, samples.by.row=TRUE)

# estimate 0 values (zCompositions)
wt.agg.n0 <- cmultRepl(t(wt.gt0), method="CZM", label=0)

# clr transform
wt.agg.n0.clr <- codaSeq.clr(wt.agg.n0)

# SVD
pcx.wt  <- prcomp(wt.agg.n0.clr)
mvar.wt.clr <- sum(pcx.wt$sdev^2)

# plot
par(mfrow=c(1,1))
biplot(pcx.wt, var.axes=F, scale=0, cex=c(2,.5))

# make a list of names to keep. found in $good
WT.g <- codaSeq.outlier(wt.agg.n0.clr, plot.me=TRUE)

# now do the same for SNF2
# get SNF indices
SNF <- grep("SNF", rownames(d.agg))
# subset
SNF.agg <- d.agg[SNF,]

SNF.gt0 <- codaSeq.filter(SNF.agg, min.reads=0, min.prop=0, min.count=0, samples.by.row=TRUE)

# estimate 0 values (zCompositions)
SNF.agg.n0 <- cmultRepl(t(SNF.gt0), method="CZM", label=0)

# clr transform
SNF.agg.n0.clr <- codaSeq.clr(SNF.agg.n0)

pcx.SNF  <- prcomp(SNF.agg.n0.clr)
mvar.SNF.clr <- sum(pcx.SNF$sdev^2)

par(mfrow=c(1,1))
biplot(pcx.SNF, var.axes=F, scale=0, cex=c(2,.5))

SNF.g <- codaSeq.outlier(SNF.agg.n0.clr, plot.me=TRUE)

```

```{r good_data_pca, message=FALSE, warning=FALSE, echo=TRUE, fig.cap='outliers', fig.height=4, fig.width=14}

# make a dataset of only the non-outlier samples
d.good <- rbind(d.agg[SNF.g$good,],d.agg[WT.g$good,])

# filter
d.good.gt0 <- codaSeq.filter(d.good, min.reads=0, min.prop=0, min.count=0, samples.by.row=TRUE)

# estimate 0 values (zCompositions)
d.good.agg.n0 <- cmultRepl(t(d.good.gt0), method="CZM", label=0)

# clr transform
d.good.agg.n0.clr <- codaSeq.clr(d.good.agg.n0)

# SVD
pcx.good  <- prcomp(d.good.agg.n0.clr)
mvar.good <- sum(pcx.good$sdev^2)

# plot and save
par(mfrow=c(1,1))
biplot(pcx.good, var.axes=F, scale=0, cex=c(2,.5))

write.table(d.good.gt0, file="data/filtered_table.txt", sep="\t", quote=F, col.names=NA)
```

```{r aldex}
d.good.gt0 <- read.table("data/filtered_table.txt",row.names=1, header=T, sep="\t", stringsAsFactors=F, comment.char="")

# check for differential abundance
library(ALDEx2)

# make a vector of conditions

# you guys to make the distribution of possible values
x <- aldex.clr(d.good.gt0)

conds <- c(rep("SNF", length(SNF.g$good)), rep("WT", length(WT.g$good)))
# me
x <- aldex.clr(d.good.gt0, conds, mc.samples=16)

x.e <- aldex.effect(x, conds)
x.t <- aldex.ttest(x, conds)
x.all <- data.frame(x.e,x.t)


# now explore BLand Altman Plot
aldex.plot(x.all, type="MA")

# volcano plot
plot(x.all$diff.btw, x.all$wi.eBH, log="y"
plot(x.all$wi.eBH,x.all$diff.btw, log="x", col="grey")
points(x.all$wi.eBH[x.all$effect > 4],x.all$diff.btw[x.all$effect > 4], col="red", pch=19)
points(x.all$wi.eBH[x.all$effect < -4],x.all$diff.btw[x.all$effect < -4], col="blue",pch=19)

# effect size plot
aldex.plot(x.all, type="MW")

# we typically use effect sizes of 2 or more

# subset to find those with very large effect sizes
rownames(x.all)[abs(x.all$effect) > 4]
```
```{r edgeR}

# https://bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

library(edgeR)

# from ALDEx above
group <- factor(conds)

y <- DGEList(counts=d.good.gt0,group=group)

y <- calcNormFactors(y)

y <- estimateDisp(y)

# plot a non-metric multidimensional scaling plot of the data
# by default uses the 500 genes with the largest difference between samples
# NMDS is essentially PCA on the rank order differences rather than actual variances
# note that in this dataset, it looks very similar to the plot from the clr values

plotMDS(y)

# how well does the data fit the model?
# looking for the blue and red lines to be coincident
plotBCV(y)

# this is the differential abundance test
# exact test using negative binomial
et <- exactTest(y)
# hack to get the data for plotting later
tt <- topTags(et, sort.by="none", n=6236)

summary(et <- decideTestsDGE(et))
detags <- rownames(y)[as.logical(et)]


# you can plot the MA plot
# compare to the MA plot from ALDEx, they are rather similar for this dataset
plotSmear(y, de.tags=detags)
abline(h=c(-1, 1), col="blue", lty=2)

```

```{ r DESeq}
# https://www.bioconductor.org/packages/3.3/bioc/vignettes/DESeq/inst/doc/DESeq.pdf

library(DESeq)

# we will use the same input table
# set up the metadata
design <- data.frame(
  row.names = colnames(d.good.gt0),
  condition = conds,
  libType="single-end"
)

# for convenince so I can cut and paste
countTable <- d.good.gt0
condition <- design$condition
condition <- factor(condition)

# make the base data table
cds = newCountDataSet( countTable, condition )

# normalize read counts
cds = estimateSizeFactors( cds )
# observe
sizeFactors( cds )

# estimate the variance and fit to a negative binomial
cds = estimateDispersions( cds )

# observe how well it fits the data
# how does this differ and what is similar to edgeR?
plotDispEsts( cds )

# find DE genes
res = nbinomTest( cds, "SNF", "WT" )

# plot the MA plot, compare to ALDEx2 and edgeR
plotMA(res)

# it is useful to look at the p value distribution
par(mfrow=c(1,3))
hist(res$pval, breaks=100, main="DESeq") #DESeq
hist(tt[[1]]$PValue, breaks=100, col=rgb(1,0,0,0.1), main="edgeR") # edgeR
hist(x.all$we.ep, breaks=100, col=rgb(0,0,1,0.1), main="ALDEx2") # ALDEx2

par(mfrow=c(1,2))
plot(res$pval, tt[[1]]$PValue, log="xy")
plot(res$pval, x.all$we.ep, log="xy")


# get lists of differential genes by tool
# use cutoff of BH of 0.05

sig.ald <- rownames(x.all)[x.all$we.eBH < 0.05]
sig.des <- res$id[res$padj < 0.05]
sig.edg <- detags

write.table(sig.ald, file="sig.ald.txt", quote=FALSE, row.names=F, col.names=F)
write.table(sig.des, file="sig.des.txt", quote=FALSE, row.names=F, col.names=F)
write.table(sig.edg, file="sig.edg.txt", quote=FALSE, row.names=F, col.names=F)


```


```{r, annotation}
# you can supply lists of genes of interest, or call directly from within edgeR for this
# see the edgeR documentation for an example of GO terms, only metazoan models
# we can do KEGG term annotation though

de.aldex <- rownames(x.all)[x.all$effect > 2]
write.table(de.aldex, file="aldex_sig.effect.txt", sep="\t", col.names=F,row.names=F, quote=F)

# http://www.kegg.jp/kegg/tool/map_pathway1.html

# kegga is going to an external database (KEGG)
kegg.aldex <- kegga(de.aldex, species.KEGG="sce")

```
