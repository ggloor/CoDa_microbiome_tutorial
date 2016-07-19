
# this is the dataset I will present in Cambridge
# for the tutorial, we will use a random subset of 50 samples from each

The dataset contains 229 genera and 1457 samples. Samples are named thusly:
# ids have leading 700 stripped off
# samples named as follows:
# td_ Tongue Dorsum
# bm_ Buccal mucosa
# ak_ Attached Keratinzed Gingiva
# hp_ Hard Palate
# pt_ Palatine Tonsils
# sa_ Saliva
# up_ Subgingival Plaque (under plaque)
# op_ Supragingival Plaque (over plaque)

# read in the example table
# this is a subset of the HMP oral microbiome with 24 TD and 18 BM samples
# no real reason to choose these except they separate

# HMP oral V3-5 microbiome dataset
mouth <- read.table("~/git/compositions/oral/data/mouth_otu.txt", header=T, row.names=1, sep="\t")
taxon <- read.table("~/git/compositions/oral/data/taxon_names.txt", header=T, row.names=1, sep="\t")

# fix unknown taxa
# all p__, c__, o__ named
notF <- grep("f__", taxon[,1], invert=T)
fillF <- "f__Unknown;g__Unknown"

for(i in 1:length(notF)){
	taxon[notF[i],1] <- ( paste(taxon[notF[i],1], fillF, sep=";") )
}
notG <- grep("g__", taxon[,1], invert=T)
fillG <- "g__Unknown"

for(i in 1:length(notG)){
	taxon[notG[i],1] <- ( paste(taxon[notG[i],1], fillG, sep=";") )
}

z <- strsplit(as.character(taxon[,1]), ';')
tax.df <- as.data.frame(do.call(rbind,z))

rownames(tax.df) <- rownames(mouth)
colnames(tax.df) <- c("root", "phylum", "class", "order", "family", "genus")
tax.df$phylum <- gsub("p__", "", tax.df$phylum)
tax.df$genus <- gsub("g__", "", tax.df$genus)


d <- data.frame(mouth[,grep("ak_", colnames(mouth))], mouth[,grep("op_", colnames(mouth))])

# keep samples with > 1000 reads
d.col <- d[, which(apply(d,2,sum) > 1000)]

# many ways to filter
# keep OTUs with mean > 0.1 reads across all samples
d.subset <- d[ which(apply(d.col,1,mean) > .1),]
taxon.1 <- tax.df[ which(apply(d.col,1,mean) > .1),]

# this is the base dataset that everyone can use
write.table(d.subset, file="~/git/Host_Microbe_2016/data/ak_vs_op.txt", sep="\t", col.names=NA, quote=F)
write.table(taxon.1, file="~/git/Host_Microbe_2016/data/taxon.txt", sep="\t", col.names=NA, quote=F)

d <- data.frame(mouth[,grep("up_", colnames(mouth))], mouth[,grep("op_", colnames(mouth))])

# keep samples with > 1000 reads
d.col <- d[, which(apply(d,2,sum) > 1000)]

# many ways to filter
# keep OTUs with mean > 0.1 reads across all samples
d.subset <- d[ which(apply(d.col,1,mean) > .1),]
taxon.1 <- tax.df[ which(apply(d.col,1,mean) > .1),]

# this is the base dataset that everyone can use
write.table(d.subset, file="~/git/Host_Microbe_2016/data/up_vs_op.txt", sep="\t", col.names=NA, quote=F)
write.table(taxon.1, file="~/git/Host_Microbe_2016/data/up_vs_op_taxon.txt", sep="\t", col.names=NA, quote=F)

#############

################ OTU LEVEL
# some initial characterization of the whole dataset


d.subset <- read.table("~/git/Host_Microbe_2016/data/ak_vs_op.txt", sep="\t", header=T, row.names=1)
taxon.1 <- read.table("~/git/Host_Microbe_2016/data/ak_vs_op_taxon.txt", sep="\t", header=T, row.names=1)

library(zCompositions)
d.n0 <- cmultRepl(t(d.subset), method="CZM", label=0)
d.clr <- t(apply(d.n0, 1, function(x) log(x) - mean(log(x))))
pcx <- prcomp(d.clr)
rownames(pcx$rotation) <- taxon.1[rownames(pcx$rotation), "genus"]

PC1 <- paste("PC1: ", round(pcx$sdev[1]^2/sum(pcx$sdev^2),3), sep="")
PC2 <- paste("PC2: ", round(pcx$sdev[2]^2/sum(pcx$sdev^2),3), sep="")


biplot(pcx, cex=c(0.5,0.4), col=c("black", rgb(1,0,0,0.2)), var.axes=F, scale=0, xlab=PC1, ylab=PC2)

d.rand <- cbind(d.subset[,sample(grep("ak", colnames(d.subset)), 25)], d.subset[,sample(grep("op", colnames(d.subset)), 25)])
d.rand.filt <- d.rand[rowSums(d.rand) > 0,]


n.rand = 15
conds <- c(rep("ak", n.rand), rep("op" , n.rand))


x <- aldex.clr(d.rand)

x.e <- aldex.effect(x, conds)
x.t <- aldex.ttest(x, conds)
x.all <- data.frame(x.e,x.t)

tax <- "genus"
all.plot <- data.frame(x.all, tax.df[rownames(x.all), tax])
colnames(all.plot)[ncol(all.plot)] <- "taxon.level"
all.plot <- droplevels(all.plot)

no.sig <- all.plot$wi.eBH > 0.01
sig.pos <- all.plot$wi.eBH < 0.01 & all.plot$effect > 0
sig.neg <- all.plot$wi.eBH < 0.01 & all.plot$effect < 0

groups <- unique(all.plot$taxon.level)
ylim<-c(length(groups) - (length(groups)+0.5), length(groups) + 0.5)

xlim = c(min(-1 * max(abs(all.plot$effect))), max(all.plot$effect))
pdf("genus_stripchart.pdf")
par(mar=c(5,15,5,1), las=1, cex=0.7)
stripchart(effect ~ taxon.level, data=all.plot[no.sig,], col=rgb(0,0,0,0.3),method="jitter", jitter=0.2, pch=19, xlim=xlim, xlab="effect", main=tax)
stripchart(effect ~ taxon.level, data=all.plot[sig.pos,], col=rgb(0,0,1,0.3),method="jitter", jitter=0.2, pch=19, xlim=xlim, xlab="effect", add=T)
stripchart(effect ~ taxon.level, data=all.plot[sig.neg,], col=rgb(1,0,0,0.3),method="jitter", jitter=0.2, pch=19, xlim=xlim, xlab="effect", add=T)
abline(v=0, lty=2, lwd=2, col=rgb(0,0,0,0.3))
#draw horizonal lines
for (j in 0.5:(length(groups)+0.5)){
	abline(h=j, lty=3, col=rgb(0,0,0,0.3))
}
dev.off()


### MAKE A FULL COMPARISON FOR FOR THE UP OP SAMPLE PAIR, and the ak op sample pair

e.subset <- read.table("data/up_vs_op.txt", row.names=1, header=T)
e.x <- aldex.clr(e.subset)
e.conds <- c(rep("up", length(grep("up", colnames(e.subset))) ), rep("op", length(grep("op", colnames(e.subset)))) )
e.eff <- aldex.effect(e.x, e.conds)
e.tt <- aldex.ttest(e.x, e.conds)
e.all <- data.frame(e.eff,e.tt)


# SAVE DATA IN up_vs_op_aldex.txt
write.table(x.all, , file="~/git/Host_Microbe_2016/data/ak_vs_op_aldex.txt", sep="\t", col.names=NA, quote=F)
write.table(e.all, , file="~/git/Host_Microbe_2016/data/up_vs_op_aldex.txt", sep="\t", col.names=NA, quote=F)

############# GENUS LEVEL
# aggregate the data by genus count
mouth.genus <- aggregate(mouth, by=list(tax.df$genus),FUN=sum)
rownames(mouth.genus) <- mouth.genus$Group.1
mouth.genus$Group.1 <- NULL
mouth.genus.subset <- mouth.genus[rowSums(mouth.genus) > 0,]
write.table(mouth.genus.subset, file="~/git/Host_Microbe_2016/data/mouth_genus.txt", sep="\t", col.names=NA, quote=F)


gen.n0 <- cmultRepl(t(mouth.genus.subset), label=0, method="CZM")

gen.clr <- t(apply(gen.n0, 1, function(x) log(x) - mean(log(x)) ))
pcx.gen <- prcomp(gen.clr)

PC1 <- paste("PC1: ", round(pcx.gen$sdev[1]^2/sum(pcx.gen$sdev^2),3), sep="")
PC2 <- paste("PC2: ", round(pcx.gen$sdev[2]^2/sum(pcx.gen$sdev^2),3), sep="")
biplot(pcx.gen, cex=c(0.6,0.5), var.axes=F, scale=0, xlab=PC1, ylab=PC2)
abline(v=0)
abline(h=0)

d.g <- data.frame(mouth.genus.subset[,grep("ak_", colnames(mouth.genus.subset))], mouth.genus.subset[,grep("op_", colnames(mouth.genus.subset))])
d.g.subset <- d.g[apply(d.g, 1, mean) > 0.1, ]
d.g.n0 <- cmultRepl(t(d.g.subset), label=0, method="CZM")
d.g.clr <- t(apply(d.g.n0, 1, function(x) log(x) - mean(log(x)) ))
pcx.d.g <- prcomp(d.g.clr)

PC1 <- paste("PC1: ", round(pcx.d.g$sdev[1]^2/sum(pcx.d.g$sdev^2),3), sep="")
PC2 <- paste("PC2: ", round(pcx.d.g$sdev[2]^2/sum(pcx.d.g$sdev^2),3), sep="")
biplot(pcx.d.g, cex=c(0.9,0.7), var.axes=F, scale=0, xlab=PC1, ylab=PC2)
abline(v=0)
abline(h=0)

conds.g <- c(rep("A", length(grep("ak", colnames(d.g.subset)))  ), rep("O", length(grep("op", colnames(d.g.subset)))) )
x.g <- aldex.clr(d.g.subset)
x.g.e <- aldex.effect(x.g, conds.g)
x.t.t <- aldex.ttest(x.g, conds.g)
x.g.all <- data.frame(x.g.e, x.t.t)


