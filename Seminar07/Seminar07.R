dat <- read.table("bottomly_count_table.tsv", header = TRUE, row.names = 1)
des <- read.table("bottomly_phenodata.tsv", header = TRUE, row.names = 1)

#### edgeR
library(edgeR)

# create a 'group' object describing which group each sample belongs to:
group <- factor(c(rep("1", 10), rep("2", 11)))

# this produces an object of type DGEList with can be manipulated in a
# similar way to any other list object in R
dge.glm <- DGEList(counts = dat, group = group)
#str(dge.glm)
#names(dge.glm)
#dge.glm[["samples"]]
nrow(dge.glm[[1]])
ncol(dge.glm[[1]])

design <- model.matrix(~group)
dge.glm.com.disp <- estimateGLMCommonDisp(dge.glm, design, verbose = TRUE)
dge.glm.trend.disp <- estimateGLMTrendedDisp(dge.glm.com.disp)
dge.glm.tag.disp <- estimateGLMTagwiseDisp(dge.glm.trend.disp, design)

# plot the tagwise dispersion against log2-CPM (counts per million)
plotBCV(dge.glm.tag.disp)

fit <- glmFit(dge.glm.tag.disp, design)
#colnames(coef(fit))

lrt <- glmLRT(fit, coef = 2)
#topTags(lrt)

tt.glm <- topTags(lrt, n = Inf)
#class(tt.glm)
#nrow(tt.glm$table[tt.glm$table$FDR < 0.01, ])

interestingSamples <- rownames(tt.glm$table[tt.glm$table$FDR < 1e-50, ])
cpm(dge.glm.tag.disp)[interestingSamples, ]

summary(de.glm <- decideTestsDGE(lrt, p = 0.05, adjust = "BH"))

# plotting the tagwise log fold changes against log-cpm
tags.glm <- rownames(dge.glm.tag.disp)[as.logical(de.glm)]
plotSmear(lrt, de.tags = tags.glm)
abline(h = c(-2, 2), col = "blue")

#### Exercises
# filter the data and remove any gene that has count equal to zero across all samples 
dge.glm.clean1 <- DGEList(counts = dat, group = group, remove.zeros=T)
str(dge.glm.clean1)

# filter the data and remove any gene that has count equal to zero in at least one sample in each genotype group
dge.glm.clean2.toremove <- apply(dge.glm, 1, function(x) {
  return(any(x[1:10] == 0) & any(x[11:21] == 0))
})
dge.glm.clean2.indices <- which(dge.glm.clean2.toremove != T)
dge.glm.clean2 <- dge.glm[dge.glm.clean2.indices, ]
str(dge.glm.clean2)

#### DESeq
# the differential expression analysis of the same dataset using DESeq.
library(DESeq)
# reading in the same count table data and grouping information
deSeqDat <- newCountDataSet(dat, group)
head(counts(deSeqDat))

# estimate the size factors to account for differences in library coverage and estimate the variance:
deSeqDat <- estimateSizeFactors(deSeqDat)
sizeFactors(deSeqDat)

deSeqDat <- estimateDispersions(deSeqDat)
# plotting the estimated dispersions against the mean normalized counts
plotDispEsts(deSeqDat)

# fit the model and examine the results
results.DESeq <- nbinomTest(deSeqDat, levels(group)[1], levels(group)[2])
str(results)
plotMA(results)

#### Voom & limma
library(limma)
norm.factor <- calcNormFactors(dat)
dat.voomed <- voom(dat, design, plot = TRUE, lib.size = colSums(dat) * norm.factor)
dat.voomed

fit <- lmFit(dat.voomed, design)
fit <- eBayes(fit)
topTable(fit)

#### Exercise
# Choose a specific threshold for the adjusted p value, find the genes identified 
# as differentially expressed using each of edgeR, DESeq and voom+limma. 
# Compare the number of genes in these 3 lists, and draw a venn digram demonstrating
# the overlap (if any!).

results.edgeR <- tt.glm$table
results.voom <- topTable(fit, coef = 2, adj = "BH", n = Inf)

limit <- 0.01
top.edgeR <- results.edgeR[which(results.edgeR$FDR < limit), ]
top.DESeq <- results.DESeq[which(results.DESeq$padj < limit), ]
top.voom <- results.voom[which(results.voom$adj.P.Val < limit), ]

gnames <- rownames(dat)
# filter out indices of top genes
dat.indices.edgeR <- which(gnames %in% rownames(top.edgeR))
dat.indices.DESeq <- which(gnames %in% top.DESeq$id)
dat.indices.voom <- which(gnames %in% rownames(top.voom))
# prepare count lists
counts.edgeR <- rep(0, length(gnames))
counts.DESeq <- rep(0, length(gnames))
counts.voom <- rep(0, length(gnames))
# fill count lists
counts.edgeR[dat.indices.edgeR] <- 1
counts.DESeq[dat.indices.DESeq] <- 1
counts.voom[dat.indices.voom] <- 1

venn.df <- data.frame(edgeR = counts.edgeR, 
                      DESeq = counts.DESeq, 
                      voom = counts.voom, 
                      row.names = gnames)
venn.counts <- vennCounts(venn.df)
vennDiagram(venn.counts)
