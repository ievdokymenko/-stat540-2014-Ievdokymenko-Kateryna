prDat <- read.table("GSE4051_data.tsv")
str(prDat, max.level = 0)
prDes <- readRDS("GSE4051_design.rds")
str(prDes)
#  extract the data for one gene
set.seed(987)
(theGene <- sample(1:nrow(prDat), 1))
pDat <- data.frame(prDes, gExp = unlist(prDat[theGene, ]))
str(pDat)
aggregate(gExp ~ gType, pDat, FUN = mean)
stripplot(gType ~ gExp, pDat)
ttRes <- t.test(gExp ~ gType, pDat)
str(ttRes)
ttRes$statistic
ttRes$p.value


(theGene2 <- sample(1:nrow(prDat), 1))
pDat2 <- data.frame(prDes, gExp = unlist(prDat[theGene2, ]))
aggregate(gExp ~ gType, pDat2, FUN = mean)
stripplot(gType ~ gExp, pDat2)
ttRes2 <- t.test(gExp ~ gType, pDat2, var.equal=FALSE)
ttRes2cv <- t.test(gExp ~ gType, pDat2, var.equal=TRUE)
wRes <- wilcox.test(gExp ~ gType, pDat2)
ksRes <- ks.test(pDat2$gType, pDat2$gExp)

ksRes <- suppressWarnings(ks.test(pDat2$gType, pDat2$gExp))

#Create a dataframe for tests results
data.frame("method"=c(ttRes2$method, ttRes2cv$method, wRes$method, ksRes$method), 
           p.value=c(ttRes2$p.value, ttRes2cv$p.value, wRes$p.value, ksRes$p.value))
# Load the data

kDat <- read.table("GSE4051_MINI.txt")
str(kDat)
kMat <- as.matrix(kDat[c('crabHammer', 'eggBomb', 'poisonFang')])
str(kMat)

#ways to get the median expression (in apply function 2 means columns) for a column
median(kMat[ , 1])
median(kMat[ , "crabHammer"])
apply(kMat, 2, median)
apply(kMat, 2, quantile, probs = 0.5)

apply(kMat, 2, quantile, probs = c(0.25, 0.75))

# minimum gene expression for each sample
apply(kMat, 1, min)
colnames(kMat)[apply(kMat, 1, which.min)]


#Computing row-wise sums
rowSums(kMat)

# check do the rowSums and apply functions give the same results
all.equal(rowSums(kMat), apply(kMat, 1, sum))
at)

#average expression of eggBomb for different levels of devStage

aggregate(eggBomb ~ devStage, kDat, FUN = mean)
aggregate(eggBomb ~ gType * devStage, kDat, FUN = mean)

# or using plyr
ddply(kDat, ~ devStage, summarize, exp = mean(eggBomb))


# choose 6 genes: 3 are interesting, 3 are boring
keepGenes <- c("1431708_a_at", "1424336_at", "1454696_at",
               "1416119_at", "1432141_x_at", "1429226_at" )

miniDat <- subset(prDat, rownames(prDat) %in% keepGenes)

#data reshaping 
as.matrix(miniDat)
t(as.matrix(miniDat))
as.vector(t(as.matrix(miniDat))
miniDat <- data.frame(gExp = as.vector(t(as.matrix(miniDat))), gene = factor(rep(rownames(miniDat), each = ncol(miniDat)), levels = keepGenes))
          
miniDat2 <- data.frame(prDes, miniDat)
head(miniDat2)
miniDat2 <- suppressWarnings(data.frame(prDes, miniDat))
          
#plot  3 'hits' and 3  boring genes
stripplot(gType ~ gExp | gene, miniDat2, scales = list(x = list(relation = "free")),group = gType, auto.key = TRUE)
          
          
# two group comparisons for each of these 6 genes
          
someDat <- droplevels(subset(miniDat2, gene == keepGenes[1]))
          t.test(gExp ~ gType, someDat)
#take-home work

#pick 100 genes at random
genesSet <-  sample(1:nrow(prDat), 100)

dataset <- rownames(prDat)[genesSet]

q <- subset(prDat, rownames(prDat) %in% dataset)

z <- data.frame(gExp = as.vector(t(as.matrix(q))),
                gene = factor(rep(rownames(q), each = ncol(q)),
                              levels = dataset))
z
y <- suppressWarnings(data.frame(prDes, z))
ttRes <- ddply(y, ~ gene, function(z) {
  zz <- t.test(gExp ~ gType, z)
  zy <- suppressWarnings(wilcox.test(gExp ~ gType, z))
  round(c(pVal_ttest = zz$p.value, pVal_W = zy$p.value ), 4)
})
ttRes