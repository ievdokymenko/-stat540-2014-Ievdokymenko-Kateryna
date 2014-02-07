# loas libraries
library("ggplot2", lib.loc="C:/Program Files/R/R-2.15.2/library")
library("lattice", lib.loc="C:/Program Files/R/R-2.15.2/library")
library("latticeExtra", lib.loc="C:/Program Files/R/R-2.15.2/library")
install.packages("RColorBrewer")

# load photoRec dataset
prDat <- read.table("GSE4051_data.tsv")
str(prDat, max.level = 0)
prDes <- readRDS("GSE4051_design.rds")
str(prDes)

# choose 50 datasets
set.seed(2)
(yo <- sample(1:nrow(prDat), size = 50))
hDat <- prDat[yo, ]

# convert prDat into matrix

hDat <- as.matrix(t(hDat))
rownames(hDat) <- with(prDes,
                       paste(devStage, gType, sidChar, sep="_"))
str(hDat)

#build heatmaps
heatmap(hDat, Rowv = NA, Colv = NA, scale="none", margins = c(5, 8))
jGraysFun <- colorRampPalette(brewer.pal(n = 9, "Greys"))
heatmap(hDat, Rowv = NA, Colv = NA, scale="none", margins = c(5, 8),
        col = jGraysFun(256))
# scaling within column
heatmap(hDat, margins = c(6, 8), scale=c("column"))
# or scaling with rows
heatmap(hDat, margins = c(6, 8), scale=c("row"))


#overplotting 

#select two samples at random
set.seed(924)
(yo <- sample(1:ncol(prDat), size = 2))

bDat <- data.frame(y = prDat[[yo[1]]], z = prDat[[yo[2]]])
str(bDat)

# plot against each other (lattice library looks like ggplot)
xyplot(y ~ z, asp = 1, par.settings = ggplot2like())


xyplot(y ~ z, asp = 1, par.settings = ggplot2like())
xyplot(y ~ z, bDat, asp = 1, par.settings = ggplot2like())

smoothScatter(y ~ z, bDat, asp = 1)
y <- prDat[[yo[1]]]
z <- prDat[[yo[2]]]

smoothScatter(y ~ z, asp = 1)
xyplot(y ~ z, asp = 1, panel = panel.smoothScatter, nbin = 150)


# select several variables
(yo <- sample(1:ncol(prDat), size = 4))
pairDat <- subset(prDat, select = yo)
str(pairDat)
# plot several variables against each other
pairs(pairDat,panel = function(...) smoothScatter(..., add=TRUE))

pairs(pairDat,panel = function(...) smoothScatter(...))
install.packages("hexbin")
library("hexbin", lib.loc="C:/Program Files/R/R-2.15.2/library")

splom(pairDat)
splom(pairDat, panel = panel.smoothScatter, raster = TRUE)
hexplom(pairDat)

