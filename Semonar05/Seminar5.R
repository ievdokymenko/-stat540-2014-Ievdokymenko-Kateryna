library(lattice) # if you don't already have this loaded ...
prDat <- read.table("GSE4051_data.tsv")
str(prDat, max.level = 0)
prDes <- readRDS("GSE4051_design.rds")
str(prDes)

###############
# takes as input the Affymetrix probeset ID(s) and gives as output a data.frame
prepareData <- function(keepGenes) {
  miniDat <- subset(prDat, rownames(prDat) %in% keepGenes)
  # reshape the data to be tall and skinny
  miniDat <- data.frame(gExp = as.vector(t(as.matrix(miniDat))),
                        gene = factor(rep(rownames(miniDat), each = ncol(miniDat)),
                                      levels = keepGenes))
  miniDat <- suppressWarnings(data.frame(prDes, miniDat))
  return(miniDat)
}

# test the function
jDat <- prepareData(c("1419655_at","1438815_at"))
str(jDat)
head(jDat)
tail(jDat)
stripplot(gExp ~ devStage | gene, jDat,
          group = gType, jitter.data = TRUE,
          auto.key = TRUE, type = c('p', 'a'), grid = TRUE)

###############
# a function to stripplot a mini-dataset
makeStripplot <- function(df) {
  stripplot(gExp ~ devStage | gene, df,
            group = gType, jitter.data = TRUE,
            auto.key = TRUE, type = c('p', 'a'), grid = TRUE)
}

# test the stripplot function
makeStripplot(jDat)
makeStripplot(newDat <- prepareData("1456341_a_at"))
str(newDat)
head(newDat)

###############
# test for a difference in expected gene expression for probeset "1456341_a_at" 
# at developmental stage P2 vs. 4 weeks post-natal 
# (ignoring genotype, i.e. lump the wild types and knockouts together).
# Let's assume a common variance in the two groups.
ttInput <- subset(prepareData("1456341_a_at"), devStage %in% c("P2","4_weeks"))
t.test(gExp ~ devStage, ttInput)

###############
# Fit a linear model with a categorical covariate
# In other words, do "one-way ANOVA". Focus on probeset "1438786_a_at".
# Here's what the data looks like:
makeStripplot(mDat <- prepareData("1438786_a_at"))

# Let's focus just on the wild type data for now. 
# Model expression as a function of the devStage factor.
mFit <- lm(formula = gExp ~ devStage, data = mDat, subset = gType == "wt")
summary(mFit)

# Vet your inferential results: does the intercept look plausible given the plot? 
# How about the devStageP2 effect, etc.?
# Yes, looks right.

###############
# Test whether the P2 and P10 effects are equal or, equivalently, whether their difference is equal to zero.
# Construct the contrast matrix to form the difference between the P2 and P10 effects.
contMat = matrix(c(0, 1, 0, -1, 0), nrow=1)
(obsDiff <- contMat %*% coef(mFit))

# check that this really is the observed difference in sample mean for the wild type mice, P2 vs. P10.
(sampMeans <- aggregate(gExp ~ devStage, mDat, FUN = mean,
                        subset = gType == "wt"))
with(sampMeans, gExp[devStage == "P2"] - gExp[devStage == "P10"])

# Now we need the (estimated) standard error for our contrast.
# The variance-covariance matrix of the parameters estimated in the original model can be obtained with vcov():
vcov(mFit)
# Let's check that this is really true. If we take the diagonal elements and take their square root, 
# they should be exactly equal to the standard errors reported for out original model.
summary(mFit)$coefficients[ , "Std. Error"]
sqrt(diag(vcov(mFit)))

# a typical matrix of inferential results can be obtained like so:
summary(mFit)$coefficients

# Returning to our test of the P2 vs. P10 contrast, recall that the variance-covariance matrix of a contrast obtained as
(estSe <- contMat %*% vcov(mFit) %*% t(contMat))
# Now we form a test statistic as an observed effect divided by its estimated standard error:
(testStat <- obsDiff/estSe)

# Under the null hypothesis that the contrast equals zero, i.e. there is no true difference in mean 
# for expression at P2 and P10 in wild type mice for this gene, the test statistic has a t distribution 
# with n−p=20−5=15 degrees of freedom. We compute a two-sided p-value and we're done.
2 * pt(abs(testStat), df = df.residual(mFit), lower.tail = FALSE)
# Not surprisingly, this p-value is rather large and we conclude there is no difference.

##################
# Fit a linear model with two categorical covariates
# Let's focus on probeset "1448690_at"
makeStripplot(oDat <- prepareData("1448690_at"))
str(oDat)
# Fit a linear model with covariates gType and devStage and include their interactions.
oFitBig <- lm(formula = gExp ~ gType * devStage, data=oDat)
summary(oFitBig)$coef
# Vet the results. Is the intercept plausible? How about the various effects? 
# Do the ones with small p-values, e.g. meeting a conventional cut-off of 0.05, look 'real' to you?
# Intercept and effects seem reasonable

# Fit a related, smaller model with the same covariates, but this time omit the interaction.
oFitSmall <- lm(formula = gExp ~ gType + devStage, data=oDat)

# Now let's determine if the interaction terms are truly necessary. 
# From the plot, the case for interaction seems very weak. This can be assessed with 
# an F test that essentially looks at the reduction in the sum of squared residuals due to
# using a larger, more complicated model and determines if it is "big enough" given 
# the number of additional parameters used. Recall the anova() function can take two fitted models, 
# one nested within the other, and perform this test.
anova(oFitSmall, oFitBig)
# With a p-value awfully close to one, we confirm that, no, there is no evidence for interaction in this particular case.

#Looking at a more interesting gene
makeStripplot(oDat2 <- prepareData("1429225_at"))
oFitBig2 <- lm(formula = gExp ~ gType * devStage, data=oDat2)
oFitSmall2 <- lm(formula = gExp ~ gType + devStage, data=oDat2)
anova(oFitSmall2, oFitBig2)
# the interaction here is highly statistically significant.

################
# further work
# use data aggregation strategies from last week to do some of the same work for small sets of genes
genes = c("1419655_at","1438815_at")
jDat <- prepareData(genes)
library(plyr)
(aggr <- ddply(jDat, ~ sidChar + devStage + gType, summarize, gExp = mean(gExp), gene = paste(genes, collapse=' & ')))
makeStripplot(aggr)

# fit linear and quadratic models to the expression data for one or several genes
