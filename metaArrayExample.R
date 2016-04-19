###################################################
### code chunk number 1: metaArray.Rnw:44-70
###################################################
library(metaArray)
library(Biobase)
library(MergeMaid)
data(mdata)
common <- intersect(rownames(chen.poe),rownames(garber.poe)) 
common <- intersect(common, rownames(lapointe.poe))
chen.poe <- chen.poe[match(common, rownames(chen.poe)),]
garber.poe <- garber.poe[match(common, rownames(garber.poe)),]
lapointe.poe <- lapointe.poe[match(common, rownames(lapointe.poe)),]
vars <- list("var1","var2")
names(vars) <- names(chen.spl)
pdata1 <- new("AnnotatedDataFrame")
pData(pdata1) <- chen.spl
varLabels(pdata1) <- vars
sample1 <- new("ExpressionSet", exprs=chen.poe, phenoData = pdata1)
names(vars) <- names(garber.spl)
pdata2 <- new("AnnotatedDataFrame")
pData(pdata2) <- garber.spl
varLabels(pdata2) <- vars
sample2 <- new("ExpressionSet", exprs=garber.poe, phenoData = pdata2)
names(vars) <- names(lapointe.spl)
pdata3 <- new("AnnotatedDataFrame")
pData(pdata3) <- lapointe.spl
varLabels(pdata3) <- vars
sample3 <- new("ExpressionSet", exprs=lapointe.poe, phenoData = pdata3)
merged <- mergeExprs(sample1,sample2,sample3)


###################################################
### code chunk number 2: metaArray.Rnw:82-88 (eval = FALSE)
###################################################
## merged.cor <- intcor(merged)$avg.cor
## mData <- exprs(intersection(merged)); mCl <- NULL
## mData <- mData[merged.cor > median(merged.cor),]
## for(i in 1:length(merged)) {
##   mCl <- c(mCl, pData(merged[i])$metastasis)
## }


###################################################
### code chunk number 3: metaArray.Rnw:103-104
###################################################
em.draw(as.numeric(chen[1,]), cl=ncol(chen))


###################################################
### code chunk number 4: metaArray.Rnw:111-112 (eval = FALSE)
###################################################
## z.stat <- Zscore(merged)