Intro <- 
"Generic workflow to identify the yeast transcriptome in a GEO data set
most similar to a query transcriptome, e.g., PsAvh172-overexpressing yeast cells.
This workflow assumes GEO data set is a yeast microarray study with biological replicates 
and can account for duplicated probes (spots) on a microarray."
cat(Intro)

# !! REPLACE THESE PATHS AND FILE NAMES AS NEEDED!!
WORKdir <- ("~/P001DrMorgan") # REPLACE with your work directory path
GSEname <- "GSE54539" # REPLACE with the GSE number of the data set to download from GEO (no spaces in name!)
GLdir <- "LabNotebooks/GeneListsYeast" # REPLACE with your gene list directory path
GLfile <- "Saccharomyces cerevisiae-KEGG_Pathway_Windows.gmt" # REPLACE with your gene list file name
# Load previously calculated romerResults from query profile; e.g., see RNAseq_workflow_PsAvh172.R
r.Query <- read.delim("RScripts/r2.PsAvh172",header = TRUE, stringsAsFactors = FALSE) # romer output using default nrot

################################################################
library(Biobase)
library(GEOquery)
library("limma", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")
library("GSEABase", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")
library("pheatmap", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")
################################################################

# Create and move to a new folder for output files
WD <- file.path(WORKdir, GSEname) # New folder name
dir.create(WD) # Create folder
setwd(WD)

# Download desired data set from GEO
gse <- getGEO(GSEname, GSEMatrix =TRUE) # This produces a "Large list" with the desired ExpressionSet as the last element
# load series and platform data from GEO
ESet <- gse[[1]]
ex <- exprs(gse[[1]]) # Extract matrix of M values
cat(paste(GSEname, "data matrix initially has", dim(ex)[1], "probes and", dim(ex)[2], "samples. "))
# convert missing data values (NA) to 0
if (anyNA(ex)) {
  ex[which(is.na(ex))] <- 0
} 
# log2 transform expression values, if not already
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T)) # examines distribution of values
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2) 
# if LogC TRUE (not yet transformed), then take log2 of each positive value and replace original matrix in the ExpressionSet
if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gse) <- log2(ex) 
} 
rm(LogC); rm(qx)  # Clean up

# uncomment following to see more information about each sample 
# head(phenoData(ESet)$geo_accession) # GSM number of each sample
head(levels(phenoData(ESet)$title)) # Long name of each sample

# extract information about features from ESet
features <- featureData(ESet) # Put featureData in its own annotated dataframe
# head(varMetadata(features)) # Gives structure of featureData
# head(features$ID) # ID of the features
# head(levels(features$ORF)) # SGD gene names

## Construct data frame with information needed for downstream analysis
ex.df <- data.frame(PROBE_ID = features$PROBE_ID, ORF = features$ORF, ex)
ex.ORF.df <- ex.df[which(ex.df$ORF != ""),] # Only keep rows with a named ORF; removes spot controls, etc.

## If each probe is duplicated on microarry, then rearrange data frame to treat duplicate spots as technical replicates
if(anyDuplicated(ex.ORF.df$PROBE_ID) > 0) {
  if(length(ex.ORF.df)/length(levels(ex.ORF.df)) < 2) {
    warning("Is the complete probe set duplicated?") # Warn if only a subset of probes are duplicated
  }
  dup.no = 2
  # Sort by probe id so duplicate probes are in consecutive rows
  ex.ordered.df <- ex.ORF.df[order(ex.ORF.df$PROBE_ID),]
  # Move values from duplicate spots to duplicate columns (to BLOCK in romer)
  ex.odd.df <- ex.ordered.df[!duplicated(ex.ordered.df$PROBE_ID),] # Move first occurence to temp df
  colnames(ex.odd.df) <- paste(colnames(ex.ordered.df), "1", sep = ".") # Add .1 to column names
  ex.even.df <- ex.ordered.df[duplicated(ex.ordered.df$PROBE_ID),] # Move duplicated occurence to temp df
  colnames(ex.even.df) <- paste(colnames(ex.ordered.df), "2", sep = ".") # Add .2 to column names
  if(!identical(ex.odd.df$PROBE_ID, ex.even.df$PROBE_ID)) {
    warning("Probe order differs between the 2 temp data frames") # Confirm probe order is identical
  }
  ex2 <- cbind(ex.odd.df[,-(1:2)],ex.even.df[,-(1:2)]) # Join the 2 temp data frames
  ex2 <- ex2[,order(colnames(ex2))] # Sort columns by name
  row.names(ex2) <- ex.odd.df$PROBE_ID.1 # Name rows using unique identifiers
  # Make a vector of the gene identifiers in ex2, using ENSG IDs (ORFs) as in the ScKEGG_pathways genset collection
  ex2.ORF <- ex.odd.df$ORF.1 # Needed for indexing later
  rm(ex.ordered.df); rm(ex.odd.df); rm(ex.even.df) # Clean up
} else { ## If probes aren't duplicated on microarry, then don't rearrange data as above
  dup.no = 1
  ex2 <- ex.ORF.df[,-(1:2)] # Matrix-like data object for lmFit
  row.names(ex2) <- ex.ORF.df$PROBE_ID # Name rows using unique identifiers
  # Make a vector of the gene identifiers in ex2, using ENSG IDs (ORFs) as in the ScKEGG_pathways genset collection
  ex2.ORF <- ex.df$ORF # Vector of corresponding gene identifiers in ex2, using ENSG IDs (ORFs) as in the ScKEGG_pathways genset collection; needed later for indexing gene sets
}
cat(paste(GSEname, "microarray has probes duplicated", dup.no, "times."))
cat(paste("Data matrix now has", dim(ex2)[1], "probes and", dim(ex2)[2], "samples."))
rm(ex); rm(ex.df); rm(ex.ORF.df) # Clean up

# Design matrix for using group-mean parameterization [limma usersguide 8.3]; Common Reference Design
Treatment <- sub("(.*)-\\d?-?[a-z]$","\\1",phenoData(ESet)$title) # Extract short description of each sample
  # This assumes biological replicates whose name name differs by an extension such as -a or -1a
targets.df <- data.frame(AssayName = colnames(ex2), SlideNum = sub("(GSM\\d*)\\.[12]", "\\1", colnames(ex2)), 
               Treatment = rep(make.names(Treatment), each = dup.no))
lev <- levels(targets.df$Treatment)
f <- factor(targets.df$Treatment, levels=lev)
design <- model.matrix(~0+f)
colnames(design) <- lev

# Estimate the correlation between samples analyzed ≥1x [limma usersguide 18.1.9]
cor <- duplicateCorrelation(ex2, design, block = targets.df$SlideNum) # SLOW
cor$consensus

# Linear modeling [limma usersguide 18.1.9]
fit <- lmFit(ex2, design, block = targets.df$SlideNum, correlation = cor$consensus)
contrasts = make.names(unique(targets.df$Treatment))
cont.matrix <- makeContrasts(contrasts = contrasts, levels = design)

## To find individual differentially expressed genes (DEGs), you could perform the following: ##
fit.1 <- contrasts.fit(fit, cont.matrix[,1])
efit.1 <- eBayes(fit.1)
allGeneValues.1 <- topTable(efit.1, number = Inf, sort.by = "none")

## Prepare for gene set enrichment analysis ##
# Convert .gmt file to Gene Set Collection in R
geneset.collection <- getGmt(file.path(GLdir, GLfile), geneIdType=GenenameIdentifier()) 
# Extract a list of the gene sets with their IDs
gene.sets.IDs <- geneIds(geneset.collection)
# Convert gene IDs in each gene set to indexed row position in ex2
index <- ids2indices(gene.sets = gene.sets.IDs, identifiers = ex2.ORF, remove.empty=TRUE)

## GSEA for linear models using rotation tests ##
# Which gene sets are differentially expressed in each treatment? Since >1 treatment, I'll gather romer results in a large list
names.contrast <- colnames(cont.matrix)
rlist <- vector(mode = "list", length = length(names.contrast))
names(rlist) <- names.contrast
for (i in seq_along(names.contrast)) {
  rlist[[i]] <- romer(ex2, index, design, contrast = cont.matrix[,i], block = targets.df$SlideNum, correlation = cor$consensus)
}
# Each element in rlist is a matrix with rows in the original gene.sets order

# Get descriptors of indexed sets from 2nd column of .gmt file
GeneSets.df <- read.delim(file.path(GLdir, GLfile), header = FALSE, stringsAsFactors = FALSE)
GeneSets.Annot <- GeneSets.df[GeneSets.df$V1 %in% names(index),1:2] # Only include pathways in index
colnames(GeneSets.Annot) <- c("Pathway.ID","Pathway.name")
if(!identical(GeneSets.Annot$Pathway.ID,row.names(rlist[[1]]))) {
  warning("Pathway IDs don't match row names") # Confirm probe order is identical
}
# Add names of each gene set using cbind (now or while looping above):
for (i in seq_along(names.contrast)) {
  rlist[[i]] <- cbind(GeneSets.Annot, rlist[[i]])
}

# Each dataframe in rlist can be written to file (now or while looping above):
dir.create(file.path(WD, "romerResults"))
setwd(file.path(WD, "romerResults"))
for (i in seq_along(names.contrast)) {
  newfilename <- paste0("r.", names(rlist[i]))
  write.table(rlist[[i]], file = newfilename, quote = FALSE, sep = "\t", row.names = FALSE) # Rows will be numbered to restore original order if needed
}

## Objective: Compare the 3 romer scores of each gene set to identify condition most similar to PsAvh172 
# Construct a dataframe of romerResults, while converting the p-values in rlist to -log10 scores
len.r <- length(rlist)
rUpScores <- data.frame(Pathway.ID = r.Query$Pathway.ID, Pathway.name = r.Query$Pathway.name, 
                        NGenes = r.Query$NGenes, PsAvh172 = round(-log10(r.Query$Up), digits = 1))
len.base <- length(rUpScores)
for (i in seq_along(names.contrast)) {
  rUpScores[i+len.base] <- round(-log10(rlist[[i]]$Up), digits = 1)
}
colnames(rUpScores)[(1+len.base):(i+len.base)] <- names.contrast # I couldn't add names during looping

rDownScores <- data.frame(Pathway.ID = r.Query$Pathway.ID, Pathway.name = r.Query$Pathway.name, NGenes = r.Query$NGenes, PsAvh172 = round(-log10(r.Query$Down), digits = 1))
len.base <- length(rDownScores)
for (i in seq_along(names.contrast)) {
  rDownScores[i+len.base] <- round(-log10(rlist[[i]]$Up), digits = 1)
}
colnames(rDownScores)[(1+len.base):(i+len.base)] <- names.contrast # I couldn't add names during looping

rMixedScores <- data.frame(Pathway.ID = r.Query$Pathway.ID, Pathway.name = r.Query$Pathway.name, NGenes = r.Query$NGenes, PsAvh172 = round(-log10(r.Query$Down), digits = 1))
len.base <- length(rMixedScores)
for (i in seq_along(names.contrast)) {
  rMixedScores[i+len.base] <- round(-log10(rlist[[i]]$Up), digits = 1)
}
colnames(rMixedScores)[(1+len.base):(i+len.base)] <- names.contrast # I couldn't add names during looping

# Next combine the 3 dataframes into a matrix, maintaining the columns
rAll3Scores <- rbind(rUpScores, rDownScores, rMixedScores) # Puts numeric data into a matrix

## Determine most similar romer profiles ##
# A. Examine similarity (correlation) of romer scores in each column
pcc <- cor(rAll3Scores[,-(1:3)]) # Pearson correlation coefficients
pheatmap(pcc) # A quick visualization using default settings
# B. Visualize and cluster by Mixed gene set scores as a heatmap
mat <- as.matrix(rMixedScores[,-(1:len.base-1)]) # Puts numeric data into a matrix
pheatmap(mat,
         # cellwidth = 15, cellheight = 12, 
         labels_row = rMixedScores$Pathway.name, show_rownames = F, 
         cluster_rows = F, #annotation_row = ann.row,
         main = paste("PsAvh172 &", GSEname), fontsize = 12)
# !! If legend = T & show_rownames = T, then legend displays on top of rownames !!

sig.set <- function(x) {length(which(x>=2))} # A short function to count scores 2 or greater (= p-val <=.01)
no.sig.sets <- apply(mat, 2, sig.set) 
rMixedScores2 <- rbind(rMixedScores, c(rep(NA,times=3),no.sig.sets)) # Add to bottom of dataframe

setwd(WD); save.image()
