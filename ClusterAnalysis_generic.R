Intro <- 
"Generic workflow to identify the yeast transcriptome in a GEO data set
most similar to a query transcriptome, e.g., PsAvh172-overexpressing yeast cells.
This workflow assumes GEO data set is a yeast microarray study with biological replicates 
and can account for duplicated probes (spots) on a microarray."
cat(Intro)

# !! REPLACE THESE PATHS AND FILE NAMES AS NEEDED!!
WORKdir <- ("~/P001DrMorgan") # REPLACE with your work directory path
GSEname <- "GSE54539" # REPLACE with the GSE number of the data set to download from GEO (no spaces in name!)
GLdir <- "~/P001DrMorganXtra/LabNotebooks/GeneListsYeast" # REPLACE with your gene list directory path
GLfile <- "Saccharomyces cerevisiae-KEGG_Pathway_Windows.gmt" # REPLACE with your gene list file name

####################################################################################
## Load previously calculated romerResults from query profile;
## e.g., see RNAseq_workflow_PsAvh172.R
## romer output using default nrot
####################################################################################
r.Query <- read.delim("~/P001DrMorganXtra/r2.PsAvh172",header = TRUE, stringsAsFactors = FALSE)

####################################################################################
library(Biobase)
library(GEOquery)
library(limma)
library(GSEABase)
library(pheatmap)
####################################################################################

# Create and move to a new folder for output files
WD <- file.path(WORKdir, GSEname) # New folder name
dir.create(WD) # Create folder
setwd(WD)

####################################################################################
# Download desired data set from GEO
####################################################################################
# This produces a "Large list" with the desired ExpressionSet as the last element
gse <- getGEO(GSEname, GSEMatrix = TRUE) 

####################################################################################
# load series and platform data from GEO
####################################################################################
ESet <- gse[[1]]
# Extract matrix of M values
ex <- exprs(gse[[1]]) 
cat(paste(GSEname, "data matrix initially has", dim(ex)[1],
                   "probes and", dim(ex)[2], "samples. "))

## convert missing data values (NA) to 0
if (anyNA(ex)) {
  ex[which(is.na(ex))] <- 0
}

####################################################################################
# log2 transform expression values, if not already
####################################################################################

## examines distribution of values
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
LogC <- (qx[5] > 100) ||
  (qx[6] - qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2) # <---- WHAT IS THIS!!?? Geo2R

## if LogC TRUE (not yet transformed), then take log2 of each positive
## value and replace original matrix in the ExpressionSet
if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gse) <- log2(ex) 
} 
rm(LogC); rm(qx)  # Clean up

## uncomment following to see more information about each sample 
# head(phenoData(ESet)$geo_accession)

# GSM number of each sample
# Long name of each sample
head(levels(phenoData(ESet)$title))

####################################################################################
# extract information about features from ESet
####################################################################################
# Put featureData in its own annotated dataframe
features <- featureData(ESet) 
# Gives structure of featureData
# head(varMetadata(features)) 
# head(features$ID) # ID of the features
# head(levels(features$ORF)) # SGD gene names

## Construct data frame with information needed for downstream analysis
# PROBE_ID is now standard, but may need to use ID (or other) for older platforms
ex.df <- data.frame(PROBE_ID = features$PROBE_ID, ORF = features$ORF, ex)
# Only keep rows with a named ORF; removes spot controls, etc.
ex.ORF.df <- ex.df[which(ex.df$ORF != ""),] 

## If each probe is duplicated on microarry, then rearrange data frame to
## treat duplicate spots as technical replicates
if (anyDuplicated(ex.ORF.df$PROBE_ID) > 0) {
  # Warn if only a subset of probes are duplicated
  if (length(ex.ORF.df$ORF)/length(levels(ex.ORF.df$ORF)) < 2) {
    warning("Is the complete probe set duplicated?") 
  }
  dup.no = 2
  
  ## Sort by probe id so duplicate probes are in consecutive rows
  ex.ordered.df <- ex.ORF.df[order(ex.ORF.df$PROBE_ID),]
  
  ## Move values from duplicate spots to duplicate columns (to BLOCK in romer)
  # Move first occurence to temp df
  ex.odd.df <- ex.ordered.df[!duplicated(ex.ordered.df$PROBE_ID),]   
  # Add .1 to column names
  colnames(ex.odd.df) <- paste(colnames(ex.ordered.df), "1", sep = ".")  
  # Move duplicated occurence to temp df
  ex.even.df <- ex.ordered.df[duplicated(ex.ordered.df$PROBE_ID),]   
  # Add .2 to column names
  colnames(ex.even.df) <- paste(colnames(ex.ordered.df), "2", sep = ".")  
  
  ## Confirm probe order is identical
  if (!identical(ex.odd.df$PROBE_ID, ex.even.df$PROBE_ID)) {
    warning("Probe order differs between the 2 temp data frames")         
  }
  
  # Join the 2 temp data frames
  ex2 <- cbind(ex.odd.df[,-(1:2)],ex.even.df[,-(1:2)]) 
  # Sort columns by name
  ex2 <- ex2[,order(colnames(ex2))]      
  # Name rows using unique identifiers
  row.names(ex2) <- ex.odd.df$PROBE_ID.1                                  
  
  ####################################################################################
  # Make a vector of the gene identifiers in ex2, using ENSG IDs
  # (ORFs) as in the ScKEGG_pathways genset collection
  ####################################################################################
  # Needed for indexing later
  ex2.ORF <- ex.odd.df$ORF.1                                              
  rm(ex.ordered.df); rm(ex.odd.df); rm(ex.even.df)  
  ## If probes aren't duplicated on microarry, then don't rearrange data as above
} else {
  dup.no = 1
  # Matrix-like data object for lmFit
  ex2 <- ex.ORF.df[,-(1:2)]    
  # Name rows using unique identifiers
  row.names(ex2) <- ex.ORF.df$PROBE_ID                                    
  
  ####################################################################################
  # Make a vector of the gene identifiers in ex2, using ENSG IDs (ORFs) as in
  # the ScKEGG_pathways genset collection, needed later for indexing gene sets
  ####################################################################################
  ex2.ORF <- ex.df$ORF
}
cat(paste(GSEname, "microarray has probes duplicated", dup.no, "times."))
cat(paste("Data matrix now has", dim(ex2)[1], "probes and", dim(ex2)[2], "samples."))
rm(ex); rm(ex.df); rm(ex.ORF.df)                                         # Clean up

####################################################################################
# Design matrix for using group-mean parameterization [limma usersguide 8.3]
####################################################################################
# Extract short description of each sample
Treatment <- sub("(.*)-\\d?-?[a-z]$","\\1",phenoData(ESet)$title)        

## This assumes biological replicates whose name differs by an extension such as -a or -1a
targets.df <- data.frame(AssayName = colnames(ex2),
                         SlideNum = sub("(GSM\\d*)\\.[12]", "\\1", colnames(ex2)), 
                         Treatment = rep(make.names(Treatment), each = dup.no))
lev <- levels(targets.df$Treatment)
f <- factor(targets.df$Treatment, levels = lev)
design <- model.matrix(~0 + f)
colnames(design) <- lev

## Estimate the correlation between samples analyzed 1x [limma usersguide 18.1.9]
cor <- duplicateCorrelation(ex2, design, block = targets.df$SlideNum) # SLOW
cor$consensus

## Linear modeling [limma usersguide 18.1.9]
fit <- lmFit(ex2, design, block = targets.df$SlideNum, correlation = cor$consensus)
contrasts = make.names(unique(targets.df$Treatment))
cont.matrix <- makeContrasts(contrasts = contrasts, levels = design)

## To find individual differentially expressed genes (DEGs), you could perform the following:
fit.1 <- contrasts.fit(fit, cont.matrix[,1])
efit.1 <- eBayes(fit.1)
allGeneValues.1 <- topTable(efit.1, number = Inf, sort.by = "none")

####################################################################################
## Prepare for gene set enrichment analysis
####################################################################################

# Convert .gmt file to Gene Set Collection in R
geneset.collection <- getGmt(file.path(GLdir, GLfile), geneIdType = GenenameIdentifier()) 

# Extract a list of the gene sets with their IDs
gene.sets.IDs <- geneIds(geneset.collection)

# Convert gene IDs in each gene set to indexed row position in ex2
index <- ids2indices(gene.sets = gene.sets.IDs, identifiers = ex2.ORF, remove.empty = TRUE)

####################################################################################
## GSEA for linear models using rotation tests
####################################################################################

# Which gene sets are differentially expressed in each treatment?
# Since >1 treatment, I'll gather romer results in a large list
names.contrast <- colnames(cont.matrix)
rlist <- vector(mode = "list", length = length(names.contrast))
names(rlist) <- names.contrast
for (i in seq_along(names.contrast)) {
  rlist[[i]] <- romer(ex2, index, design, contrast = cont.matrix[,i],
                      block = targets.df$SlideNum, correlation = cor$consensus)
}

# Each element in rlist is a matrix with rows in the original gene.sets order

# Get descriptors of indexed sets from 2nd column of .gmt file
GeneSets.df <- read.delim(file.path(GLdir, GLfile), header = FALSE, stringsAsFactors = FALSE)
# Only include pathways in index
GeneSets.Annot <- GeneSets.df[GeneSets.df$V1 %in% names(index),1:2]         
colnames(GeneSets.Annot) <- c("Pathway.ID","Pathway.name")

# Confirm probe order is identical
if (!identical(GeneSets.Annot$Pathway.ID,row.names(rlist[[1]]))) {
  warning("Pathway IDs don't match row names")                              
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
  # Rows will be numbered to restore original order if needed
  write.table(rlist[[i]], file = newfilename, quote = FALSE, sep = "\t", row.names = FALSE)
}

####################################################################################
## Objective: Compare the 3 romer scores of each gene set to
## identify condition most similar to PsAvh172 
####################################################################################

# Construct a dataframe of romerResults, while converting the p-values in rlist to -log10 scores
len.r <- length(rlist)
rUpScores <- data.frame(Pathway.ID = r.Query$Pathway.ID, Pathway.name = r.Query$Pathway.name, 
                        NGenes = r.Query$NGenes, PsAvh172 = round(-log10(r.Query$Up), digits = 1))
len.base <- length(rUpScores)
for (i in seq_along(names.contrast)) {
  rUpScores[i + len.base] <- round(-log10(rlist[[i]]$Up), digits = 1)
}
# I couldn't add names during looping
colnames(rUpScores)[(1 + len.base):(i + len.base)] <- names.contrast        

rDownScores <- data.frame(Pathway.ID = r.Query$Pathway.ID, Pathway.name = r.Query$Pathway.name,
                          NGenes = r.Query$NGenes, PsAvh172 = round(-log10(r.Query$Down), digits = 1))
len.base <- length(rDownScores)
for (i in seq_along(names.contrast)) {
  rDownScores[i + len.base] <- round(-log10(rlist[[i]]$Up), digits = 1)
}
# I couldn't add names during looping
colnames(rDownScores)[(1 + len.base):(i + len.base)] <- names.contrast      

rMixedScores <- data.frame(Pathway.ID = r.Query$Pathway.ID, Pathway.name = r.Query$Pathway.name,
                           NGenes = r.Query$NGenes, PsAvh172 = round(-log10(r.Query$Down), digits = 1))
len.base <- length(rMixedScores)
for (i in seq_along(names.contrast)) {
  rMixedScores[i + len.base] <- round(-log10(rlist[[i]]$Up), digits = 1)
}
# I couldn't add names during looping
colnames(rMixedScores)[(1 + len.base):(i + len.base)] <- names.contrast     

# Next combine the 3 dataframes into a matrix, maintaining the columns
# Puts numeric data into a matrix
rAll3Scores <- rbind(rUpScores, rDownScores, rMixedScores)                  

####################################################################################
## Determine most similar romer profiles
####################################################################################

# A. Examine similarity (correlation) of romer scores in each column
pcc <- cor(rAll3Scores[,-(1:3)]) # Pearson correlation coefficients
pheatmap(pcc) # A quick visualization using default settings

# B. Visualize and cluster by Mixed gene set scores as a heatmap
# Puts numeric data into a matrix
mat <- as.matrix(rMixedScores[,-(1:len.base - 1)])                          
pheatmap(mat,
         # cellwidth = 15, cellheight = 12, 
         labels_row = rMixedScores$Pathway.name, show_rownames = F, 
         cluster_rows = F, #annotation_row = ann.row,
         main = paste("PsAvh172 &", GSEname), fontsize = 12)

# !! If legend = T & show_rownames = T, then legend displays on top of rownames !!

# A short function to count scores 2 or greater (= p-val <=.01)
sig.set <- function(x) {length(which(x >= 2))}                             
no.sig.sets <- apply(mat, 2, sig.set) 
# Add to bottom of dataframe
rMixedScores2 <- rbind(rMixedScores, c(rep(NA,times = 3), no.sig.sets))     

setwd(WORKdir); save.image()
