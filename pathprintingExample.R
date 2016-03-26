# Use fingerprints to analyze the ALL dataset
require(ALL)
data(ALL)
annotation(ALL)

# The chip used was the Affymetrix Human Genome U95 Version 2 Array
# The correspending GEO ID is GPL8300

# Extract portion of the expression matrix
ALL.exprs<-exprs(ALL)
ALL.exprs.sub<-ALL.exprs[,1:5]

# Process fingerprints
ALL.fingerprint<-exprs2fingerprint(exprs = ALL.exprs.sub,
                                   platform = "GPL8300",
                                   species = "human",
                                   progressBar = TRUE
)

head(ALL.fingerprint)


####
# Construct consensus fingerprint based on pluripotent records
# Use this consensus to find similar arrays

pluripotent.consensus<-consensusFingerprint(
  GEO.fingerprint.matrix[,pluripotents.frame$GSM], threshold=0.9)

# calculate distance from the pluripotent consensus
geo.pluripotentDistance<-consensusDistance(
  pluripotent.consensus, GEO.fingerprint.matrix)

# plot histograms
par(mfcol = c(2,1), mar = c(0, 4, 4, 2))
geo.pluripotentDistance.hist<-hist(geo.pluripotentDistance[,"distance"],
                                   nclass = 50, xlim = c(0,1), main = "Distance from pluripotent consensus")

par(mar = c(7, 4, 4, 2))
hist(geo.pluripotentDistance[pluripotents.frame$GSM, "distance"],
     breaks = geo.pluripotentDistance.hist$breaks, xlim = c(0,1), 
     main = "", xlab = "above: all GEO, below: pluripotent samples")


# annotate top 100 matches not in original seed with metadata
geo.pluripotentDistance.noSeed<-geo.pluripotentDistance[
  !(rownames(geo.pluripotentDistance) %in% pluripotents.frame$GSM),
  ]

top.noSeed.meta<-GEO.metadata.matrix[
  match(head(rownames(geo.pluripotentDistance.noSeed), 1000),
        GEO.metadata.matrix$GSM),
  ]
head(top.noSeed.meta[,c("GSM", "GPL", "Source")],10)
