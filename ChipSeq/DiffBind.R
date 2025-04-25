####CHIPQC####
## Load libraries
library(ChIPQC)
library(TxDb.Rnorvegicus.UCSC.rn6.refGene)
TxDb.Rnorvegicus.UCSC.rn4.ensGene <- TxDb.Rnorvegicus.UCSC.rn6.refGene # This is bc there is no rn6 annotation in ChIPQC, so we overwrite the rn4 annotation with the rn6 genome annotation and then call rn4 in CHIPQC

# Load samplesheet
samples <- read.csv("samplesheet.csv")
summary(samples)

# Call the ChIPQC
chipObj <- ChIPQC(samples, annotation = "rn4")
# Save the report
ChIPQCreport(
  chipObj,
  reportName = "ChIP QC report: P24016",
  reportFolder = "ChIPQCreport",
  facetBy = c("Factor"),
  colourBy = c("Sample")
)


#### DIFFBIND #####

# Load packages
library(DiffBind)
library(tidyverse)

compute <- FALSE

# Load samplesheet
samples <- read.csv("samplesheet.csv") %>%
  dplyr::select(-Condition) # Remove condition column since we don't need it

# Load the samplesheet to the dba
dbaObj <- dba(sampleSheet = samples)

# Set the analysis as Paired End
dbaObj$config$singleEnd <- FALSE

#
if (compute) {
  dbaObj <- dba.count(
    dbaObj,
    bUseSummarizeOverlaps = TRUE,
    bParallel = TRUE,
    summits = 50
  )
  dba.save(dbaObj, file = "dbaObj", dir = "DiffBind/") # save the dba objects so we don't loose the computation time in case a error is encountered
} else {
  dbaObj <- dba.load(file = "dbaObj", dir = "DiffBind/")
}

# Compute PCA
# Save PCA
jpeg(file = "DiffBind/PCA_T_vs_R.jpeg")
dba.plotPCA(dbaObj, attributes = DBA_FACTOR, label = DBA_ID)
dev.off()

# Heatmap
jpeg(file = "DiffBind/Correlation_Heatmap.jpeg")
dba.plotHeatmap(dbaObj)
dev.off()


# Normalization
dbaObj <- dba.normalize(dbaObj, normalize = DBA_NORM_DEFAULT)

# Contrasts
dbaObj <- dba.contrast(
  dbaObj,
  design = "~Factor",
  contrast = c("Factor", "T", "R")
)

# Analysis
dbaObj <- dba.analyze(
  dbaObj,
  method = DBA_ALL_METHODS, # Test edgeR and DEseq2
  bParallel = TRUE,
  bGreylist = FALSE
)

# Save the contrast summary results in a table
dbaObj.table <- tibble(dba.show(dbaObj, bContrasts = TRUE))

write.csv(dbaObj.table, file = "DBA_samples")


# PCA edger
try(
  {
    jpeg(file = "DiffBind/Correlation_Heatmap.jpeg")
    dba.plotPCA(dbaObj, contrast = 1, label = DBA_ID, method = DBA_EDGER)
    dev.off()
  },
  silent = TRUE
)
######
