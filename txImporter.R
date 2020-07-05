#!/usr/bin/R

##########
# This is the second part of the RNA-Seq processing pipeline
#=========
# AUTHOR: Katie Metz Reed
# UPDATED: 2019-02-07
#=========
# This is what the pipeline does. (Update this at some point) # UPDATE
##########

# Load any required libraries
library(tximport)
library(GenomicFeatures)



#-------------------------------------- I N I T I A L I Z E --------------------------------------#

cat("\n\n")
cat("------------\n")
cat("Initializing\n")
cat("------------\n")
cat("\n\n")


# USER INPUT
# ==========

# Inputs for (1) project name and (2) the path to the relevant transcriptome annotation file
args <- commandArgs(trailingOnly=T)
config <- (args[1])
txannt <- (args[2])
txiout <- (args[3])
bsname <- (args[4])


# Errors for improper input
if (length(args)!=4) 
{
  stop("Please provide (1) a path to a tsv config file containing all samples used, (2) a transcriptome annotation file path, (3) a path to a tximport output directory, and (4) a base name")
}



# READ IN SALMON FILES
# ====================

# Read in config file
config <- read.delim(config, sep = "\t", header=T)

# Find salmon files from config
samples <- config$Name

files <- file.path(txiout, samples, "quant.sf")

cat("Accessing quant files:", files, sep="\n")

names(files) <- samples



#--------------------------------------------- R U N ---------------------------------------------#

cat("----------------\n")
cat("Running tximport\n")
cat("----------------\n")
cat("\n\n\n")


# Make TxDb from GENCODE Annotation file 
txdb <- makeTxDbFromGFF(txannt, format="gtf")
k <- keys(txdb, keytype="GENEID")
df <- select(txdb, keys=k, keytype="GENEID", columns="TXNAME")

# Assign proper TxDb columns to tx2gene object
tx2gene <- df[, 2:1]

# Run tximport
txi <- tximport(files, type='salmon', tx2gene=tx2gene)

# Write resulting txi object to RDS file
rdsFile <- paste0(txiout, "/", bsname, "_txi.rds")
saveRDS(txi, file=rdsFile)
cat(paste0("Saving txi object: ", rdsFile))




