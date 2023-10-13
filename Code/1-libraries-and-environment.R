# load R libraries
library(Matrix)
library(Seurat)
library(sctransform)
library(Giotto)
library(scRNAseq)
library(scran)
library(arrow)
library(feather)
library(reshape2)
library(ggplot2)
library(dplyr)
library(scales)
library(zellkonverter)

# check and/or set environment options
options(Seurat.object.assay.version = "v5")
genv_exists = checkGiottoEnvironment()
if(!genv_exists){
  installGiottoEnvironment()
}

# specify path to rds input file
path_to_rds_file <- "namz_sce_final_input_qc.rds"

# load rds input file
rds_file <- readRDS(file = path_to_rds_file)

# save rds input file as h5ad (necessary to upload data into Python)
writeH5AD(rds_file, "rds_file.h5ad")


# extract and download gene and spot names as csv file (necessary to upload data into Python) 
names_of_genes <- rds_file@rowRanges@elementMetadata@listData$ENSEMBL
names_of_spots <- rds_file@colData@rownames
write.csv(names_of_genes, "names_of_genes.csv", row.names = FALSE)
write.csv(names_of_spots, "names_of_spots.csv", row.names = FALSE)
