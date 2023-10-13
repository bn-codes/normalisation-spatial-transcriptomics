# run Giotto, Seurat and scran normalisation algorithms
# load Scanpy and stLearn normalised data obtained in Python 
# load clusters of spots

# save raw matrix obtained from the first input file
raw_matrix <- rds_file@assays@data@listData$counts
rownames(raw_matrix) <- names_of_genes
colnames(raw_matrix) <- names_of_spots

# specify order in which the normalisation algorithms are run
order_of_normalisation_algorithms <- c("Raw", "Seurat", "Giotto", "scran", "Scanpy", "stLearn")

#--------------------------------------run Seurat normalisation algorithm-------------------------------------

seurat_matrix <- as.Seurat(rds_file, counts = "counts", data = NULL)

seurat_normalise_matrix <- SCTransform(seurat_matrix, assay = "originalexp")
seurat_normalise_matrix <- seurat_normalise_matrix@assays$SCT@data
rownames(seurat_normalise_matrix) <- names_of_genes
colnames(seurat_normalise_matrix) <- names_of_spots

#--------------------------------------run Giotto normalisation algorithm-------------------------------------

giotto_normalise_matrix <- normalizeGiotto(createGiottoObject(seurat_matrix@assays$originalexp@data))
giotto_normalise_matrix <- getExpression(giotto_normalise_matrix, values = "normalized")
giotto_normalise_matrix <- giotto_normalise_matrix@exprMat
rownames(giotto_normalise_matrix) <- names_of_genes
colnames(giotto_normalise_matrix) <- names_of_spots

#--------------------------------------run scran normalisation algorithm-------------------------------------

clusters_scran <- quickCluster(rds_file)
scran_normalise_matrix <- computeSumFactors(rds_file, clusters=clusters_scran)
scran_normalise_matrix <- logNormCounts(scran_normalise_matrix)
scran_normalise_matrix <- scran_normalise_matrix@assays@data@listData$logcounts
rownames(scran_normalise_matrix) <- names_of_genes
colnames(scran_normalise_matrix) <- names_of_spots

#--------------------------------------upload Scanpy normalisation data-------------------------------------

scanpy_normalise_matrix <- arrow::read_feather(path_to_scanpy_normalise_file)
scanpy_normalise_matrix <- scanpy_normalise_matrix[,-2210]
scanpy_normalise_matrix <- as.matrix(scanpy_normalise_matrix)
scanpy_normalise_matrix <- as(scanpy_normalise_matrix, "dgCMatrix")
rownames(scanpy_normalise_matrix) <- names_of_genes
colnames(scanpy_normalise_matrix) <- names_of_spots

#--------------------------------------upload stLearn normalisation data-------------------------------------

stlearn_normalise_matrix <- arrow::read_feather(path_to_stlearn_normalise_file)
stlearn_normalise_matrix <- data.matrix(stlearn_normalise_matrix)
stlearn_normalise_matrix <- as(stlearn_normalise_matrix, "dgCMatrix")
rownames(stlearn_normalise_matrix) <- names_of_genes
colnames(stlearn_normalise_matrix) <- names_of_spots

#--------------------------------------make a list of normalised matrices-------------------------------------

normalise_matrices <- list(raw_matrix, seurat_normalise_matrix, giotto_normalise_matrix, scran_normalise_matrix, 
                           scanpy_normalise_matrix, stlearn_normalise_matrix)

#--------------------------------------------load clusters of spots-------------------------------------------

cluster1 <- read.csv(path_to_cluster1, header = FALSE)
cluster2 <- read.csv(path_to_cluster2, header = FALSE)
cluster3 <- read.csv(path_to_cluster3, header = FALSE)
clusters_list <- list(cluster1, cluster2, cluster3)
