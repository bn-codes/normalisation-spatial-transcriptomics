# divide genes and spots into their respective groups

# part of the code taken and modified from "Hafemeister, C., Satija, R. Normalization and variance stabilization of
# single-cell RNA-seq data using regularized negative binomial regression. Genome Biol 20, 296 (2019).
# https://doi.org/10.1186/s13059-019-1874-1"

#-----------------------------------------------Gene groups----------------------------------------------------------

# divide genes into groups depending on their factors
# factors are determined by the logarithmic geometric mean of genes

gene.gmean <- sctransform:::row_gmean(raw_matrix)
raw_matrix_row_means <- log10(gene.gmean)
gene.grp.breaks <- seq(from=min(log10(gene.gmean))-10*.Machine$double.eps, 
                       to=max(log10(gene.gmean)), length.out = 7)

gene_groups_factors_array <- c(0)

for (i in names_of_genes){
  index <- which(rownames(raw_matrix) == i)
  index_row_mean <- raw_matrix_row_means[index]
  
  if(index_row_mean > gene.grp.breaks[6]){
    gene_groups_factors_array[index] <- 1
  }
  else if(index_row_mean <= gene.grp.breaks[6] & index_row_mean > gene.grp.breaks[5]){
    gene_groups_factors_array[index] <- 2
  }
  else if(index_row_mean <= gene.grp.breaks[5] & index_row_mean > gene.grp.breaks[4]){
    gene_groups_factors_array[index] <- 3
  }
  else if(index_row_mean <= gene.grp.breaks[4] & index_row_mean > gene.grp.breaks[3]){
    gene_groups_factors_array[index] <- 4
  }
  else if(index_row_mean <= gene.grp.breaks[3] & index_row_mean > gene.grp.breaks[2]){
    gene_groups_factors_array[index] <- 5
  }
  else if(index_row_mean <= gene.grp.breaks[2]){
    gene_groups_factors_array[index] <- 6
  }
}

#-------------------------------------------------Spot groups------------------------------------------------------

# divide spots into groups depending on their factors
# factors are determined by which cluster the spot is located in

spot_groups_factors_array <- c()

for (i in names_of_spots){
  index <- which(colnames(raw_matrix) == i)
  if(i %in% clusters_list[[1]]$V1){
    spot_groups_factors_array[index] <- 1
  }
  else if(i %in% clusters_list[[2]]$V1){
    spot_groups_factors_array[index] <- 2
  }
  else if(i %in% clusters_list[[3]]$V1){
    spot_groups_factors_array[index] <- 3
  }
}

spot_groups_factors_array <- append(spot_groups_factors_array, "xxx")
