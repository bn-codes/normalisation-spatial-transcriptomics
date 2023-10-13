# create a visualisation of variance contribution of spots within gene groups for raw data and all normalisation algorithms

#------Function to create sub-matrices, each one containing data about genes of one group and spots of one group------

make_submatrices <- function(matrix){
  
  list_submatrices <- list()
  
  for (i in 1:number_of_gene_groups){
    for (j in 1:number_of_spot_groups){
      submatrix <- current_matrix[current_matrix$spot_groups_factors_array == j, current_matrix["gene_groups_factors_array",] == i]
      submatrix <- list(submatrix)
      list_submatrices <- append(list_submatrices, submatrix)
    }
  }
  
  return(list_submatrices)
  
}

#-------------Function to count variance contribution of a spot group within a gene group----------------------

count_variance_contribution <- function(current_submatrix){
  
  sub_gene_means <- colMeans(current_submatrix)
  
  variance_contribution <- 0
  
  for (i in 1:ncol(current_submatrix)){
    index_col <- i
    for (j in 1:nrow(current_submatrix)){
      index_row <- j
      sub_value <- current_submatrix[index_row, index_col]
      sub_value <- sub_value - sub_gene_means[index_col]
      sub_value <- sub_value**2
      variance_contribution <- variance_contribution + sub_value
    }
  }
  return(variance_contribution)
}

#-Count and visualise variance contribution of spots within gene groups for raw data and all normalisation algorithms-

number_of_gene_groups <- 6
number_of_spot_groups <- 3

for (i in 1:length(normalise_matrices)){
  current_matrix <- normalise_matrices[[i]]
  current_matrix <- as.matrix(current_matrix)
  rownames(current_matrix) <- names_of_genes
  colnames(current_matrix) <- names_of_spots
  
  current_matrix <- cbind(current_matrix, gene_groups_factors_array)
  current_matrix <- rbind(current_matrix, spot_groups_factors_array)
  
  current_matrix <- t(current_matrix)
  current_matrix <- as.data.frame(current_matrix)
  current_matrix <- current_matrix[rowSums(is.na(current_matrix)) == 0, ]
  current_matrix <- current_matrix[order(current_matrix$spot_groups_factors_array), 
                                   order(current_matrix["gene_groups_factors_array",])]
  
  submatrices <- list()
  submatrices <- append(submatrices, make_submatrices(current_matrix))
  
  variance_contribution_list <- list()
  for (j in 1:length(submatrices)){
    submatrix <- submatrices[[j]]
    submatrix <- as.data.frame(lapply(submatrix, as.numeric))
    variance_contribution <- count_variance_contribution(submatrix)
    variance_contribution_list <- append(variance_contribution_list, variance_contribution)
  }
  
  sub_data_frame <- as.data.frame(matrix(nrow=18, ncol = 3))
  colnames(sub_data_frame) <- c("nor_var_con", "gene_group", "spot_group")
  sub_data_frame$gene_group <- rep(1:6, times=1, each=3)
  sub_data_frame$spot_group <- rep(1:3, times=6, each=1)
  variance_contribution_values <- unname(unlist(variance_contribution_list))
  
  sums_variance_contribution <- c()
  for (j in 0:5){
    value <- variance_contribution_values[1+j*3] + variance_contribution_values[2+j*3] + variance_contribution_values[3+j*3]
    sums_variance_contribution <- append(sums_variance_contribution, value)
  }
  sums_variance_contribution <- rep(sums_variance_contribution, times=1, each=3)
  
  for (j in 1:18){
    sub_data_frame[j, 1] <- variance_contribution_values[j] / sums_variance_contribution[j]
  }
  
  sub_data_frame$spot_group <- as.factor(sub_data_frame$spot_group)  
  
  if (i == 1){
    second_part_of_name <- "data"
  }
  else{
    second_part_of_name <- "normalised data"
  }
  
  graph_raw <- ggplot(sub_data_frame, aes(fill=spot_group, y=nor_var_con, x=gene_group)) + 
    geom_bar(position="fill", stat="identity") + guides(fill = guide_legend(title = "Spot groups")) +
    scale_fill_manual(values=c('#000099', '#006699', '#00CCCC')) 
  
  graph_raw <- graph_raw + xlab("Gene groups") + ylab("% variance contribution") + ggtitle(paste(order_of_normalisation_algorithms[i], second_part_of_name, sep=" "))
  
  ggsave(paste(order_of_normalisation_algorithms[i], "-spot-groups.pdf", sep = ""), plot = graph_raw, 
         path = path_to_output_folder, width = 3.59, height = 3.74, units = "in")
}

  
