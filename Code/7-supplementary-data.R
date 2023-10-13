# create a visualisation of supplementary data 
# 1) graph showing the number of genes in groups
# 2) graph showing the number of spots in groups
# 3) graph showing the variance contribution of spots within the sample

#---------------------------------Graph showing the number of genes in groups-----------------------------------------

graph <- ggplot(as.data.frame(gene_groups_factors_array), aes(x=gene_groups_factors_array, fill=..x..)) + 
  geom_bar(show.legend = FALSE) +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  scale_fill_gradient(low = "black", high = "blue") + 
  ggtitle("Number of genes in gene groups") + xlab("Gene groups") + ylab("Number of genes")

ggsave("number_of_genes.pdf", plot = graph, 
       path = path_to_output_folder, width = 3.59, height = 3.74, units = "in")

#---------------------------------Graph showing the number of spots in groups-----------------------------------------

number_of_spots <- as.data.frame(spot_groups_factors_array)
number_of_spots <- number_of_spots[rowSums(is.na(number_of_spots))==0, ]
number_of_spots <- as.data.frame(number_of_spots)
number_of_spots <- number_of_spots[-nrow(number_of_spots),]
number_of_spots <- as.data.frame(number_of_spots)

graph <- ggplot(number_of_spots, aes(x=number_of_spots, fill=..x..)) + 
  geom_bar(show.legend = FALSE) +
  stat_count(aes(y=..count.., label=..count..), geom="text", vjust=-.5) +
  scale_fill_gradient(low = "black", high = "blue") + 
  ggtitle("Number of spots in spot groups") + xlab("Spot groups") + ylab("Number of spots")

ggsave("number_of_spots.pdf", plot = graph, 
       path = path_to_output_folder, width = 3.59, height = 3.74, units = "in")

#--------------------------Graph showing the variance contribution of spots within the sample-------------------------

number_of_spots_contribution <- rep(0,3)
number_of_spots <- as.array(number_of_spots$number_of_spots)
for (i in number_of_spots){
  i = as.numeric(i)
  if (is.na(i)){
    
  }
  else if (i == "1" | i == "2" | i == "3"){
    number_of_spots_contribution[as.numeric(i)] <- number_of_spots_contribution[as.numeric(i)] + 1
  }
}

number_of_spots_sum <- sum(number_of_spots_contribution)
number_of_spots_contribution <- number_of_spots_contribution / number_of_spots_sum
number_of_spots_contribution <- as.data.frame(number_of_spots_contribution)
group <- rep(1,3)
number_of_spots_contribution <- cbind(number_of_spots_contribution, group)
spot <- rep(1:3)
number_of_spots_contribution <- cbind(number_of_spots_contribution, spot)
number_of_spots_contribution$spot <- as.factor(number_of_spots_contribution$spot)

graph <- ggplot(number_of_spots_contribution, aes(fill=spot, y=number_of_spots_contribution, x=group)) + 
  geom_bar(position="fill", stat="identity") + guides(fill = guide_legend(title = "Spot groups")) +
  scale_fill_manual(values=c('#000099', '#006699', '#00CCCC')) 

graph <- graph + xlab("Sample") + ylab("% variance contribution") + ggtitle("Variance contribution of spots in the sample")

ggsave("sample-spot-groups.pdf", plot = graph, 
       path = path_to_output_folder, width = 3.59, height = 3.74, units = "in")
