# create a visualisation of (not)-scaled normalised data based on logarithmic spot UMI count
# for raw data and all normalisation algorithms


# part of the code taken and modified from "Hafemeister, C., Satija, R. Normalization and variance stabilization of
# single-cell RNA-seq data using regularized negative binomial regression. Genome Biol 20, 296 (2019).
# https://doi.org/10.1186/s13059-019-1874-1"

#-------------------------------Function to visualise not-scaled data----------------------------------------------

visualise_not_scaled_data <- function(raw_expression, norm_expression, name_algorithm){
  
  cm <- raw_expression
  cm <- cm[rowSums(cm > 0) > 4, ]
  
  gene.gmean <- sctransform:::row_gmean(cm)
  gene.grp.breaks <- seq(from=min(log10(gene.gmean))-10*.Machine$double.eps, 
                         to=max(log10(gene.gmean)), length.out = 7)
  gene.grp <- cut(log10(gene.gmean), breaks=gene.grp.breaks, ordered_result = TRUE)
  gene.grp <- factor(gene.grp, ordered=TRUE, levels=rev(levels(gene.grp)))
  names(gene.grp) <- names(gene.gmean)
  
  norm.expression <- as.matrix(norm_expression)
  
  x <- log10(colSums(cm))
  x.out <- seq(from=quantile(x, 0.05), to=quantile(x, 0.95), length.out = 200)
  bw <- bw.SJ(x) * 20
  norm.expression.smooth <- t(apply(norm.expression, 1, function(y) 
   ksmooth(y = y, x = x, x.points = x.out, bandwidth = bw)$y))
  
  # create a data frame for plotting
  df <- melt(t(norm.expression.smooth), varnames = c('cell', 'gene'), value.name = 'expr')
  # add total UMI per cell to data frame
  df$log.umi <- x.out
  # add gene group label
  df$grp <- gene.grp[df$gene]
  
  # create a data frame for the labels
  labdf <- df %>% 
    group_by(grp) %>% 
    summarise(log.umi = median(log.umi), 
              expr = quantile(df$expr, 0.999), 
              lab = paste('Gene group', as.numeric(grp))[1])
  
  theme_set(theme_classic(base_size=9))
  
  if (name_algorithm == "Raw"){
    second_part_of_name <- "expression"
  }
  else{
    second_part_of_name <- "normalised expression"
  }
  
  plot_new <- ggplot(df, aes(log.umi, expr)) + 
    stat_summary(aes(fill=grp, group=grp), geom="ribbon", 
                 fun.ymin = function(z) { quantile(z,0.25) },
                 fun.ymax = function(z) { quantile(z,0.75) }, 
                 alpha=1.0) +
    scale_fill_manual(values = rev(brewer_pal(palette='YlOrRd')(7))) +
    stat_summary(fun.y = "mean", size = 1.5, geom = "line", color='black') +
    facet_wrap(~ grp, ncol = 2, scales = 'free_y') +
    coord_cartesian(ylim=quantile(df$expr, c(0.001, 0.999)), expand = FALSE) +
    geom_text(data=labdf, aes(label = lab), vjust = 1, size = 2.2, color='gray25') +
    xlab('Total spot UMI [log10]') + ylab(paste(name_algorithm, 'not-scaled', second_part_of_name, sep = " ")) +
    theme(strip.text = element_blank()) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_x_continuous(labels = function(x) {format(x, digits=2)})
  
  return(plot_new)
}

#-------------------------------Function to visualise scaled data----------------------------------------------

visualize_scaled_data <- function(raw_expression, norm_expression, name_algorithm){
  
  cm <- raw_expression
  cm <- cm[rowSums(cm > 0) > 4, ]
  
  gene.gmean <- sctransform:::row_gmean(cm)
  gene.grp.breaks <- seq(from=min(log10(gene.gmean))-10*.Machine$double.eps, 
                         to=max(log10(gene.gmean)), length.out = 7)
  gene.grp <- cut(log10(gene.gmean), breaks=gene.grp.breaks, ordered_result = TRUE)
  gene.grp <- factor(gene.grp, ordered=TRUE, levels=rev(levels(gene.grp)))
  names(gene.grp) <- names(gene.gmean)
  
  norm.expression <- as.matrix(norm_expression)
  
  x <- log10(colSums(cm))
  x.out <- seq(from=quantile(x, 0.05), to=quantile(x, 0.95), length.out = 200)
  bw <- bw.SJ(x) * 20
  norm.expression.smooth <- t(apply(norm.expression, 1, function(y) 
    ksmooth(y = scale(y), x = x, x.points = x.out, bandwidth = bw)$y))
  
  # create a data frame for plotting
  df <- melt(t(norm.expression.smooth), varnames = c('cell', 'gene'), value.name = 'expr')
  # add total UMI per cell to data frame
  df$log.umi <- x.out
  # add gene group label
  df$grp <- gene.grp[df$gene]
  
  # create a data frame for the labels
  labdf <- df %>% 
    group_by(grp) %>% 
    summarise(log.umi = median(log.umi), 
              expr = quantile(df$expr, 0.999), 
              lab = paste('Gene group', as.numeric(grp))[1])
  
  theme_set(theme_classic(base_size=9))
  
  if (name_algorithm == "Raw"){
    second_part_of_name <- "expression"
  }
  else{
    second_part_of_name <- "normalised expression"
  }
  
  plot_new <- ggplot(df, aes(log.umi, expr)) + 
    stat_summary(aes(fill=grp, group=grp), geom="ribbon", 
                 fun.ymin = function(z) { quantile(z,0.25) },
                 fun.ymax = function(z) { quantile(z,0.75) }, 
                 alpha=1.0) +
    scale_fill_manual(values = rev(brewer_pal(palette='YlOrRd')(7))) +
    stat_summary(fun.y = "mean", size = 1.5, geom = "line", color='black') +
    facet_wrap(~ grp, ncol = 2, scales = 'free_y') +
    coord_cartesian(ylim=quantile(df$expr, c(0.001, 0.999)), expand = FALSE) +
    geom_text(data=labdf, aes(label = lab), vjust = 1, size = 2.2, color='gray25') +
    xlab('Total spot UMI [log10]') + ylab(paste(name_algorithm, 'scaled', second_part_of_name, sep = " ")) +
    theme(strip.text = element_blank()) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_x_continuous(labels = function(x) {format(x, digits=2)})
  
  return(plot_new)
}

#---------------------------Run both functions for raw data and all normalisation algorithms--------------------

for (i in 1:length(normalise_matrices)){
  current_matrix <- normalise_matrices[[i]]
  
  plot_not_scaled <- visualise_not_scaled_data(raw_matrix, current_matrix, order_of_normalisation_algorithms[i])
  ggsave(paste(order_of_normalisation_algorithms[i], "-not-scaled.pdf", sep = ""), plot = plot_not_scaled, 
         path = path_to_output_folder, width = 3.59, height = 3.74, units = "in")
  
  plot_scaled <- visualize_scaled_data(raw_matrix, current_matrix, order_of_normalisation_algorithms[i])
  ggsave(paste(order_of_normalisation_algorithms[i], "-scaled.pdf", sep = ""), plot = plot_scaled, 
         path = path_to_output_folder, width = 3.59, height = 3.74, units = "in")
}

