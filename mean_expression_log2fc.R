# --------------------------------------------
# Calculate mean expression for each gene
# --------------------------------------------

# This function computes the average expression value for each gene across all samples in the provided expression dataset.
calculate_gene_means <- function(expression_data, gene_list) {
  results <- lapply(gene_list, function(gene_of_interest) {
    # Subset the rows corresponding to the current gene
    gene_rows <- expression_data[expression_data$gene_symbol == gene_of_interest, , drop = FALSE]

    # Calculate the mean expression across all samples
    # na.rm = TRUE removes missing values before computing the mean
    mean_value <- rowMeans(gene_rows[, -1, drop = FALSE], na.rm = TRUE)

    # Store the gene symbol and its calculated mean expression
    data.frame(gene_symbol = gene_of_interest, mean_value = mean_value, stringsAsFactors = FALSE)
  })
  # Combine the list of results into a single dataframe
  bind_rows(results)
}

# Example usage:
# gene_means_cancer <- calculate_gene_means(cancer_normalized_data, genes_of_interest)
# gene_means_normal <- calculate_gene_means(normal_normalized_data, genes_of_interest)

# --------------------------------------------
# Calculate log2 fold change between groups
# --------------------------------------------

# This function merges mean expression values from two groups (e.g., cancer vs normal) and 
# calculates the log2 fold change as the difference between their mean expression values.
calculate_log2fc <- function(group1_means, group2_means, group1_name = "cancer", group2_name = "normal") {
  # Merge the two dataframes by gene_symbol so the same genes align. Suffixes are added to distinguish mean values from each group.
  merged_df <- merge(group1_means, group2_means, by = "gene_symbol",
    suffixes = c(paste0("_", group1_name), paste0("_", group2_name))
  )

  # Calculate log2 fold change as the difference in mean expression (assuming the input data is already log2-transformed).
  merged_df$log2_fc <- merged_df[[paste0("mean_value_", group1_name)]] - merged_df[[paste0("mean_value_", group2_name)]]

  # Return the dataframe containing gene means and log2 fold change
  merged_df
}

# Example usage:
# log2fc_results <- calculate_log2fc(gene_means_cancer, gene_means_normal)
