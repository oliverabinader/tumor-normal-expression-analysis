library(corrplot)

# --------------------------------------------
# Prepare expression data for correlation analysis
# --------------------------------------------
prepare_correlation_input <- function(expression_data) {
  corr_df <- expression_data # Start from the original expression data
  rownames(corr_df) <- corr_df$gene_symbol
  corr_df <- corr_df[, -1, drop = FALSE]
  corr_df <- t(corr_df)
  corr_df <- as.data.frame(corr_df, stringsAsFactors = FALSE)

  # Convert all columns to numeric
  corr_df[] <- lapply(corr_df, function(x) as.numeric(as.character(x)))
  return(corr_df)
}

# --------------------------------------------
# Compute Pearson correlation matrix
# --------------------------------------------

# Compute Pearson correlation between all columns (genes)
# use = "complete.obs" removes rows with any NA when computing each pair
compute_correlation_matrix <- function(corr_input) {
  cor(corr_input, method = "pearson", use = "complete.obs")
}

# --------------------------------------------
# Plot correlation matrix
# --------------------------------------------

# Visualize the correlation matrix with corrplot
# cl.cex argument controls the size of the color legend labels
plot_correlation_matrix <- function(cor_matrix) {
  corrplot(cor_matrix, cl.cex = 1.2)
}

# --------------------------------------------
# Compute pairwise gene-gene Pearson correlations
# --------------------------------------------
compute_pairwise_correlations <- function(expression_data) {
  genes <- unique(expression_data$gene_symbol)
  results <- list()
  idx <- 1

  # Loop over all unique gene pairs (i < j) to avoid duplicates and self-pairs
  for (i in 1:(length(genes) - 1)) {
    for (j in (i + 1):length(genes)) {
      gene1 <- genes[i]
      gene2 <- genes[j]

      # Extract expression values for gene1 across all samples (drop gene_symbol column)
      gene1_data <- as.numeric(expression_data[expression_data$gene_symbol == gene1, -1, drop = TRUE])

      # Extract expression values for gene2 across all samples
      gene2_data <- as.numeric(expression_data[expression_data$gene_symbol == gene2, -1, drop = TRUE])

      # Perform Pearson correlation test between the two genes
      res <- cor.test(gene1_data, gene2_data, method = "pearson")

      results[[idx]] <- data.frame(
        Gene1 = gene1,
        Gene2 = gene2,
        Correlation = unname(res$estimate), # Numeric correlation value for this gene pair
        P_Value = res$p.value, # P-value
        stringsAsFactors = FALSE
      )

      # Move to the next slot (gene-pair) in the list
      idx <- idx + 1
    }
  }

  # Combine all rows into one data.frame with one row per gene pair
  bind_rows(results)
}

# Example usage:
# corr_input <- prepare_correlation_input(cancer_normalized_data)
# cor_matrix <- compute_correlation_matrix(corr_input)
# plot_correlation_matrix(cor_matrix)
# pairwise_results <- compute_pairwise_correlations(cancer_normalized_data)
