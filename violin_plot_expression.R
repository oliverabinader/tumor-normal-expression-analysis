library(dplyr)
library(tidyr)
library(ggplot2)

# --------------------------------------------
# Prepare expression data for violin plotting
# --------------------------------------------

# Purpose: Format gene expression data (normal vs cancer) into a long format suitable for violin plots
prepare_violin_data <- function(normal_data, cancer_data) {

  # Process the normal dataset
  normal_df <- normal_data
  rownames(normal_df) <- normal_df$gene_symbol
  normal_df <- normal_df[, -1, drop = FALSE]
  normal_df <- t(normal_df) %>% as.data.frame(stringsAsFactors = FALSE)
  # label samples as normal
  normal_df$condition <- "normal"

  # Process the cancer dataset
  cancer_df <- cancer_data
  rownames(cancer_df) <- cancer_df$gene_symbol
  cancer_df <- cancer_df[, -1, drop = FALSE]
  cancer_df <- t(cancer_df) %>% as.data.frame(stringsAsFactors = FALSE)
  # label samples as cancer
  cancer_df$condition <- "cancer"

  # Combine both datasets into one
  combined_data <- bind_rows(normal_df, cancer_df)

  # Identify gene columns (everything except the condition column in the dataset)
  gene_cols <- setdiff(colnames(combined_data), "condition")

  # Convert expression values to numeric
  combined_data[gene_cols] <- lapply(combined_data[gene_cols], function(x) as.numeric(as.character(x)))

  # Reshape from wide to long format so that each row is one gene/sample pair
  counts_long <- combined_data %>%
    pivot_longer(
      cols = -condition, # all gene columns
      names_to = "gene", # new column for gene names
      values_to = "expression" # new column for expression values
    )

  return(counts_long)
}

# --------------------------------------------
# Generate violin plot
# --------------------------------------------

# Purpose: Produce a violin plot showing expression distributions for each gene across two conditions (normal vs cancer)
plot_violin_expression <- function(counts_long) {
  p <- ggplot(counts_long, aes(x = gene, y = expression, fill = condition)) +
    geom_violin(trim = FALSE, na.rm = FALSE) +
  
    # Add median bars for each group
    stat_summary(
      fun = median,
      geom = "crossbar",
      position = position_dodge(width = 0.9)
    ) +
    theme_bw() +

    # Format axis text and legend
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 22),
      axis.text.y = element_text(size = 22),
      axis.title.x = element_text(size = 22),
      axis.title.y = element_text(size = 22),
      legend.text = element_text(size = 22),
      legend.title = element_text(size = 22)
    ) +

    # Label axes 
    labs(
      x = "Genes",
      y = "Normalized Expression (log2)"
    ) +

    # Fill colors for the two sample types
    scale_fill_manual(values = c("normal" = "skyblue", "cancer" = "salmon")) +
    scale_y_continuous(
      breaks = round(
        seq(
          min(counts_long$expression, na.rm = TRUE),
          max(counts_long$expression, na.rm = TRUE),
          by = 5
        ),
        1
      )
    )

  return(p)
}

# Example usage:
# counts_long <- prepare_violin_data(normal_normalized_data, cancer_normalized_data)
# p1 <- plot_violin_expression(counts_long)
# print(p1)
