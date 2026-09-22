# ============================================================
# Reusable Volcano Plot Function
#
# Creates a volcano plot for edgeR or DESeq2 results.
#
# Colors:
#   Red   = upregulated
#   Blue  = downregulated
#   Black = not significant
#
# Y-axis:
#   -log2(p-value or adjusted p-value)
# ============================================================

create_volcano_plot <- function(
    data,
    gene_column = "Gene_ID",
    logfc_column,
    pvalue_column,
    logfc_threshold = 1,
    pvalue_threshold = 0.05,
    top_n = 10,
    bottom_n = 5,
    y_max = NULL,
    title = "Volcano Plot",
    x_label = "Log2 Fold Change",
    y_label = "-Log2 P-Value") {
  
  # ----------------------------------------------------------
  # 1. Prepare the plotting data
  # ----------------------------------------------------------
  
  plot_data <- data |>
    dplyr::mutate(
      genes = .data[[gene_column]],
      logFC = .data[[logfc_column]],
      P.Value = .data[[pvalue_column]],
      
      # Prevent log2(0) from producing an infinite value
      negLog2P = -log2(
        pmax(P.Value, .Machine$double.xmin)
      )
    ) |>
    dplyr::filter(
      !is.na(genes),
      !is.na(logFC),
      !is.na(P.Value)
    )
  
  # ----------------------------------------------------------
  # 2. Classify gene regulation
  # ----------------------------------------------------------
  
  plot_data <- plot_data |>
    dplyr::mutate(
      color_group = dplyr::case_when(
        P.Value < pvalue_threshold &
          logFC > logfc_threshold ~ "up-regulated",
        
        P.Value < pvalue_threshold &
          logFC < -logfc_threshold ~ "down-regulated",
        
        TRUE ~ "non-significant"
      )
    )
  
  # ----------------------------------------------------------
  # 3. Select labels from the positive fold-change side
  # ----------------------------------------------------------
  
  right_labels <- plot_data |>
    dplyr::filter(logFC > 0) |>
    dplyr::arrange(P.Value) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::pull(genes)
  
  # ----------------------------------------------------------
  # 4. Select labels from the negative fold-change side
  # ----------------------------------------------------------
  
  left_labels <- plot_data |>
    dplyr::filter(logFC < 0) |>
    dplyr::arrange(P.Value) |>
    dplyr::slice_head(n = bottom_n) |>
    dplyr::pull(genes)
  
  # ----------------------------------------------------------
  # 5. Identify statistically significant genes
  # ----------------------------------------------------------
  
  significant_labels <- plot_data |>
    dplyr::filter(
      P.Value < pvalue_threshold,
      abs(logFC) > logfc_threshold
    ) |>
    dplyr::pull(genes)
  
  # Combine all selected gene labels
  combined_labels <- unique(
    c(
      right_labels,
      left_labels,
      significant_labels
    )
  )
  
  # Mark genes that will receive text labels
  plot_data <- plot_data |>
    dplyr::mutate(
      label_flag = genes %in% combined_labels
    )
  
  # ----------------------------------------------------------
  # 6. Define plot colors
  # ----------------------------------------------------------
  
  color_palette <- c(
    "up-regulated" = "red",
    "down-regulated" = "blue",
    "non-significant" = "black"
  )
  
  # ----------------------------------------------------------
  # 7. Create the volcano plot
  # ----------------------------------------------------------
  
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = logFC,
      y = negLog2P
    )
  ) +
    
    # Plot each gene
    ggplot2::geom_point(
      ggplot2::aes(color = color_group),
      size = 1
    ) +
    
    # Add gene labels
    ggrepel::geom_text_repel(
      data = plot_data |>
        dplyr::filter(label_flag),
      ggplot2::aes(label = genes),
      size = 3,
      max.overlaps = 10,
      color = "black"
    ) +
    
    # Apply the custom color palette
    ggplot2::scale_color_manual(
      values = color_palette,
      name = "Expression Regulation"
    ) +
    
    # Add the horizontal p-value threshold
    ggplot2::geom_hline(
      yintercept = -log2(pvalue_threshold),
      linetype = "dashed",
      color = "black",
      linewidth = 0.5
    ) +
    
    # Add the vertical fold-change thresholds
    ggplot2::geom_vline(
      xintercept = c(
        -logfc_threshold,
        logfc_threshold
      ),
      linetype = "dashed",
      color = "black",
      linewidth = 0.5
    ) +
    
    # Add plot labels
    ggplot2::labs(
      title = title,
      x = x_label,
      y = y_label
    ) +
    
    # Apply the plot theme
    ggplot2::theme_minimal() +
    
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(
        color = "black"
      ),
      legend.position = "bottom"
    )
  
  # ----------------------------------------------------------
  # 8. Apply a custom y-axis limit when provided
  # ----------------------------------------------------------
  
  if (!is.null(y_max)) {
    p <- p +
      ggplot2::coord_cartesian(
        ylim = c(0, y_max)
      )
  }
  
  # ----------------------------------------------------------
  # 9. Return the completed plot
  # ----------------------------------------------------------
  
  return(p)
}