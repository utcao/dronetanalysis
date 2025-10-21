suppressPackageStartupMessages({
    library(purrr)
    library(data.table)
    library(glue)
    library(stringr)
    library(ggplot2)
    library(tibble)  # for rownames_to_column / column_to_rownames
    library(futile.logger)
})

#' Render plot for deciding soft-thresholding power for WGCNA
#' 
#' This function generates and saves a plot to help decide the appropriate soft-thresholding power
#' for constructing a scale-free network using WGCNA. The plot includes the scale-free topology fit index
#' and mean connectivity as functions of the soft-thresholding power.
#'
#' @param sft Output from WGCNA::pickSoftThreshold() function containing fit indices
#' @param output_file Full path to output PDF file
#' @param powers Vector of powers used in soft threshold analysis (optional, extracted from sft if not provided)
#' @return NULL (saves plot to file)
#' @export

sft_plot <- function(sft, output_file, powers = c(1:20)){
    # Extract powers from sft if not provided
    if(is.null(powers)) {
        powers <- sft$fitIndices[,1]
    }
    
    # Create PDF with 2x2 layout
    pdf(file = output_file, width = 12, height = 10)
    par(mfrow = c(2,2), mar = c(5, 5, 4, 2) + 0.1) 
    cex1 = 1.0

    # Scale-free topology fit index as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
        xlab="Soft Threshold (power)",
        ylab="Scale Free Topology Model Fit,signed R^2",
        type="n", main = "Scale independence")
        text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
            labels=powers, cex=cex1, col="red")
        abline(h=0.9, col="red")

        # Mean connectivity as a function of the soft-thresholding power
        plot(sft$fitIndices[,1], sft$fitIndices[,5],
            xlab="Soft Threshold (power)",
            ylab="Mean Connectivity", type="n",
            main = "Mean connectivity")
        text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")

        plot(sft$fitIndices[,1],sft$fitIndices[,6],
            xlab="Soft Threshold (power)",
            ylab="Median Connectivity",type="n",
            main = "Median connectivity")
        text(sft$fitIndices[,1], sft$fitIndices[,6], labels=powers,cex=cex1,col="red")

        plot(sft$fitIndices[,1],sft$fitIndices[,7],
            xlab="Soft Threshold (power)",
            ylab="Max Connectivity",type="n",
            main = "Max connectivity")
        text(sft$fitIndices[,1], sft$fitIndices[,7], labels=powers,cex=cex1,col="red")
        dev.off()
}
#' Calculate network metrics for a given matrix
#'
#' This function computes various network metrics for a provided gene-gene matrix,
#' including network-level statistics (mean, median, max, min correlations) and 
#' node-level statistics (connectivity measures). It uses 0.1 as the threshold for degree calculation,
#' which is appropriate for global modifier identification.
#' 
#' @param matrix_data A numeric matrix (gene x gene) for which to calculate network metrics
#' @param matrix_name A string name for the matrix (used in output messages)
#' @return A list containing various network metrics
#' @export
#'
#' 
calculate_network_metrics <- function(matrix_data, matrix_name) {
    cat("\n=== Network Metrics for", matrix_name, "===\n")
    upper_tri_values <- matrix_data[upper.tri(matrix_data)]
    
    # Network-level metrics
    mean_abs_corr <- mean(abs(upper_tri_values), na.rm = TRUE)
    median_abs_corr <- median(abs(upper_tri_values), na.rm = TRUE)
    max_corr <- max(upper_tri_values, na.rm = TRUE)
    min_corr <- min(upper_tri_values, na.rm = TRUE)
    
    # Sum of correlations (connectivity measure)
    sum_correlations <- rowSums(abs(matrix_data), na.rm = TRUE)
    
    # Weighted connectivity (sum of all edge weights)
    weighted_connectivity <- rowSums(matrix_data, na.rm = TRUE)
    
    # Degree (discrete connections above threshold)
    connection_threshold <- 0.3 ## needs adjusting to appropriate value depending on adjacency type
    degree <- rowSums(matrix_data > connection_threshold, na.rm = TRUE)
    
    # Max connectivity measures
    max_connectivity <- max(weighted_connectivity, na.rm = TRUE)
    max_degree <- max(degree, na.rm = TRUE)
    
    # Print results
    cat("Network-level metrics:\n")
    cat("  Mean absolute correlation:", round(mean_abs_corr, 4), "\n")
    cat("  Median absolute correlation:", round(median_abs_corr, 4), "\n")
    cat("  Max correlation:", round(max_corr, 4), "\n")
    cat("  Min correlation:", round(min_corr, 4), "\n")
    
    cat("Node-level metrics:\n")
    cat("  Mean weighted connectivity:", round(mean(weighted_connectivity, na.rm = TRUE), 4), "\n")
    cat("  Mean degree (threshold=", connection_threshold, "):", round(mean(degree, na.rm = TRUE), 2), "\n")
    cat("  Max weighted connectivity:", round(max_connectivity, 4), "\n")
    cat("  Max degree:", max_degree, "\n")
    
    # Return metrics as a list for potential further analysis
    return(list(
        matrix_name = matrix_name,
        mean_abs_corr = mean_abs_corr,
        median_abs_corr = median_abs_corr,
        max_corr = max_corr,
        min_corr = min_corr,
        weighted_connectivity = weighted_connectivity,
        degree = degree,
        sum_correlations = sum_correlations
    ))
}