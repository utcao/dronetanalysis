suppressPackageStartupMessages({
    library(purrr)
    library(data.table)
    library(glue)
    library(stringr)
    library(ggplot2)
    library(tibble)  # for rownames_to_column / column_to_rownames
    library(futile.logger)
})

#' Generate diagnostic plots for choosing connectivity threshold
#' 
#' This function creates multiple plots to help decide on an appropriate
#' connectivity threshold for network analysis. It shows how network properties
#' change across different threshold values.
#'
#' @param matrix_data Symmetric matrix (correlation or adjacency)
#' @param output_file Path to save the PDF plots
#' @param threshold_range Vector of thresholds to test (default: seq(0.05, 0.95, 0.05))
#' @param matrix_name Name for the matrix (used in plot titles)
#' @return Invisible data.frame with threshold analysis results
#' @export

plot_threshold_analysis <- function(matrix_data, output_file, 
                                   threshold_range = seq(0.05, 0.95, 0.05),
                                   matrix_name = "Network") {
    
    # Extract upper triangle values
    upper_tri_values <- matrix_data[upper.tri(matrix_data)]
    upper_tri_values <- upper_tri_values[!is.na(upper_tri_values)]
    
    # Calculate network properties at each threshold
    results <- data.frame(
        threshold = threshold_range,
        n_connections = sapply(threshold_range, function(t) sum(upper_tri_values > t)),
        network_density = sapply(threshold_range, function(t) sum(upper_tri_values > t) / length(upper_tri_values)),
        mean_degree = sapply(threshold_range, function(t) {
            adj_binary <- matrix_data > t
            mean(rowSums(adj_binary, na.rm = TRUE))
        }),
        max_degree = sapply(threshold_range, function(t) {
            adj_binary <- matrix_data > t
            max(rowSums(adj_binary, na.rm = TRUE))
        }),
        mean_edge_weight = sapply(threshold_range, function(t) {
            surviving_edges <- upper_tri_values[upper_tri_values > t]
            if(length(surviving_edges) > 0) mean(surviving_edges) else NA
        })
    )
    
    # Start PDF device
    pdf(file = output_file, width = 16, height = 12)
    par(mfrow = c(3, 2), mar = c(5, 5, 4, 2) + 0.1, cex.main = 1.3, cex.lab = 1.2)
    
    # Plot 1: Edge weight distribution with threshold lines
    hist(upper_tri_values, breaks = 50, main = paste("Edge Weight Distribution -", matrix_name),
         xlab = "Edge Weight", ylab = "Frequency", col = "lightblue", border = "white")
    
    # Add vertical lines for key thresholds
    abline(v = quantile(upper_tri_values, 0.95, na.rm = TRUE), col = "red", lwd = 2, lty = 2)
    abline(v = quantile(upper_tri_values, 0.90, na.rm = TRUE), col = "orange", lwd = 2, lty = 2)
    
    legend("topright", legend = c("95th percentile", "90th percentile"),
           col = c("red", "orange"), lty = 2, lwd = 2, cex = 1.1)

    # Plot 2: Number of connections vs threshold
    plot(results$threshold, results$n_connections, type = "b", pch = 16,
         main = "Number of Connections vs Threshold", 
         xlab = "Threshold", ylab = "Number of Connections",
         col = "darkblue", lwd = 2)
    grid(col = "gray", lty = 3)
    
    # Plot 3: Network density vs threshold
    plot(results$threshold, results$network_density, type = "b", pch = 16,
         main = "Network Density vs Threshold",
         xlab = "Threshold", ylab = "Network Density", 
         col = "darkgreen", lwd = 2)
    grid(col = "gray", lty = 3)
    abline(h = c(0.01, 0.05, 0.10), col = "pink", lty = 2)
    text(0.8, 0.01, "1% density", col = "red", cex = 0.9)
    text(0.8, 0.05, "5% density", col = "red", cex = 0.9)
    text(0.8, 0.10, "10% density", col = "red", cex = 0.9)
    
    # Plot 4: Mean degree vs threshold
    plot(results$threshold, results$mean_degree, type = "b", pch = 16,
         main = "Mean Node Degree vs Threshold",
         xlab = "Threshold", ylab = "Mean Degree",
         col = "purple", lwd = 2)
    grid(col = "gray", lty = 3)
    
    # Plot 5: Mean edge weight of surviving edges
    plot(results$threshold, results$mean_edge_weight, type = "b", pch = 16,
         main = "Mean Weight of Surviving Edges vs Threshold",
         xlab = "Threshold", ylab = "Mean Edge Weight",
         col = "darkorange", lwd = 2)
    grid(col = "gray", lty = 3)
    
    # Plot 6: Log-scale visualization for better resolution
    plot(results$threshold, log10(results$n_connections + 1), type = "b", pch = 16,
         main = "Log10(Connections) vs Threshold",
         xlab = "Threshold", ylab = "Log10(Number of Connections + 1)",
         col = "darkred", lwd = 2)
    grid(col = "gray", lty = 3)
    
    dev.off()
    
    # Print summary statistics
    cat("\n=== Threshold Analysis Results ===\n")
    cat("Matrix:", matrix_name, "\n")
    cat("Data range:", round(range(upper_tri_values, na.rm = TRUE), 3), "\n")
    cat("Total possible connections:", length(upper_tri_values), "\n")
    
    # Suggest some thresholds
    cat("\n=== Suggested Thresholds ===\n")
    
    # Find thresholds for different network densities
    density_1pct <- results$threshold[which.min(abs(results$network_density - 0.01))]
    density_5pct <- results$threshold[which.min(abs(results$network_density - 0.05))]
    density_10pct <- results$threshold[which.min(abs(results$network_density - 0.10))]
    
    cat("For 1% network density: threshold ≈", round(density_1pct, 3), "\n")
    cat("For 5% network density: threshold ≈", round(density_5pct, 3), "\n")
    cat("For 10% network density: threshold ≈", round(density_10pct, 3), "\n")
    
    # Percentile-based suggestions
    cat("95th percentile threshold:", round(quantile(upper_tri_values, 0.95, na.rm = TRUE), 3), "\n")
    cat("90th percentile threshold:", round(quantile(upper_tri_values, 0.90, na.rm = TRUE), 3), "\n")
    
    cat("\nPlots saved to:", output_file, "\n")
    
    return(invisible(results))
}

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

#' Create modules using WGCNA pipeline
#'
#' This function performs the complete WGCNA module detection pipeline including
#' TOM calculation, hierarchical clustering, dynamic tree cutting, and module merging.
#'
#' @param adjacency_matrix A numeric matrix (gene x gene) adjacency matrix
#' @param network_type Character: "signed" or "unsigned" for TOM calculation
#' @param min_module_size Integer: minimum module size (default: 30)
#' @param merge_threshold Numeric: height threshold for merging similar modules (default: 0.25)
#' @param deep_split Integer: sensitivity of module detection (0-4, default: 2)
#' @param output_dir Character: directory to save plots and results
#' @param save_plots Logical: whether to save diagnostic plots (default: TRUE)
#' @return List containing module assignments, eigengenes, and other results
#' @export
create_modules <- function(adjacency_matrix, 
                          network_type = "signed",
                          min_module_size = 30,
                          merge_threshold = 0.25,
                          deep_split = 2,
                          output_dir = NULL,
                          save_plots = TRUE) {
    
    # Validate inputs
    if (!network_type %in% c("signed", "unsigned")) {
        stop("network_type must be 'signed' or 'unsigned'")
    }
    
    if (!is.matrix(adjacency_matrix)) {
        stop("adjacency_matrix must be a matrix")
    }
    
    cat("=== WGCNA Module Creation for", network_type, "network ===\n")
    cat("Matrix dimensions:", dim(adjacency_matrix), "\n")
    
    # Step 1: Calculate TOM and TOM dissimilarity
    cat("Step 1: Calculating TOM similarity matrix...\n")
    tom_matrix <- TOMsimilarity(adjacency_matrix, TOMType = network_type)
    dimnames(tom_matrix) <- list(rownames(adjacency_matrix), colnames(adjacency_matrix))
    
    tomd_matrix <- 1 - tom_matrix
    dimnames(tomd_matrix) <- list(rownames(adjacency_matrix), colnames(adjacency_matrix))
    
    # Step 2: Hierarchical clustering
    cat("Step 2: Performing hierarchical clustering...\n")
    gene_tree <- hclust(as.dist(tomd_matrix), method = "average")
    
    # Step 3: Dynamic tree cutting for initial modules
    cat("Step 3: Dynamic tree cutting (min module size:", min_module_size, ")...\n")
    dynamic_modules <- cutreeDynamic(dendro = gene_tree, 
                                   distM = tomd_matrix,
                                   deepSplit = deep_split, 
                                   pamRespectsDendro = FALSE,
                                   minClusterSize = min_module_size)
    
    # Convert to colours
    module_colours <- labels2colors(dynamic_modules)
    
    cat("Initial modules detected:", length(unique(dynamic_modules)), "\n")
    cat("Genes in grey module (unassigned):", sum(module_colours == "grey"), "\n")
    
    # Step 4: Calculate module eigengenes
    cat("Step 4: Calculating module eigengenes...\n")
    MEs <- moduleEigengenes(adjacency_matrix, colors = module_colours)$eigengenes

    # Step 5: Calculate module eigengene dissimilarity
    MEDiss <- 1 - cor(MEs)
    METree <- hclust(as.dist(MEDiss), method = "average")
    
    # Step 6: Plot module eigengene clustering if requested
    if (save_plots && !is.null(output_dir)) {
        cat("Step 5: Saving module eigengene clustering plot...\n")
        pdf(file = file.path(output_dir, paste0(network_type, "_module_eigengene_clustering.pdf")), 
            width = 10, height = 6)
        plot(METree, main = paste(stringr::str_to_title(network_type), "Module Eigengene Clustering"), 
             xlab = "", sub = "", cex = 0.9)
        abline(h = merge_threshold, col = "red", lwd = 2)
        text(x = par("usr")[1] + 0.1 * (par("usr")[2] - par("usr")[1]), 
             y = merge_threshold + 0.02, 
             labels = paste("Merge threshold =", merge_threshold), 
             col = "red", cex = 0.8)
        dev.off()
    }
    
    # Step 7: Merge close modules
    cat("Step 6: Merging similar modules (threshold:", merge_threshold, ")...\n")
    merge_result <- mergeCloseModules(adjacency_matrix, 
                                    module_colours, 
                                    cutHeight = merge_threshold, 
                                    verbose = 1)
    
    # Extract final results
    final_colours <- merge_result$colours
    final_MEs <- merge_result$newMEs

    cat("Final modules after merging:", length(unique(final_colours[final_colours != "grey"])), "\n")
    cat("Genes in grey module (unassigned):", sum(final_colours == "grey"), "\n")

    # Return comprehensive results
    results <- list(
        network_type = network_type,
        adjacency_matrix = adjacency_matrix,
        tom_matrix = tom_matrix,
        tomd_matrix = tomd_matrix,
        gene_tree = gene_tree,
        initial_modules = dynamic_modules,
        initial_colours = module_colours,
        initial_MEs = MEs,
        final_colours = final_colours,
        final_MEs = final_MEs,
        merge_result = merge_result,
        parameters = list(
            min_module_size = min_module_size,
            merge_threshold = merge_threshold,
            deep_split = deep_split
        )
    )
    
    cat("Module creation complete!\n\n")
    return(results)
}

#' Identify hub genes from WGCNA modules
#'
#' This function identifies hub genes using multiple approaches: single top hub per module,
#' top N hubs per module by connectivity, and comprehensive connectivity statistics.
#'
#' @param module_results List: output from create_modules() function
#' @param top_n_hubs Integer: number of top hub genes per module (default: 5)
#' @param output_dir Character: directory to save results (optional)
#' @param save_files Logical: whether to save results to CSV files (default: TRUE)
#' @return List containing various hub gene identification results
#' @export
identify_hubs <- function(module_results, 
                         top_n_hubs = 5,
                         output_dir = NULL,
                         save_files = TRUE) {
    
    # Extract necessary components
    adjacency_matrix <- module_results$adjacency_matrix
    final_colours <- module_results$final_colours
    final_MEs <- module_results$final_MEs
    network_type <- module_results$network_type
    
    cat("=== Hub Gene Identification for", network_type, "network ===\n")
    
    # Step 1: Calculate intramodular connectivity
    cat("Step 1: Calculating intramodular connectivity...\n")
    connectivity <- intramodularConnectivity(adjacency_matrix, final_colours)
    connectivity$gene <- rownames(adjacency_matrix)
    connectivity$module <- final_colours
    
    # Step 2: Calculate module membership (kME)
    cat("Step 2: Calculating module membership (kME)...\n")
    MM <- as.data.frame(cor(adjacency_matrix, final_MEs, use = "p"))
    MM$gene <- rownames(adjacency_matrix)
    
    # Step 3: Single top hub per module
    cat("Step 3: Identifying single top hub per module...\n")
    single_hubs <- chooseTopHubInEachModule(adjacency_matrix, final_colours)
    single_hubs_df <- data.frame(
        module = names(single_hubs),
        hub_gene = single_hubs,
        stringsAsFactors = FALSE
    )
    
    # Step 4: Top N hubs per module by connectivity
    cat("Step 4: Identifying top", top_n_hubs, "hubs per module...\n")
    
    # Remove grey module for hub analysis (unassigned genes)
    connectivity_no_grey <- connectivity[connectivity$module != "grey", ]
    
    if (nrow(connectivity_no_grey) > 0) {
        top_hubs <- connectivity_no_grey %>%
            group_by(module) %>%
            arrange(desc(kWithin)) %>%
            slice_head(n = top_n_hubs) %>%
            ungroup()
    } else {
        top_hubs <- data.frame()
        cat("Warning: No genes assigned to modules (all genes in grey module)\n")
    }
    
    # Step 5: Module summary statistics
    cat("Step 5: Calculating module summary statistics...\n")
    module_summary <- connectivity %>%
        filter(module != "grey") %>%
        group_by(module) %>%
        summarise(
            module_size = n(),
            mean_kWithin = mean(kWithin, na.rm = TRUE),
            max_kWithin = max(kWithin, na.rm = TRUE),
            mean_kOut = mean(kOut, na.rm = TRUE),
            mean_kDiff = mean(kDiff, na.rm = TRUE),
            .groups = 'drop'
        )
    
    # Print summary
    cat("\n=== Hub Gene Summary ===\n")
    cat("Total modules (excluding grey):", nrow(module_summary), "\n")
    cat("Single hubs identified:", nrow(single_hubs_df), "\n")
    cat("Top hubs identified:", nrow(top_hubs), "\n")
    
    if (nrow(module_summary) > 0) {
        cat("Module size range:", min(module_summary$module_size), "-", max(module_summary$module_size), "\n")
        print(module_summary)
    }
    
    # Step 6: Save results if requested
    if (save_files && !is.null(output_dir)) {
        cat("\nStep 6: Saving results to files...\n")
        
        write.csv(single_hubs_df, 
                 file.path(output_dir, paste0(network_type, "_single_hub_genes.csv")), 
                 row.names = FALSE)
        
        write.csv(top_hubs, 
                 file.path(output_dir, paste0(network_type, "_top_hub_genes.csv")), 
                 row.names = FALSE)
        
        write.csv(connectivity, 
                 file.path(output_dir, paste0(network_type, "_gene_connectivity.csv")), 
                 row.names = FALSE)
        
        write.csv(MM, 
                 file.path(output_dir, paste0(network_type, "_module_membership.csv")), 
                 row.names = FALSE)
        
        write.csv(module_summary, 
                 file.path(output_dir, paste0(network_type, "_module_summary.csv")), 
                 row.names = FALSE)
        
        cat("Results saved to:", output_dir, "\n")
    }
    
    # Return comprehensive results
    results <- list(
        network_type = network_type,
        connectivity = connectivity,
        module_membership = MM,
        single_hubs = single_hubs_df,
        top_hubs = top_hubs,
        module_summary = module_summary,
        parameters = list(
            top_n_hubs = top_n_hubs
        )
    )
    
    cat("Hub gene identification complete!\n\n")
    return(results)
}
