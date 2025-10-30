#' Calculate comprehensive gene-level network metrics
#'
#' This function calculates various centrality and connectivity measures for each gene
#' in a network, including degree, weighted connectivity, betweenness centrality,
#' clustering coefficient, and other network properties.
#'
#' @param matrix_data Symmetric adjacency or correlation matrix (gene x gene)
#' @param matrix_name Character string name for the matrix
#' @param threshold Numeric threshold for binary network conversion
#' @param matrix_type Character: "correlation" or "adjacency" 
#' @return Data frame with gene-level network metrics
#' @export
calculate_gene_level_metrics <- function(matrix_data, matrix_name, threshold, matrix_type = "adjacency") {
    
    require(igraph)
    
    cat("Calculating metrics for", matrix_name, "(threshold =", threshold, ")\n")
    
    # Diagnostic information
    cat("  Matrix diagnostics:\n")
    cat("    Matrix dimensions:", dim(matrix_data), "\n")
    cat("    Matrix range:", range(matrix_data, na.rm = TRUE), "\n")
    cat("    Diagonal values (first 5):", diag(matrix_data)[1:5], "\n")
    cat("    Off-diagonal range:", range(matrix_data[upper.tri(matrix_data)], na.rm = TRUE), "\n")
    
    gene_names <- rownames(matrix_data)
    n_genes <- nrow(matrix_data)
    
    # Initialize results data frame
    results <- data.frame(
        gene = gene_names,
        matrix_type = matrix_name,
        stringsAsFactors = FALSE
    )
    
    # 1. Weighted connectivity (sum of edge weights)
    if (matrix_type == "correlation") {
        # For correlation matrices, use absolute values for connectivity
        results$weighted_connectivity <- rowSums(abs(matrix_data), na.rm = TRUE)
    } else {
        # For adjacency matrices, use raw values
        results$weighted_connectivity <- rowSums(matrix_data, na.rm = TRUE)
    }
    
    # 2. Degree (number of connections above threshold)
    if (matrix_type == "correlation") {
        # For correlations, count absolute values above threshold
        binary_matrix <- abs(matrix_data) > threshold
    } else {
        # For adjacency, count values above threshold
        binary_matrix <- matrix_data > threshold
    }
    
    results$degree <- rowSums(binary_matrix, na.rm = TRUE)
    
    # Threshold diagnostics
    cat("  Threshold diagnostics:\n")
    cat("    Threshold:", threshold, "\n")
    cat("    Connections above threshold:", sum(binary_matrix), "\n")
    cat("    Percentage above threshold:", round(100 * mean(binary_matrix, na.rm = TRUE), 2), "%\n")
    cat("    Mean degree:", round(mean(results$degree), 2), "\n")
    
    # 3. Create igraph object for centrality measures
    # Use absolute values and remove self-loops
    if (matrix_type == "correlation") {
        graph_matrix <- abs(matrix_data)
        # For correlation matrices, diagonal should be 1, but we set to 0 for network analysis
        diag(graph_matrix) <- 0
        # Apply threshold
        graph_matrix[graph_matrix <= threshold] <- 0
    } else {
        # For adjacency matrices, use raw values
        graph_matrix <- matrix_data
        # For adjacency matrices, ensure diagonal is 0
        diag(graph_matrix) <- 0
        # Apply threshold
        graph_matrix[graph_matrix <= threshold] <- 0
    }
    
    cat("  After thresholding:\n")
    cat("    Non-zero edges:", sum(graph_matrix > 0), "\n")
    cat("    Network density:", round(sum(graph_matrix > 0) / (n_genes * (n_genes - 1)), 4), "\n")
    
    # Create weighted graph
    g <- graph_from_adjacency_matrix(graph_matrix, 
                                   mode = "undirected", 
                                   weighted = TRUE,
                                   diag = FALSE)
    
    # Ensure vertex names are preserved
    V(g)$name <- gene_names
    
    # Graph diagnostics
    cat("  Graph diagnostics:\n")
    cat("    Vertices:", vcount(g), "\n")
    cat("    Edges:", ecount(g), "\n")
    cat("    Density:", edge_density(g), "\n")
    if (ecount(g) > 0) {
        cat("    Edge weight range:", range(E(g)$weight, na.rm = TRUE), "\n")
    }
    
    # 4. Betweenness centrality
    cat("  Calculating betweenness centrality...\n")
    results$betweenness_centrality <- betweenness(g, weights = E(g)$weight, normalized = TRUE)
    
    # 5. Closeness centrality  
    cat("  Calculating closeness centrality...\n")
    results$closeness_centrality <- closeness(g, weights = E(g)$weight, normalized = TRUE)
    
    # 6. Eigenvector centrality
    cat("  Calculating eigenvector centrality...\n")
    eig_cent <- eigen_centrality(g, weights = E(g)$weight)
    results$eigenvector_centrality <- eig_cent$vector
    
    # 7. Clustering coefficient (local transitivity)
    cat("  Calculating clustering coefficient...\n")
    results$clustering_coefficient <- transitivity(g, type = "local", weights = E(g)$weight)
    
    # 8. PageRank centrality
    cat("  Calculating PageRank centrality...\n") 
    results$pagerank_centrality <- page_rank(g, weights = E(g)$weight)$vector
    
    # 9. Hub and authority scores (for directed interpretation)
    cat("  Calculating hub/authority scores...\n")
    hub_auth <- hits_scores(g, weights = E(g)$weight)
    results$hub_score <- hub_auth$vector
    
    # 10. Mean edge weight (average weight of connections above threshold)
    results$mean_edge_weight <- sapply(1:n_genes, function(i) {
        connected_weights <- graph_matrix[i, graph_matrix[i,] > 0]  # Use cleaned matrix
        if (length(connected_weights) > 0) {
            mean(connected_weights, na.rm = TRUE)
        } else {
            0
        }
    })
    
    # 11. Network efficiency (inverse of average shortest path)
    cat("  Calculating local efficiency...\n")
    results$local_efficiency <- sapply(1:n_genes, function(i) {
        neighbors <- which(binary_matrix[i,])
        if (length(neighbors) < 2) return(0)
        
        subgraph_nodes <- c(i, neighbors)
        subgraph <- induced_subgraph(g, subgraph_nodes)
        
        if (vcount(subgraph) < 2) return(0)
        
        distances <- distances(subgraph, weights = E(subgraph)$weight)
        distances[distances == Inf] <- 0
        distances[distances == 0] <- Inf
        diag(distances) <- 0
        
        if (all(is.infinite(distances))) return(0)
        
        efficiency_sum <- sum(1/distances[distances != Inf], na.rm = TRUE)
        n_pairs <- length(neighbors) * (length(neighbors) - 1)
        
        if (n_pairs > 0) efficiency_sum / n_pairs else 0
    })
    
    # 12. Add some basic statistics (use cleaned graph matrix, not original)
    results$max_edge_weight <- apply(graph_matrix, 1, max, na.rm = TRUE)
    results$min_edge_weight <- apply(graph_matrix, 1, min, na.rm = TRUE)
    results$std_edge_weight <- apply(graph_matrix, 1, sd, na.rm = TRUE)
    
    # Handle any NAs or infinite values
    results[is.na(results)] <- 0
    results[sapply(results, is.infinite)] <- 0
    
    cat("  Metrics calculated for", n_genes, "genes\n")
    
    return(results)
}

#' Print summary of gene-level metrics
#'
#' @param gene_metrics Data frame output from calculate_gene_level_metrics
#' @param matrix_name Name of the matrix for reporting
#' @export
summarize_gene_metrics <- function(gene_metrics, matrix_name) {
    
    cat("\n=== Gene-Level Metrics Summary:", matrix_name, "===\n")
    
    numeric_cols <- sapply(gene_metrics, is.numeric)
    
    for (col in names(gene_metrics)[numeric_cols]) {
        if (col == "gene") next
        
        values <- gene_metrics[[col]]
        cat(sprintf("%-25s: Mean = %8.4f, SD = %8.4f, Range = [%8.4f, %8.4f]\n",
                   col,
                   mean(values, na.rm = TRUE),
                   sd(values, na.rm = TRUE), 
                   min(values, na.rm = TRUE),
                   max(values, na.rm = TRUE)))
    }
    
    # Identify top genes by different metrics
    cat("\nTop 5 genes by weighted connectivity:\n")
    top_connectivity <- gene_metrics[order(-gene_metrics$weighted_connectivity), ][1:5, c("gene", "weighted_connectivity")]
    print(top_connectivity, row.names = FALSE)
    
    cat("\nTop 5 genes by betweenness centrality:\n")  
    top_betweenness <- gene_metrics[order(-gene_metrics$betweenness_centrality), ][1:5, c("gene", "betweenness_centrality")]
    print(top_betweenness, row.names = FALSE)
}