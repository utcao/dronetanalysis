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
    
    # 11. Local efficiency: average inverse path length in neighborhood
    # Local efficiency measures how well-connected a node's neighbors are to each other
    cat("  Calculating local efficiency...\n")
    results$local_efficiency <- sapply(1:vcount(g), function(i) {
        # Get neighbors of node i from the actual graph object (not binary_matrix)
        neighbors <- neighbors(g, v = i)
        
        if (length(neighbors) < 2) return(0)  # Need at least 2 neighbors for local efficiency
        
        # Create subgraph of neighbors only (not including node i)
        subgraph <- induced_subgraph(g, neighbors)
        
        if (ecount(subgraph) == 0) return(0)  # If no edges between neighbors, efficiency = 0
        
        # Compute shortest path distances in the subgraph
        dist_matrix <- distances(subgraph, weights = E(subgraph)$weight)
        
        # Convert distances to efficiencies (1/distance)
        # Inf or missing paths should be treated as 0 efficiency
        efficiency_matrix <- matrix(0, nrow = nrow(dist_matrix), ncol = ncol(dist_matrix))
        for (r in 1:nrow(dist_matrix)) {
            for (c in 1:ncol(dist_matrix)) {
                if (r != c && dist_matrix[r,c] > 0 && !is.infinite(dist_matrix[r,c])) {
                    efficiency_matrix[r,c] <- 1 / dist_matrix[r,c]
                }
            }
        }
        
        # Local efficiency = mean pairwise efficiency among neighbors
        n_neighbors <- length(neighbors)
        n_pairs <- n_neighbors * (n_neighbors - 1) / 2  # number of neighbor pairs
        
        if (n_pairs > 0) {
            sum(efficiency_matrix[upper.tri(efficiency_matrix)]) / n_pairs
        } else {
            0
        }
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