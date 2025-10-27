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
    
    # 3. Create igraph object for centrality measures
    # Use absolute values and remove self-loops
    if (matrix_type == "correlation") {
        graph_matrix <- abs(matrix_data)
    } else {
        graph_matrix <- matrix_data
    }
    
    # Set diagonal to 0 (remove self-loops)
    diag(graph_matrix) <- 0
    
    # Create weighted graph
    g <- graph_from_adjacency_matrix(graph_matrix, 
                                   mode = "undirected", 
                                   weighted = TRUE,
                                   diag = FALSE)
    
    # Ensure vertex names are preserved
    V(g)$name <- gene_names
    
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
    hub_auth <- hub_score(g, weights = E(g)$weight)
    results$hub_score <- hub_auth$vector
    
    # 10. Mean edge weight (average weight of connections)
    results$mean_edge_weight <- sapply(1:n_genes, function(i) {
        connected_weights <- matrix_data[i, matrix_data[i,] > threshold]
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
    
    # 12. Add some basic statistics
    results$max_edge_weight <- apply(matrix_data, 1, max, na.rm = TRUE)
    results$min_edge_weight <- apply(matrix_data, 1, min, na.rm = TRUE)
    results$std_edge_weight <- apply(matrix_data, 1, sd, na.rm = TRUE)
    
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