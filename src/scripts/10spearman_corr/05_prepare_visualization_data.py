#!/usr/bin/env python3
"""
Stage 5: Prepare Visualization Data for R/igraph.

Exports network data from HDF5 to TSV/GraphML format for visualization
in R using igraph, or in Cytoscape.

Pipeline position
-----------------
Stage 1-4  ... → differential_network.h5
Stage 5    THIS SCRIPT → visualization_data/

Output Structure
----------------
visualization_data/
├── edges_low.tsv         - Low network edges (from, to, r, sig)
├── edges_high.tsv        - High network edges (from, to, r, sig)
├── edges_diff.tsv        - Differential edges (from, to, delta, r_low, r_high, qual)
├── nodes.tsv             - Node attributes (gene_id, degree_low, degree_high, degree_diff, ...)
├── ego_networks/         - GraphML files for top rewiring genes
│   ├── ego_gene_0.graphml
│   └── ...
└── visualize_networks.R  - R script template for visualization

Usage
-----
python 05_prepare_visualization_data.py \\
    --diff-h5 results/differential_network.h5 \\
    --out-dir visualization_data \\
    --top-n 10 \\
    --gene-ids data/gene_ids.txt
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np

try:
    import networkx as nx
    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False


def load_differential_network(h5_path: Path) -> dict:
    """Load differential network from HDF5."""
    print(f"Loading differential network from {h5_path}...")

    with h5py.File(h5_path, "r") as h5:
        # Metadata
        meta = h5["meta"]
        n_genes = meta.attrs["n_genes"]
        n_significant = meta.attrs["n_significant"]

        if n_significant == 0:
            print("  No significant edges in network.")
            if "gene_names" in h5:
                gene_names = [x.decode() if isinstance(x, bytes) else x for x in h5["gene_names"][:]]
            else:
                gene_names = None
            return {
                "n_genes": n_genes,
                "n_edges": 0,
                "gene_i": np.array([], dtype=np.int32),
                "gene_j": np.array([], dtype=np.int32),
                "r_low": np.array([], dtype=np.float32),
                "r_high": np.array([], dtype=np.float32),
                "delta": np.array([], dtype=np.float32),
                "qual_score": np.array([], dtype=np.int32),
                "qual_label": np.array([], dtype="U20"),
                "degrees_low": np.zeros(n_genes, dtype=np.int32),
                "degrees_high": np.zeros(n_genes, dtype=np.int32),
                "degrees_diff": np.zeros(n_genes, dtype=np.int32),
                "focus_gene": None,
                "gene_names": gene_names,
            }

        # Edge data
        edges = h5["edges"]
        gene_i = edges["gene_i"][:]
        gene_j = edges["gene_j"][:]
        r_low = edges["r_low"][:]
        r_high = edges["r_high"][:]
        delta = edges["delta_base"][:]
        qual_score = edges["qual_score"][:]
        qual_label = edges["qual_label"][:]
        if len(qual_label) > 0 and isinstance(qual_label[0], bytes):
            qual_label = np.array([x.decode() for x in qual_label])

        # Topology - degrees from each network
        topo = h5["topology"]
        degrees_low = topo["global_low/degrees"][:]
        degrees_high = topo["global_high/degrees"][:]
        degrees_diff = topo["global_diff/degrees"][:]

        # Gene names (propagated through HDF5 chain)
        if "gene_names" in h5:
            gene_names = [x.decode() if isinstance(x, bytes) else x for x in h5["gene_names"][:]]
        else:
            gene_names = None

        # Focus gene info if available
        focus_gene = None
        if "focus_gene" in h5:
            focus_gene = h5["focus_gene"].attrs["gene_index"]

    print(f"  {n_significant:,} edges loaded")
    if gene_names:
        print(f"  Gene names: {gene_names[0]} ... {gene_names[-1]}")

    return {
        "n_genes": n_genes,
        "n_edges": n_significant,
        "gene_i": gene_i,
        "gene_j": gene_j,
        "r_low": r_low,
        "r_high": r_high,
        "delta": delta,
        "qual_score": qual_score,
        "qual_label": qual_label,
        "degrees_low": degrees_low,
        "degrees_high": degrees_high,
        "degrees_diff": degrees_diff,
        "focus_gene": focus_gene,
        "gene_names": gene_names,
    }


def export_edge_lists(data: dict, out_dir: Path, gene_ids: list = None, corr_threshold: float = 0.0001) -> None:
    """Export edge lists as TSV files."""
    print("\nExporting edge lists...")

    gene_i = data["gene_i"]
    gene_j = data["gene_j"]
    r_low = data["r_low"]
    r_high = data["r_high"]
    delta = data["delta"]
    qual_label = data["qual_label"]

    def get_id(idx):
        return gene_ids[idx] if gene_ids else f"gene_{idx}"

    # Low network edges
    low_path = out_dir / "edges_low.tsv"
    with open(low_path, "w") as f:
        f.write("from\tto\tr\tabs_r\n")
        for i, (gi, gj, r) in enumerate(zip(gene_i, gene_j, r_low)):
            if np.abs(r) >= corr_threshold:
                f.write(f"{get_id(gi)}\t{get_id(gj)}\t{r:.4f}\t{np.abs(r):.4f}\n")
    print(f"  Saved {low_path}")

    # High network edges
    high_path = out_dir / "edges_high.tsv"
    with open(high_path, "w") as f:
        f.write("from\tto\tr\tabs_r\n")
        for i, (gi, gj, r) in enumerate(zip(gene_i, gene_j, r_high)):
            if np.abs(r) >= corr_threshold:
                f.write(f"{get_id(gi)}\t{get_id(gj)}\t{r:.4f}\t{np.abs(r):.4f}\n")
    print(f"  Saved {high_path}")

    # Differential edges (all significant)
    diff_path = out_dir / "edges_diff.tsv"
    with open(diff_path, "w") as f:
        f.write("from\tto\tdelta\tr_low\tr_high\tabs_delta\tqual_label\n")
        for gi, gj, d, rl, rh, ql in zip(gene_i, gene_j, delta, r_low, r_high, qual_label):
            f.write(f"{get_id(gi)}\t{get_id(gj)}\t{d:.4f}\t{rl:.4f}\t{rh:.4f}\t{np.abs(d):.4f}\t{ql}\n")
    print(f"  Saved {diff_path}")


def export_node_attributes(data: dict, out_dir: Path, gene_ids: list = None) -> None:
    """Export node attributes as TSV."""
    print("\nExporting node attributes...")

    node_path = out_dir / "nodes.tsv"
    n_genes = data["n_genes"]

    # Calculate rewiring metrics per gene from edges
    gene_i = data["gene_i"]
    gene_j = data["gene_j"]
    qual_score = data["qual_score"]
    delta = data["delta"]

    n_disappear = np.zeros(n_genes, dtype=np.int32)
    n_new = np.zeros(n_genes, dtype=np.int32)
    n_sign_change = np.zeros(n_genes, dtype=np.int32)
    sum_abs_delta = np.zeros(n_genes, dtype=np.float32)

    for gi, gj, qs, d in zip(gene_i, gene_j, qual_score, delta):
        for g in [gi, gj]:
            if qs == 1:  # QUAL_DISAPPEAR
                n_disappear[g] += 1
            elif qs == 2:  # QUAL_NEW
                n_new[g] += 1
            elif qs == 3:  # QUAL_SIGN_CHANGE
                n_sign_change[g] += 1
            sum_abs_delta[g] += np.abs(d)

    rewiring_score = n_disappear + n_new + n_sign_change

    def get_id(idx):
        return gene_ids[idx] if gene_ids else f"gene_{idx}"

    with open(node_path, "w") as f:
        f.write("gene_id\tgene_idx\tdegree_low\tdegree_high\tdegree_diff\t"
                "n_disappear\tn_new\tn_sign_change\trewiring_score\tsum_abs_delta\n")
        for i in range(n_genes):
            # Only include genes with at least one edge
            if data["degrees_diff"][i] > 0:
                f.write(f"{get_id(i)}\t{i}\t"
                        f"{data['degrees_low'][i]}\t{data['degrees_high'][i]}\t{data['degrees_diff'][i]}\t"
                        f"{n_disappear[i]}\t{n_new[i]}\t{n_sign_change[i]}\t"
                        f"{rewiring_score[i]}\t{sum_abs_delta[i]:.4f}\n")

    print(f"  Saved {node_path}")


def export_ego_networks(
    data: dict,
    out_dir: Path,
    top_n: int = 10,
    gene_ids: list = None,
) -> None:
    """Export ego networks for top rewiring genes as GraphML."""
    if not HAS_NETWORKX:
        print("\nSkipping ego network export (networkx not installed)")
        return

    print(f"\nExporting top {top_n} ego networks as GraphML...")

    ego_dir = out_dir / "ego_networks"
    ego_dir.mkdir(exist_ok=True)

    gene_i = data["gene_i"]
    gene_j = data["gene_j"]
    r_low = data["r_low"]
    r_high = data["r_high"]
    delta = data["delta"]
    qual_label = data["qual_label"]
    degrees_diff = data["degrees_diff"]

    def get_id(idx):
        return gene_ids[idx] if gene_ids else f"gene_{idx}"

    # Find top rewiring genes
    top_genes = np.argsort(degrees_diff)[::-1][:top_n]

    # Build edge dict for fast lookup
    edge_dict = {}
    for idx, (gi, gj) in enumerate(zip(gene_i, gene_j)):
        key = (min(gi, gj), max(gi, gj))
        edge_dict[key] = {
            "r_low": r_low[idx],
            "r_high": r_high[idx],
            "delta": delta[idx],
            "qual_label": qual_label[idx],
        }

    # Build adjacency list
    adj = {i: set() for i in range(data["n_genes"])}
    for gi, gj in zip(gene_i, gene_j):
        adj[gi].add(gj)
        adj[gj].add(gi)

    for ego in top_genes:
        ego_id = get_id(ego)

        # Get 1-hop neighbors
        neighbors = adj[ego]
        nodes = {ego} | neighbors

        # Create graph
        G = nx.Graph()

        # Add nodes with attributes
        for n in nodes:
            G.add_node(
                get_id(n),
                is_ego=(n == ego),
                degree_low=int(data["degrees_low"][n]),
                degree_high=int(data["degrees_high"][n]),
                degree_diff=int(data["degrees_diff"][n]),
            )

        # Add edges
        for ni in nodes:
            for nj in adj[ni]:
                if nj in nodes and ni < nj:
                    key = (ni, nj)
                    if key in edge_dict:
                        e = edge_dict[key]
                        G.add_edge(
                            get_id(ni), get_id(nj),
                            r_low=float(e["r_low"]),
                            r_high=float(e["r_high"]),
                            delta=float(e["delta"]),
                            qual_label=str(e["qual_label"]),
                        )

        # Save as GraphML
        out_path = ego_dir / f"ego_{ego_id}.graphml"
        nx.write_graphml(G, out_path)
        print(f"  Saved {out_path} ({len(G.nodes())} nodes, {len(G.edges())} edges)")


def generate_r_script(out_dir: Path) -> None:
    """Generate R visualization script template."""
    print("\nGenerating R visualization script...")

    r_script = '''#!/usr/bin/env Rscript
# Differential Network Visualization with igraph
# Generated by 05_prepare_visualization_data.py

library(igraph)
library(ggplot2)
library(dplyr)
library(cowplot)

# =============================================================================
# Load Data
# =============================================================================

# Edge lists
edges_low <- read.delim("edges_low.tsv", stringsAsFactors = FALSE)
edges_high <- read.delim("edges_high.tsv", stringsAsFactors = FALSE)
edges_diff <- read.delim("edges_diff.tsv", stringsAsFactors = FALSE)

# Node attributes
nodes <- read.delim("nodes.tsv", stringsAsFactors = FALSE)

cat("Loaded data:\\n")
cat("  Low network:", nrow(edges_low), "edges\\n")
cat("  High network:", nrow(edges_high), "edges\\n")
cat("  Differential:", nrow(edges_diff), "edges\\n")
cat("  Nodes:", nrow(nodes), "genes\\n")

# =============================================================================
# Create Networks
# =============================================================================

# Differential network
g_diff <- graph_from_data_frame(
  edges_diff[, c("from", "to")],
  directed = FALSE,
  vertices = nodes
)
E(g_diff)$delta <- edges_diff$delta
E(g_diff)$qual_label <- edges_diff$qual_label
E(g_diff)$r_low <- edges_diff$r_low
E(g_diff)$r_high <- edges_diff$r_high

# =============================================================================
# Identify Top Rewiring Genes
# =============================================================================

top_rewiring <- nodes %>%
  arrange(desc(rewiring_score)) %>%
  head(20)

cat("\\nTop 20 rewiring genes:\\n")
print(top_rewiring[, c("gene_id", "degree_low", "degree_high", "rewiring_score")])

# =============================================================================
# Visualization: Ego Networks
# =============================================================================

plot_ego_network <- function(gene_id, g, title_prefix = "") {
  # Extract ego network (1-hop neighborhood)
  ego_nodes <- neighbors(g, gene_id)
  ego_subgraph <- induced_subgraph(g, c(gene_id, names(ego_nodes)))

  # Color edges by qualitative change
  qual_colors <- c(
    "unchanged" = "gray80",
    "disappear" = "#d73027",  # red
    "new" = "#1a9850",        # green
    "sign_change" = "#7570b3", # purple
    "strengthen" = "#fdae61", # orange
    "weaken" = "#abd9e9"      # light blue
  )

  edge_qual <- E(ego_subgraph)$qual_label
  edge_colors <- qual_colors[edge_qual]
  edge_colors[is.na(edge_colors)] <- "gray50"

  # Node sizes by degree_diff
  node_sizes <- V(ego_subgraph)$degree_diff / max(V(ego_subgraph)$degree_diff) * 15 + 5

  # Highlight ego node
  node_colors <- rep("lightblue", vcount(ego_subgraph))
  node_colors[V(ego_subgraph)$name == gene_id] <- "red"

  # Plot
  layout <- layout_with_fr(ego_subgraph)
  plot(ego_subgraph,
       layout = layout,
       vertex.size = node_sizes,
       vertex.color = node_colors,
       vertex.label.cex = 0.7,
       edge.color = edge_colors,
       edge.width = abs(E(ego_subgraph)$delta) * 5 + 1,
       main = paste0(title_prefix, gene_id))

  # Legend
  legend("bottomright",
         legend = names(qual_colors),
         col = qual_colors,
         lwd = 2,
         cex = 0.6,
         bty = "n")
}

# Plot top 4 rewiring genes
pdf("ego_networks_top4.pdf", width = 12, height = 12)
par(mfrow = c(2, 2))
for (gene_id in top_rewiring$gene_id[1:4]) {
  if (gene_id %in% V(g_diff)$name) {
    plot_ego_network(gene_id, g_diff, "Ego: ")
  }
}
dev.off()
cat("\\nSaved ego_networks_top4.pdf\\n")

# =============================================================================
# Visualization: Degree Distribution Comparison
# =============================================================================

pdf("degree_distribution.pdf", width = 10, height = 6)

# Prepare data for ggplot
degree_data <- data.frame(
  degree = c(nodes$degree_low, nodes$degree_high),
  network = rep(c("Low", "High"), each = nrow(nodes))
)

p1 <- ggplot(degree_data, aes(x = degree, fill = network)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("Low" = "#3182bd", "High" = "#e6550d")) +
  labs(title = "Degree Distribution: Low vs High Networks",
       x = "Degree", y = "Count") +
  theme_minimal()

# Degree change scatter
p2 <- ggplot(nodes, aes(x = degree_low, y = degree_high, color = rewiring_score)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_viridis_c() +
  labs(title = "Degree Change: Low vs High",
       x = "Degree (Low)", y = "Degree (High)") +
  theme_minimal()

print(plot_grid(p1, p2, ncol = 2))
dev.off()
cat("Saved degree_distribution.pdf\\n")

# =============================================================================
# Visualization: Qualitative Changes Summary
# =============================================================================

pdf("qualitative_changes.pdf", width = 8, height = 6)

qual_summary <- edges_diff %>%
  group_by(qual_label) %>%
  summarize(count = n(), mean_delta = mean(abs(delta)))

p <- ggplot(qual_summary, aes(x = reorder(qual_label, -count), y = count, fill = qual_label)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "unchanged" = "gray80",
    "disappear" = "#d73027",
    "new" = "#1a9850",
    "sign_change" = "#7570b3",
    "strengthen" = "#fdae61",
    "weaken" = "#abd9e9"
  )) +
  labs(title = "Qualitative Edge Changes (Low -> High)",
       x = "Change Type", y = "Number of Edges") +
  theme_minimal() +
  theme(legend.position = "none")

print(p)
dev.off()
cat("Saved qualitative_changes.pdf\\n")

# =============================================================================
# Save Summary Statistics
# =============================================================================

summary_stats <- list(
  n_edges_low = nrow(edges_low),
  n_edges_high = nrow(edges_high),
  n_edges_diff = nrow(edges_diff),
  n_nodes = nrow(nodes),
  mean_degree_low = mean(nodes$degree_low),
  mean_degree_high = mean(nodes$degree_high),
  n_disappear = sum(edges_diff$qual_label == "disappear"),
  n_new = sum(edges_diff$qual_label == "new"),
  n_sign_change = sum(edges_diff$qual_label == "sign_change")
)

cat("\\n=== Summary Statistics ===\\n")
for (name in names(summary_stats)) {
  cat(sprintf("  %s: %s\\n", name, summary_stats[[name]]))
}

cat("\\nVisualization complete!\\n")
'''

    r_path = out_dir / "visualize_networks.R"
    with open(r_path, "w") as f:
        f.write(r_script)

    print(f"  Saved {r_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stage 5: Prepare visualization data for R/igraph."
    )
    parser.add_argument(
        "--diff-h5", type=str, required=True,
        help="Path to differential_network.h5 from Stage 3.",
    )
    parser.add_argument(
        "--out-dir", type=str, default="visualization_data",
        help="Output directory for visualization files.",
    )
    parser.add_argument(
        "--top-n", type=int, default=10,
        help="Number of top rewiring genes for ego network export.",
    )
    parser.add_argument(
        "--gene-ids", type=str, default=None,
        help="Path to file with gene IDs (one per line).",
    )
    parser.add_argument(
        "--corr-threshold", type=float, default=0.0001,
        help="Minimum |r| for low/high network edges (default: 0.0001).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load differential network
    data = load_differential_network(Path(args.diff_h5))

    # Gene IDs: prefer propagated names from HDF5 chain, fall back to --gene-ids file
    gene_ids = data.get("gene_names")
    if gene_ids is None and args.gene_ids:
        with open(args.gene_ids) as f:
            gene_ids = [line.strip() for line in f]

    if data["n_edges"] == 0:
        print("\nNo significant edges to visualize. Writing empty outputs.")
        # Still generate the R script template so the directory structure is complete
        generate_r_script(out_dir)
        print(f"\nDone. Output directory: {out_dir}")
        return

    # Export edge lists
    export_edge_lists(data, out_dir, gene_ids, args.corr_threshold)

    # Export node attributes
    export_node_attributes(data, out_dir, gene_ids)

    # Export ego networks
    export_ego_networks(data, out_dir, args.top_n, gene_ids)

    # Generate R script
    generate_r_script(out_dir)

    print(f"\n{'=' * 60}")
    print("VISUALIZATION DATA EXPORT COMPLETE")
    print(f"{'=' * 60}")
    print(f"\nOutput directory: {out_dir}")
    print(f"\nTo visualize in R:")
    print(f"  cd {out_dir}")
    print(f"  Rscript visualize_networks.R")


if __name__ == "__main__":
    main()
