# CONCEPT-01: Network Analysis Hierarchy

> **Purpose:** Explain the conceptual hierarchy from raw expression data to biological
> interpretation, clarifying the role of each analytical layer and where the new
> co-expression network metrics pipeline fits.

---

## The Six-Level Hierarchy

```
Level 0 ─ Gene Expression Matrix
Level 1 ─ Co-expression Network  (what the network metrics pipeline builds)
Level 2 ─ Network Metrics         (what the network metrics pipeline computes)
Level 3 ─ Module Detection        (NEXT: community structure algorithms)
Level 4 ─ Eigenvalue-based Methods (spectral, WGCNA, Laplacian)
Level 5 ─ Gene Module Membership   (output of levels 3-4)
Level 6 ─ Biological Interpretation (hubs, pathways, TF enrichment)
```

Each level depends on all levels below it. Levels 3-6 are **not yet implemented**
and are described here as a roadmap.

---

## Level 0: Gene Expression Matrix

**What it is:**
A matrix of shape `(n_genes, n_samples)` containing normalised expression values
(e.g. log-CPM, VOOM-normalised). Each entry `expr[g, s]` is the expression level
of gene `g` in sample `s`.

**How it connects to Level 1:**
The expression matrix is the sole input to correlation-based network construction.
Two genes that tend to rise and fall together across samples are co-expressed.

**In this project:** stored as `expression.h5`, produced by `00convert_expr_to_hdf5.py`.

---

## Level 1: Co-expression Network

**What it is:**
An undirected weighted graph `G = (V, E)` where:
- **Nodes `V`** = genes (|V| = n_genes)
- **Edges `E`** = statistically significant pairwise co-expression relationships
- **Edge weight** = Spearman correlation coefficient `r ∈ [-1, 1]`

A positive edge (`r > 0`) means the two genes tend to increase/decrease together.
A negative edge (`r < 0`) means one increases when the other decreases (anti-correlation).

**How edges are defined (in this pipeline):**
1. Compute all-vs-all Spearman correlations on the full expression matrix.
2. Convert correlations to p-values via the t-distribution (df = n-2).
3. Apply Benjamini-Hochberg FDR correction globally.
4. Retain edges with q-value < α (default 0.05) and |r| ≥ threshold (default 0.0001).

**Relationship to the differential co-expression pipeline:**
The differential pipeline builds THREE networks (Low / High / Differential) and asks
*what changes between conditions*. The network metrics pipeline builds ONE network
from all samples and asks *what structure exists*.

**Key distinction:**
| Question | Pipeline | Method |
|----------|----------|--------|
| What co-expression structure exists? | Network metrics (this pipeline) | Single all-vs-all Spearman + one FDR |
| What changes between low/high? | Differential (existing pipeline) | Two subsets + Fisher-Z test + 3× FDR |

---

## Level 2: Network Metrics

**What they are:**
Quantitative summaries of the graph topology — how it is organised, which genes are
central, and how genes relate to their local neighbourhood.

Metrics fall into three categories:

### 2a. Global Metrics (whole-network properties)

| Metric | Definition | Biological meaning |
|--------|-----------|-------------------|
| `n_edges`, `density` | Edges / possible edges | How connected the network is overall |
| `avg_degree`, `max_degree` | Mean/max edges per gene | Typical vs hub connectivity |
| `power_law_exponent γ` | Slope of log(degree) vs log(P(degree)) | Scale-free if γ ∈ [2, 3]: hubs exist |
| `global_clustering_coefficient` | Mean fraction of neighbour pairs that are connected | Tendency to form local modules |
| `n_components` | Number of disconnected subgraphs | Network fragmentation |
| `assortativity` | Pearson r between degrees at edge endpoints | Positive: hubs connect to hubs; negative: hubs connect to periphery |
| `positive_edge_fraction` | Fraction of edges with r > 0 | Activation-dominated vs mixed network |

### 2b. Per-Gene Centrality Metrics

| Metric | Definition | Biological meaning |
|--------|-----------|-------------------|
| `degree` | Number of co-expression partners | Gene's direct influence radius |
| `weighted_degree` | Sum of |r| values | WGCNA-style connectivity (k) |
| `mean_abs_corr` | Mean |r| of edges | Average co-expression strength |
| `eigenvector_centrality` | Eigenvector of the largest eigenvalue | Importance through important neighbours (hub-of-hubs) |
| `approx_betweenness` | Fraction of shortest paths passing through gene | Bridge / bottleneck genes |
| `pagerank` | Random-walk steady-state probability | Influence considering the full network topology |
| `local_clustering_coefficient` | Fraction of the gene's partner pairs that are connected | Gene embedded in a tight cluster vs sparse hub |
| `positive_fraction` | Fraction of edges with r > 0 | Activator vs repressor profile |

### 2c. Topological Distance Metrics (L1/L2 Neighbourhood)

These measure how the neighbourhood around a gene is structured at two hops:

**Layer definitions (single network):**
- **L1 (direct)**: Genes with a significant edge to the focus gene. Count = `degree`.
- **L2 (second-hop)**: Partners of L1 partners that are NOT the focus gene and NOT in L1.
- **L1→L1 edges**: Co-expression edges within the L1 neighbourhood.
- **L1→L2 edges**: Co-expression edges from L1 nodes to L2 nodes.

| Metric | Formula | Interpretation |
|--------|---------|----------------|
| `L1_n_nodes` | = degree | Size of the immediate neighbourhood |
| `L1_conn_sum` | Σ |r| over L1 edges | Total co-expression strength to direct partners |
| `L1_clique_density` | n_L1-L1 edges / (L1*(L1-1)/2) | How clique-like the direct neighbourhood is (1.0 = all L1 genes are connected to each other) |
| `L2_n_nodes` | # pure 2nd-hop nodes | How far the gene's influence reaches |
| `L2_n_edges` | L1-L1 edges + L1-L2 edges | Connectivity of the outer ring |
| `L2L1_deg` | L2_n_nodes / L1_n_nodes | **Expansion ratio**: < 1 means tightly clustered hub; > 1 means sparse hub in a large neighbourhood |
| `L2L1_conn` | L2_conn_sum / L1_conn_sum | Relative connectivity strength: outer ring vs inner ring |

**Relationship to differential pipeline L2/L1 metrics:**

In the differential pipeline, `L2L1_rewire = L2_rewire / L1_rewire` measures the
*ratio of rewiring events* (how much more or less the outer ring changes than the inner
ring between conditions). This metric is not applicable in the single-network pipeline
because there is no comparison condition and no rewiring concept.

In the network metrics pipeline, `L2L1_deg` and `L2L1_conn` measure *topological
expansion* — the same geometric quantities, but interpreted as topology rather than
change. The formula for `L2L1_deg` is identical in both pipelines.

---

## Level 3: Module Detection / Community Structure

> **Status: PLANNED, not yet implemented.**

**What it is:**
Partition the network nodes (genes) into groups ("modules", "communities", "clusters")
such that genes within a module are more densely connected to each other than to genes
in other modules.

**Why it matters:**
Co-expressed modules often correspond to biological pathways, protein complexes, or
co-regulated gene sets. Detecting modules is the bridge from "this gene is a hub" to
"this gene is a hub within the glycolysis pathway".

**Prerequisite:** A constructed network (Level 1) and its metrics (Level 2). Specifically,
the degree distribution (Level 2a) informs which algorithm is appropriate (scale-free
networks favour algorithms that handle hubs differently from random networks).

**Common algorithms and their relationships to network metrics:**

| Algorithm | Input needed from Level 2 | Key property |
|-----------|--------------------------|-------------|
| Louvain / Leiden | Adjacency + edge weights | Maximises modularity Q; fast for large networks |
| Infomap | Adjacency + edge weights | Random-walk based; good for overlapping paths |
| Walktrap | Adjacency | Short random walks; captures dense subgraphs |
| Label propagation | Adjacency | Fastest; less stable |
| WGCNA hierarchical | Correlation matrix | Based on topological overlap matrix (TOM); standard in genomics |

---

## Level 4: Eigenvalue-Based Clustering Methods

> **Status: PLANNED, not yet implemented.**

**What they are:**
Methods that use the eigenvalues and eigenvectors of a matrix derived from the network
(adjacency, Laplacian, or TOM) to find structure.

### 4a. Spectral Clustering

1. Build the graph Laplacian: `L = D - A` (where D is the diagonal degree matrix, A the adjacency)
2. Compute the first k eigenvectors of L (corresponding to the k smallest eigenvalues)
3. Embed genes as points in R^k using these eigenvectors
4. Cluster the embedded points with k-means

**Why eigenvalues?** The eigenvalue gap (where values jump) indicates how many modules
exist. The eigenvectors encode community structure algebraically.

**Relationship to Level 2:** The degree matrix (D) uses per-gene degree from Level 2.
The scale-free property (Level 2a) means the Laplacian spectrum has a characteristic shape.

### 4b. WGCNA (Weighted Gene Co-expression Network Analysis)

1. Raise the correlation matrix to a soft-threshold power `β` (chosen to approximate scale-free topology):
   `w_ij = |r_ij|^β`
2. Compute the Topological Overlap Matrix (TOM):
   `TOM_ij = (Σ_k w_ik * w_kj + w_ij) / (min(k_i, k_j) + 1 - w_ij)`
3. Hierarchical clustering on `1 - TOM`
4. Cut the dendrogram to define modules

**How Level 2 connects:** The soft-threshold β is chosen as the smallest β where the
network becomes approximately scale-free (i.e. where `power_law_r_squared > 0.8`). This
can be directly read from `global_topology.power_law_r_squared` computed in Level 2.

### 4c. Laplacian Eigenmaps / Diffusion Maps

More general manifold-learning methods that use the graph Laplacian to embed genes in a
low-dimensional space. The embedding coordinates can then be clustered. Useful when
modules are not well-separated (non-convex shapes in gene expression space).

---

## Level 5: Gene Module Membership

> **Status: PLANNED, not yet implemented.**

**What it is:**
The assignment of each gene to one or more modules, plus quantitative membership scores.

### 5a. Hard Membership (one module per gene)
Each gene belongs to exactly one module (typical output of Louvain, k-means, WGCNA).
- `module_id`: integer module label for each gene
- `module_size`: number of genes in the module

### 5b. Soft Membership (graded membership)
A gene may have partial membership in multiple modules. Important for hub genes that
bridge modules.

In WGCNA this is the **Module Membership (kME)**:
`kME_g = Pearson(expr_g, ME_m)` where `ME_m` is the module eigengene
(first principal component of all genes in module m).

| Metric | Meaning |
|--------|---------|
| `kME` | Correlation of a gene with the module eigengene; high = core member |
| `intramodular_connectivity` | Sum of edge weights to genes within the same module |
| `hub_score` | Rank within module by kME or intramodular connectivity |
| `gateway_score` | Connectivity to other modules (high = bridge gene) |

**How Level 2 connects:** Genes with high `eigenvector_centrality` or `betweenness` at
Level 2 are likely to be module hubs or inter-module bridges at Level 5. The `L2L1_deg`
ratio helps distinguish dense-cluster hubs (low L2L1_deg) from sparse broadcast hubs
(high L2L1_deg) before modules are formally detected.

---

## Level 6: Biological Interpretation

> **Status: PARTIALLY IMPLEMENTED** (pathway enrichment script exists for differential pipeline).

**What it is:**
Translate network topology and module membership into biology.

### 6a. Hub Gene Analysis
Genes with the highest degree / weighted_degree / eigenvector_centrality are **hub genes**.
In biological networks, hubs are typically:
- Master regulators or transcription factors
- Genes in core metabolic pathways
- Disease-associated genes with many interaction partners

### 6b. Module-Level Pathway Enrichment
For each module, test whether its member genes are enriched in GO terms, KEGG pathways,
or other gene sets (Fisher's exact test or GSEA).

**In this project:** `pathway_enrichment_hubs.R` (in the differential pipeline) runs GO/KEGG
enrichment on top/bottom N rewiring hub genes. An equivalent can be run on Level 5 modules.

### 6c. Transcription Factor (TF) Enrichment
Are hub genes or module cores enriched for known TF binding motifs? Methods:
- MotifDb + FIMO for sequence-level enrichment
- ChEA / ENCODE for experimental TF-target data
- GENIE3 / GRNBoost2 for co-expression-based GRN inference

### 6d. Inter-condition Comparison (connecting to the differential pipeline)
Once modules are defined from the full-matrix network (Level 3-5), the differential
co-expression pipeline results can be re-interpreted:
- Which modules show the most rewiring between low and high expression conditions?
- Are hub genes identified in the network metrics pipeline the same genes identified as
  rewiring hubs in the differential pipeline?
- Do modules with high L2L1_rewire (differential) correspond to modules with high
  L2L1_deg (network metrics)?

---

## Summary Table: How Levels Connect

| Level | Output | Used by |
|-------|--------|---------|
| 0: Expression matrix | expression.h5 | Level 1 |
| 1: Network (edges + weights) | corr_significant.h5, network_edges.h5 | Level 2, 3, 4 |
| 2: Per-gene metrics | network_metrics_summary.h5, metrics.tsv | Level 3 (β selection), Level 6 (hub ranking) |
| 3: Module assignments | module_labels.tsv | Level 5, 6 |
| 4: Eigenvectors / TOM | (internal to algorithm) | Level 3 (input to clustering) |
| 5: Module membership | kME_matrix.tsv, hub_genes.tsv | Level 6 |
| 6: Biology | Enrichment tables, GRN | Final output |

---

## What to Implement Next (After Network Metrics Pipeline)

The recommended implementation order:

1. **Network metrics pipeline** ← current work (this document describes the context)
2. **Soft-threshold selection** — use `power_law_r_squared` curve from Level 2 to pick β for WGCNA
3. **WGCNA module detection** — build TOM, hierarchical clustering, define modules
4. **Module eigengenes + kME** — compute soft membership scores
5. **Pathway enrichment per module** — adapt `pathway_enrichment_hubs.R`
6. **Cross-pipeline integration** — compare differential rewiring hubs vs module hubs

---

*Document version: 2026-04-08*
*Related pipelines:*
- *Differential co-expression: `src/pipelines/Snakefile_bootstrap`*
- *Network metrics: `network_metrics/Snakefile`*
