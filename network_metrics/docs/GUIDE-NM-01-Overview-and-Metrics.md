# GUIDE-NM-01: Co-Expression Network Metrics Pipeline

> **Purpose:** Build a single co-expression network from a gene expression matrix and
> compute per-gene and global topology metrics. Separate from the differential
> co-expression pipeline (`src/pipelines/Snakefile_bootstrap`).

---

## 1. Directory Structure

```
network_metrics/            ← top-level, independent pipeline
├── Snakefile               ← workflow definition
├── config/
│   └── example.yaml        ← example configuration
├── docs/
│   ├── GUIDE-NM-01-Overview-and-Metrics.md    ← this file
│   └── PLAN-NM-01-Multi-Dataset-Support.md   ← future: multi-dataset refactor plan
└── scripts/
    ├── 01_compute_correlations.py    ← Stage 1: Spearman + BH-FDR
    ├── 02_build_network.py           ← Stage 2: global topology + CSR adjacency
    ├── 03_per_gene_metrics.py        ← Stage 3: per-gene L1/L2 metrics (per-gene job)
    ├── 03b_collect_metrics.py        ← Stage 3b: aggregate per-gene HDF5s
    └── 04_global_metrics.py          ← Stage 4: eigenvector centrality (optional)
```

---

## 2. What This Pipeline Does

**Input:** Gene expression matrix (`expression.h5` — same format as the differential pipeline)

**Output:**
- `corr_significant.h5` — sparse edge list (gene_i, gene_j, r, q-value)
- `network_edges.h5` — global topology + CSR adjacency for efficient per-gene access
- `per_gene/{gi}_{gene_id}.h5` — per-gene L1/L2 metric files (temporary)
- `network_metrics_summary.h5` — aggregated per-gene metrics
- `network_metrics.tsv` — TSV table for downstream R/Python analysis
- `global_topology.h5` — optional: eigenvector centrality, approximate betweenness

**This is NOT the differential pipeline.** The differential pipeline:
- Splits samples into LOW / HIGH subpopulations
- Builds THREE networks (low / high / differential)
- Asks: *what changes between conditions?*

This pipeline:
- Uses ALL samples as one group
- Builds ONE network
- Asks: *what topology structure exists?*

---

## 3. Pipeline DAG

```
expression.h5
      │
      ▼
[Stage 1]  01_compute_correlations.py
      │        All-vs-all Spearman + single BH-FDR
      ▼
corr_significant.h5  (significant edges only, sparse)
      │
      ▼
[Stage 2]  02_build_network.py
      │        Global topology metrics + CSR adjacency index
      ▼
network_edges.h5
      │
      ├──────────────────────┐
      │  (one job per gene,  │
      │   parallelised by    │
      │   Snakemake)         │
      ▼                      │
[Stage 3]  03_per_gene_metrics.py  (per gene)
      │        L1/L2 neighbourhood metrics
      ▼
per_gene/{gi}_{gene_id}.h5  (×n_genes, marked temp)
      │
      ▼
[Stage 3b]  03b_collect_metrics.py
      │         Aggregate all per-gene files
      ▼
network_metrics_summary.h5
network_metrics.tsv
      │
      ▼  (optional)
[Stage 4]  04_global_metrics.py
               Eigenvector centrality + approx. betweenness
               → global_topology.h5
```

---

## 4. Quick Start

### 4.1 Toy data (local, fast test)

```bash
snakemake -s network_metrics/Snakefile \
    --config toy=true out_dir=tmp/nm_toy \
    -j 4 --use-conda
```

### 4.2 Real data (using existing expression.h5)

```bash
snakemake -s network_metrics/Snakefile \
    --configfile network_metrics/config/example.yaml \
    -j 100 --use-conda
```

### 4.3 SGE cluster

```bash
snakemake -s network_metrics/Snakefile \
    --configfile network_metrics/config/example.yaml \
    --profile config/sge_profile \
    -j 200
```

### 4.4 Dry run

```bash
snakemake -s network_metrics/Snakefile \
    --configfile network_metrics/config/example.yaml -n
```

### 4.5 Use existing expression.h5 (skip Stage 0)

```yaml
# in config YAML:
expr_h5: "path/to/existing/expression.h5"
```

---

## 5. Configuration Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `out_dir` | `"results/network_metrics"` | Output directory |
| `expr_tsv` | `""` | Path to expression TSV (alternative to expr_h5) |
| `expr_h5` | `""` | Path to existing expression.h5 (skips preprocessing) |
| `toy` | `false` | Use random toy data (20 genes, 80 samples) |
| `fdr_alpha` | `0.05` | BH-FDR significance threshold for edges |
| `corr_threshold` | `0.0001` | Minimum |r| to store (filters near-zero edges) |
| `chunk_size` | `null` | Row chunk size for Stage 1 (null = single pass, fastest) |
| `compression_level` | `6` | Gzip compression level (1=fast, 9=max) |
| `gene_subset` | `[]` | Optional: restrict Stage 3 to named genes |
| `skip_global_metrics` | `true` | Disable Stage 4 (eigenvector centrality, betweenness) |
| `compute_betweenness` | `false` | Enable approximate betweenness in Stage 4 |
| `betweenness_n_samples` | `500` | Number of pivot nodes for approximate betweenness |

---

## 6. Resource Estimates

| Stage | Rule | Memory (MB) | Runtime (min) | Threads |
|-------|------|-------------|---------------|---------|
| 1 | `compute_correlations` | 8000 | 60 | 8 |
| 2 | `build_network` | 4000 | 10 | 4 |
| 3 | `per_gene_metrics` | 512 | 5 | 1 |
| 3b | `collect_metrics` | 8000 | 30 | 4 |
| 4 | `global_metrics` | 8000 | 60 | 4 |

**Memory analysis for 18,000 genes × 100 samples (Stage 1):**
- z @ z.T intermediate matrix: 18K × 18K × 4B ≈ 1.2 GB (freed immediately after triu extraction)
- corr_triu: 162M × 4B = 648 MB
- pval_triu: 162M × 4B = 648 MB
- qval_triu: 162M × 4B = 648 MB
- **Peak: ~3.8 GB** → `mem_mb=8000` in Snakefile provides comfortable headroom

**Low-RAM mode** (activate with `chunk_size: 1000` in config):
- Avoids holding the full 18K×18K matrix simultaneously
- Reduces Stage 1 peak RAM to ~2 GB at the cost of ~2× runtime

---

## 7. Output Files Reference

### `corr_significant.h5`
```
meta/     n_genes, n_samples, n_tests, n_significant, fdr_alpha, corr_threshold
edges/    gene_i (int32), gene_j (int32), corr (float32), qval (float32)
gene_names   (n_genes,) str  — dataset, not attribute (avoids 64 KB HDF5 limit)
sample_names (n_samples,) str
```

### `network_edges.h5`
```
meta/     n_genes, n_edges, n_genes_active, density, fdr_alpha, corr_threshold
edges/    gene_i, gene_j, corr, qval   (copied from corr_significant.h5)
gene_names  (n_genes,) str
adjacency/
    indptr   (n_genes+1,) int64   — CSR row pointers
    indices  (2*n_edges,) int32   — neighbour gene indices (sorted within each row)
    data     (2*n_edges,) float32 — correlation values
global_topology/
    (all global metric scalars as HDF5 attributes)
    degrees  (n_genes,) int32
```

**Why CSR adjacency?** Each Stage 3 per-gene job loads only O(degree) bytes from the
HDF5 file, instead of loading all 5–50M edges. For 18K genes with 5M significant edges,
the full edge list is ~60 MB; but a single gene with degree 500 needs only ~6 KB.

### `network_metrics_summary.h5` + `network_metrics.tsv`
```
meta/     n_genes, n_ref_genes, n_edges
gene_names  (n_genes,) str
per_gene/
    gene_index, gene_name
    degree, weighted_degree, mean_abs_corr, max_abs_corr, std_abs_corr
    positive_edge_count, negative_edge_count, positive_fraction
    L1_n_nodes, L1_conn_sum, L1_conn_mean, L1_clique_density, n_l1_to_l1_edges
    local_clustering_coefficient
    L2_n_nodes, L2_n_edges, L2_conn_sum, L2_conn_mean, n_l1_to_l2_edges
    L2L1_deg, L2L1_conn
    full_n_edges, full_conn_sum
```

### `global_topology.h5` (Stage 4, optional)
```
per_gene/
    eigenvector_centrality  (n_genes,) float32
    approx_betweenness      (n_genes,) float32    [if compute_betweenness=true]
    pagerank                (n_genes,) float32
```

---

## 8. Complete Metrics Reference

### 8.1 Global Metrics (Stage 2)

| Metric key | Type | Definition | Biological meaning |
|-----------|------|-----------|-------------------|
| `n_genes` | int | Total genes in dataset | Network size |
| `n_genes_active` | int | Genes with ≥1 significant edge | Active participants |
| `fraction_active` | float | n_genes_active / n_genes | Coverage of the network |
| `n_edges` | int | Significant co-expression edges | Network connectivity |
| `density` | float | 2E / (V(V-1)) | Overall connectedness |
| `avg_degree` | float | 2E / V | Mean connectivity |
| `median_degree` | float | Median of active degrees | Typical hub size |
| `max_degree` | int | Highest degree in network | Top hub gene connectivity |
| `std_degree` | float | Degree standard deviation | Heterogeneity of connectivity |
| `degree_p90/95/99` | float | 90th/95th/99th percentile | Hub threshold candidates |
| `power_law_exponent` γ | float | −slope of log(k) vs log(P(k)) | Scale-free topology if γ ∈ [2, 3.5] |
| `power_law_r_squared` | float | R² of log-log regression | Quality of scale-free fit; also used to choose WGCNA β |
| `is_scale_free` | bool | γ ∈ [2, 3.5] AND R² > 0.8 | Whether network has hub architecture |
| `global_clustering_coefficient` | float | Mean of local clustering | Tendency to form modules |
| `n_components` | int | Number of disconnected subgraphs | Network fragmentation |
| `largest_component_size` | int | Genes in largest subgraph | Dominant connected component |
| `largest_component_fraction` | float | Largest component / n_genes | Connectivity of the main network |
| `assortativity` | float | Pearson r of degrees at edge endpoints | Positive = hubs connect to hubs |
| `mean_abs_corr` | float | Mean |r| of all significant edges | Average co-expression strength |
| `median_abs_corr` | float | Median |r| | Typical edge strength |
| `std_abs_corr` | float | Std |r| | Edge weight variability |
| `positive_edge_fraction` | float | Fraction of edges with r > 0 | Activation vs repression balance |
| `negative_edge_fraction` | float | Fraction of edges with r < 0 | Anti-correlation fraction |

### 8.2 Per-Gene Connectivity Metrics (Stage 3)

| Metric key | Type | Definition | Biological meaning |
|-----------|------|-----------|-------------------|
| `degree` | int | Number of significant co-expression partners | Direct connectivity |
| `weighted_degree` | float | Σ|r| over all edges | WGCNA-style soft connectivity (k) |
| `mean_abs_corr` | float | Mean |r| of edges | Average co-expression strength per partner |
| `max_abs_corr` | float | Max |r| of edges | Strongest co-expression relationship |
| `std_abs_corr` | float | Std |r| | Homogeneity of co-expression strength |
| `positive_edge_count` | int | Edges with r > 0 | Activator-like profile |
| `negative_edge_count` | int | Edges with r < 0 | Repressor-like profile |
| `positive_fraction` | float | positive_edge_count / degree | Positive correlation bias |

### 8.3 Per-Gene L1 Neighbourhood Metrics (Stage 3)

| Metric key | Type | Definition | Biological meaning |
|-----------|------|-----------|-------------------|
| `L1_n_nodes` | int | Number of direct partners (= degree) | Immediate influence radius |
| `L1_conn_sum` | float | Σ|r| of L1 edges | Total co-expression energy to direct partners |
| `L1_conn_mean` | float | Mean |r| of L1 edges | Average edge strength to direct partners |
| `L1_clique_density` | float | (L1-L1 edges) / max possible L1-L1 edges | How tightly the direct neighbourhood is internally connected; 1.0 = all partners are also co-expressed with each other |
| `n_l1_to_l1_edges` | int | Edges within L1 neighbourhood | Triangle count supporting clustering |
| `local_clustering_coefficient` | float | = L1_clique_density | Standard graph theory clustering (fraction of triangles through focus gene) |

### 8.4 Per-Gene L2 Neighbourhood Metrics (Stage 3)

| Metric key | Type | Definition | Biological meaning |
|-----------|------|-----------|-------------------|
| `L2_n_nodes` | int | Pure 2nd-hop nodes (partners of partners, not L1) | Network expansion at 2 hops |
| `L2_n_edges` | int | L1-L1 + L1-L2 edges | Connectivity of the outer ring |
| `L2_conn_sum` | float | Σ|r| of outer-ring edges | Co-expression energy in the outer ring |
| `L2_conn_mean` | float | Mean |r| of outer-ring edges | Average outer-ring edge strength |
| `n_l1_to_l2_edges` | int | Edges from L1 nodes to L2 nodes | Bridges to the extended neighbourhood |

### 8.5 Topological Ratio Metrics (Stage 3)

| Metric key | Formula | Biological meaning |
|-----------|---------|-------------------|
| `L2L1_deg` | L2_n_nodes / L1_n_nodes | **Expansion ratio**: < 1 → dense local cluster; > 1 → sparse hub with large outer ring |
| `L2L1_conn` | L2_conn_sum / L1_conn_sum | Relative connectivity of outer vs inner ring; > 1 → outer ring is strongly connected |
| `full_n_edges` | L1_n_nodes + L2_n_edges | Total edges in the 2-layer neighbourhood |
| `full_conn_sum` | L1_conn_sum + L2_conn_sum | Total co-expression energy in 2-layer neighbourhood |

### 8.6 Optional Advanced Centrality (Stage 4)

| Metric key | Method | Biological meaning | Complexity |
|-----------|--------|-------------------|------------|
| `eigenvector_centrality` | Power iteration on CSR adjacency matrix | Hub-of-hubs: high if connected to other high-degree genes | O(k × E), k≈100 iterations; ~minutes |
| `approx_betweenness` | Brandes sampling, n_samples pivot nodes | Bottleneck/bridge genes: how often a gene lies on shortest paths | O(n_samples × E); ~10-30 min with n_samples=500 |
| `pagerank` | Random walk steady state | Influence accounting for the full network topology | O(iterations × E) |

**Why eigenvector centrality instead of exact betweenness by default?**
For 18K genes and ~5M edges, exact betweenness requires O(V × E) = O(9×10¹⁰) operations
— infeasible in Python. Eigenvector centrality captures similar information (important
neighbours → high score) at O(100 × 5M) = O(5×10⁸), which runs in minutes.

---

## 9. Key Differences from Differential Pipeline

| Aspect | Differential pipeline | Network metrics pipeline |
|--------|----------------------|-------------------------|
| **Purpose** | What changes between LOW/HIGH conditions? | What topology structure exists? |
| **Input subsets** | LOW (bottom 20%) + HIGH (top 20%) | Full matrix, all samples |
| **FDR corrections** | 3× independent (low, high, Fisher-Z diff) | 1× global |
| **Fisher-Z test** | Yes (compares r_low vs r_high) | Not applicable |
| **Bootstrap** | Yes (confidence intervals on delta-r) | Not needed |
| **Adjacency count** | 3 (adj_low, adj_high, adj_diff) | 1 (adj_net) |
| **Edge classification** | disappear / new / sign_change / strengthen / weaken | Not applicable |
| **L2L1_rewire** | L2_rewire / L1_rewire (rewiring ratio) | **Not computed** (no comparison condition) |
| **L2L1_deg** | L2_n_nodes / L1_n_nodes | **Same formula** |
| **L2L1_conn** | L2_conn_diff / L1_conn_diff (uses Δr) | L2_conn_sum / L1_conn_sum (uses \|r\|) |
| **Betweenness** | Not computed | Approximate (optional, Stage 4) |
| **Eigenvector centrality** | Not computed | Default Stage 4 output |

---

## 10. Sanity Checks

After running, verify:

1. **`L1_n_nodes == degree`** for every gene (identity check — must always hold)
2. **`n_l1_to_l1_edges <= L1_n_nodes * (L1_n_nodes - 1) / 2`** (clique bound)
3. **`L1_clique_density ∈ [0, 1]`** for all genes
4. **`local_clustering_coefficient == L1_clique_density`** (same quantity, two names)
5. **`network_metrics.tsv` row count == n_genes_active** (all active genes present)
6. **`global_clustering_coefficient ∈ [0, 1]`**
7. **`density` consistent with `n_edges`**: `density ≈ 2*n_edges / (n_genes*(n_genes-1))`

---

## 11. Future Work

| ID | Description | Document |
|----|-------------|----------|
| FTD-NM-01 | Multi-dataset support: run the full pipeline over multiple expression matrices in a single Snakemake invocation, defined via a `datasets:` dict in the config | [PLAN-NM-01-Multi-Dataset-Support.md](PLAN-NM-01-Multi-Dataset-Support.md) |

---

## 12. Files to Delete After Refactor

The directory `src/scripts/20network_metrics/` was created in error and should be deleted:

```bash
git checkout HEAD -- src/scripts/  # restore if tracked
# or
rm -rf src/scripts/20network_metrics/
```

---

*Document version: 2026-04-08*
*Pipeline: `network_metrics/`*
*Concept reference: `docs/CONCEPT-01-Network-Analysis-Hierarchy.md`*
