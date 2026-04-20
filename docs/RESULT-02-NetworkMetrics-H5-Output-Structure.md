# Network Metrics Pipeline — HDF5 Output File Structures

## Overview

This document describes the structure and contents of every HDF5 (and TSV) file produced by the `network_metrics` sub-pipeline. The pipeline operates on a **single expression matrix** (all samples, no condition split) and produces a co-expression network with per-gene topology metrics.

**Pipeline stages and their outputs:**

```
Stage 1  01_compute_correlations.py   →  corr_significant.h5
Stage 2  02_build_network.py          →  network_edges.h5
Stage 3  03_per_gene_metrics.py       →  per_gene/{gi}_{gene_id}.h5  (×n_genes)
Stage 3b 03b_collect_metrics.py       →  network_metrics_summary.h5
                                          network_metrics.tsv
Stage 4  04_global_metrics.py         →  global_topology.h5          (optional)
```

> **Relationship to the differential pipeline:** These outputs describe a single co-expression network (all samples). They do not contain HIGH/LOW condition splits, rewiring metrics, delta-r, or qualitative change categories. For those, see [RESULT-01-H5-Output-File-Structure.md](RESULT-01-H5-Output-File-Structure.md).

---

## Stage 1 Output: `corr_significant.h5`

All-vs-all Spearman correlations on the full expression matrix, filtered to significant edges after a single global BH-FDR correction.

**Top-level structure:**

| Path | Content |
|------|---------|
| `meta/` | Run parameters as HDF5 group attributes |
| `edges/` | Significant edge arrays |
| `gene_names` | (n_genes,) str dataset |
| `sample_names` | (n_samples,) str dataset |

**`meta/` attributes:**

| Attribute | Type | Description |
|-----------|------|-------------|
| `n_genes` | int | Total genes in expression matrix |
| `n_samples` | int | Total samples |
| `n_tests` | int | Total gene-pair tests = n_genes×(n_genes−1)/2 |
| `n_significant` | int | Edges surviving FDR + corr threshold filters |
| `fdr_alpha` | float | BH-FDR significance threshold (default: 0.05) |
| `corr_threshold` | float | Minimum \|r\| retained (default: 1e-4) |

**`edges/` datasets** — one entry per significant edge (upper-triangle only, gene_i < gene_j):

| Dataset | Type | Description |
|---------|------|-------------|
| `gene_i` | int32 (n_sig,) | Row gene index (always < gene_j) |
| `gene_j` | int32 (n_sig,) | Column gene index (always > gene_i) |
| `corr` | float32 (n_sig,) | Spearman r |
| `qval` | float32 (n_sig,) | BH-FDR q-value |

---

## Stage 2 Output: `network_edges.h5`

Builds on Stage 1 by adding a CSR-format undirected adjacency index and global network topology metrics. This file is the primary input for all Stage 3 per-gene jobs.

**Top-level structure:**

| Path | Content |
|------|---------|
| `meta/` | Run parameters as HDF5 group attributes |
| `edges/` | Edge arrays (copied from Stage 1) |
| `adjacency/` | CSR adjacency for O(degree) per-gene reads |
| `global_topology/` | Network-wide topology metrics |
| `gene_names` | (n_genes,) str dataset |

**`meta/` attributes:**

| Attribute | Type | Description |
|-----------|------|-------------|
| `n_genes` | int | Total genes |
| `n_edges` | int | Significant edges |
| `n_genes_active` | int | Genes with at least one edge |
| `density` | float | Network density = 2E / (N×(N−1)) |
| `fdr_alpha` | float | BH-FDR threshold inherited from Stage 1 |
| `corr_threshold` | float | Minimum \|r\| inherited from Stage 1 |

**`edges/` datasets** — same as Stage 1:

| Dataset | Type | Description |
|---------|------|-------------|
| `gene_i` | int32 (n_edges,) | Gene pair index (upper triangle) |
| `gene_j` | int32 (n_edges,) | Gene pair index |
| `corr` | float32 (n_edges,) | Spearman r |
| `qval` | float32 (n_edges,) | BH-FDR q-value |

**`adjacency/` datasets** — CSR-format undirected adjacency (both directions stored):

| Dataset | Type | Description |
|---------|------|-------------|
| `indptr` | int64 (n_genes+1,) | Row pointers; `indptr[g]:indptr[g+1]` gives the slice for gene g |
| `indices` | int32 (2×n_edges,) | Neighbour indices, sorted within each row |
| `data` | float32 (2×n_edges,) | Spearman r values aligned to `indices` |

> **Usage:** To retrieve all neighbours of gene `g`: `nb = indices[indptr[g]:indptr[g+1]]`, `r = data[indptr[g]:indptr[g+1]]`. This is the O(degree) read pattern used by Stage 3.

**`global_topology/` attributes** — all stored as HDF5 group attributes except `degrees`:

| Attribute | Type | Description |
|-----------|------|-------------|
| `n_genes` | int | Total gene count |
| `n_genes_active` | int | Genes with degree > 0 |
| `fraction_active` | float | n_genes_active / n_genes |
| `n_edges` | int | Total edge count |
| `density` | float | Network density |
| `avg_degree` | float | Mean degree (over all genes) |
| `median_degree` | float | Median degree (active genes only) |
| `max_degree` | int | Maximum degree |
| `min_degree` | int | Minimum degree (active genes only) |
| `std_degree` | float | Std of degree (active genes only) |
| `degree_p25/50/75/90/95/99` | float | Degree percentiles (active genes) |
| `power_law_exponent` | float | γ from log-log power-law fit; γ=2–3 → scale-free |
| `power_law_r_squared` | float | R² of power-law fit |
| `is_scale_free` | int (bool) | 1 if γ ∈ [2.0, 3.5] and R² > 0.8 |
| `global_clustering` | float | Mean local clustering coefficient (nodes with degree ≥ 2) |
| `n_components` | int | Number of connected components |
| `largest_component_size` | int | Node count in largest component |
| `largest_component_fraction` | float | largest_component_size / n_genes |
| `assortativity` | float | Degree assortativity (Pearson r of edge endpoint degrees) |
| `mean_abs_corr` | float | Mean \|r\| across all edges |
| `median_abs_corr` | float | Median \|r\| across all edges |
| `std_abs_corr` | float | Std \|r\| across all edges |
| `positive_edge_fraction` | float | Fraction of edges with r > 0 |
| `negative_edge_fraction` | float | Fraction of edges with r < 0 |

**`global_topology/` dataset:**

| Dataset | Type | Description |
|---------|------|-------------|
| `degrees` | int32 (n_genes,) | Per-gene degree array |

---

## Stage 3 Output: `per_gene/{gi}_{gene_id}.h5`

One file per focus gene, named `{gene_index:04d}_{gene_id}.h5` (e.g. `0042_FBgn0001234.h5`). Contains L1/L2 neighbourhood metrics and partner gene indices.

**Top-level structure:**

| Path | Content |
|------|---------|
| `meta/` | Gene identity and network context as HDF5 group attributes |
| `metrics/` | All per-gene metrics as HDF5 group attributes |
| `partners/` | Direct and indirect partner gene indices |

**`meta/` attributes:**

| Attribute | Type | Description |
|-----------|------|-------------|
| `gene_index` | int | Zero-based gene index |
| `gene_name` | str | Gene identifier (e.g. FBgn ID) |
| `n_genes` | int | Total genes in network |
| `n_edges` | int | Total edges in network |

**`metrics/` attributes** — all stored as HDF5 group attributes:

| Attribute | Type | Description |
|-----------|------|-------------|
| `degree` | int32 | # direct edges (= L1 neighbour count) |
| `weighted_degree` | float32 | Sum of \|r\| over direct edges |
| `mean_abs_corr` | float32 | Mean \|r\| over direct edges |
| `max_abs_corr` | float32 | Max \|r\| over direct edges |
| `std_abs_corr` | float32 | Std of \|r\| over direct edges |
| `positive_edge_count` | int32 | # direct edges with r > 0 |
| `negative_edge_count` | int32 | # direct edges with r < 0 |
| `positive_fraction` | float32 | positive_edge_count / degree |
| `L1_n_nodes` | int32 | # L1 (direct) neighbours (= degree) |
| `L1_conn_sum` | float32 | Sum \|r\| over L1 edges (= weighted_degree) |
| `L1_conn_mean` | float32 | Mean \|r\| over L1 edges (= mean_abs_corr) |
| `L1_clique_density` | float32 | Fraction of L1 pairs that are connected to each other |
| `n_l1_to_l1_edges` | int32 | # edges among L1 neighbours (triangles through focus) |
| `local_clustering_coefficient` | float32 | Same as L1_clique_density |
| `L2_n_nodes` | int32 | # L2 (2-hop) neighbours (excludes focus and L1) |
| `L2_n_edges` | int32 | # edges in L2 ring = n_l1_to_l1_edges + n_l1_to_l2_edges |
| `L2_conn_sum` | float32 | Sum \|r\| over L2-ring edges |
| `L2_conn_mean` | float32 | Mean \|r\| over L2-ring edges |
| `n_l1_to_l2_edges` | int32 | # edges from L1 nodes to L2 nodes |
| `L2L1_deg` | float32 | L2_n_nodes / L1_n_nodes — neighbourhood expansion ratio |
| `L2L1_conn` | float32 | L2_conn_sum / L1_conn_sum — connectivity expansion ratio |
| `full_n_edges` | int32 | L1_n_nodes + L2_n_edges (total edges in 2-layer subgraph) |
| `full_conn_sum` | float32 | L1_conn_sum + L2_conn_sum |

**`partners/` datasets:**

| Dataset | Type | Description |
|---------|------|-------------|
| `direct` | int32 (L1_n_nodes,) | Gene indices of L1 direct neighbours (sorted) |
| `indirect` | int32 (L2_n_nodes,) | Gene indices of L2 indirect neighbours (sorted) |

> **Note:** `gene_names` is not stored in the per-gene file. Use the gene index to look up the name from `network_edges.h5 / gene_names`.

---

## Stage 3b Output: `network_metrics_summary.h5` + `network_metrics.tsv`

Aggregates all per-gene files from Stage 3 into a single genome-wide summary. Mirrors the structure used by the differential pipeline's `differential_network_summary_*.h5`.

### `network_metrics_summary.h5`

**Top-level structure:**

| Path | Content |
|------|---------|
| `meta/` | Counts as HDF5 group attributes |
| `gene_names` | (n_genes,) str — full gene list from `network_edges.h5` |
| `per_gene/` | One array per metric across all reference genes |

**`meta/` attributes:**

| Attribute | Type | Description |
|-----------|------|-------------|
| `n_genes` | int | Total genes in expression matrix |
| `n_ref_genes` | int | Genes for which a per_gene .h5 file was found |
| `n_edges` | int | Total edges in the network |

**`per_gene/` datasets** — shape (n_ref_genes,) each:

| Dataset | Type | Description |
|---------|------|-------------|
| `gene_index` | int32 | Gene index in global `gene_names` array |
| `gene_name` | str | Gene identifier parsed from filename |
| `degree` | int32 | Direct neighbour count |
| `weighted_degree` | float32 | Sum \|r\| (L1) |
| `mean_abs_corr` | float32 | Mean \|r\| (L1) |
| `max_abs_corr` | float32 | Max \|r\| (L1) |
| `std_abs_corr` | float32 | Std \|r\| (L1) |
| `positive_edge_count` | int32 | # positive L1 edges |
| `negative_edge_count` | int32 | # negative L1 edges |
| `positive_fraction` | float32 | Fraction of positive L1 edges |
| `L1_n_nodes` | int32 | L1 neighbour count |
| `L1_conn_sum` | float32 | Sum \|r\| over L1 |
| `L1_conn_mean` | float32 | Mean \|r\| over L1 |
| `L1_clique_density` | float32 | Fraction of L1 pairs connected |
| `n_l1_to_l1_edges` | int32 | L1–L1 edge count |
| `local_clustering_coefficient` | float32 | = L1_clique_density |
| `L2_n_nodes` | int32 | L2 neighbour count |
| `L2_n_edges` | int32 | Edges in L2 ring |
| `L2_conn_sum` | float32 | Sum \|r\| over L2 ring |
| `L2_conn_mean` | float32 | Mean \|r\| over L2 ring |
| `n_l1_to_l2_edges` | int32 | L1→L2 edge count |
| `L2L1_deg` | float32 | L2/L1 neighbourhood expansion |
| `L2L1_conn` | float32 | L2/L1 connectivity expansion |
| `full_n_edges` | int32 | Total edges in 2-layer subgraph |
| `full_conn_sum` | float32 | Total \|r\| sum in 2-layer subgraph |

### `network_metrics.tsv`

Flat tab-separated version of `network_metrics_summary.h5 / per_gene/`, sorted by `degree` descending.

**Columns:** `gene_idx`, `gene_id`, then all metric keys in order:
`degree`, `weighted_degree`, `mean_abs_corr`, `max_abs_corr`, `std_abs_corr`, `positive_edge_count`, `negative_edge_count`, `positive_fraction`, `L1_n_nodes`, `L1_conn_sum`, `L1_conn_mean`, `L1_clique_density`, `n_l1_to_l1_edges`, `local_clustering_coefficient`, `L2_n_nodes`, `L2_n_edges`, `L2_conn_sum`, `L2_conn_mean`, `n_l1_to_l2_edges`, `L2L1_deg`, `L2L1_conn`, `full_n_edges`, `full_conn_sum`

Integer columns written as integers; float columns written with 6 significant figures.

---

## Stage 4 Output: `global_topology.h5` (optional)

Eigenvector centrality and PageRank for every gene. Approximate betweenness is computed only when `--compute-betweenness` is passed (requires `networkit` package; ~10–30 min for 18K genes / 5M edges).

**Top-level structure:**

| Path | Content |
|------|---------|
| `meta/` | Counts and computation flags as HDF5 group attributes |
| `per_gene/` | Centrality score arrays |

**`meta/` attributes:**

| Attribute | Type | Description |
|-----------|------|-------------|
| `n_genes` | int | Total genes |
| `n_edges` | int | Total edges |
| `computed_betweenness` | int (bool) | 1 if approximate betweenness was computed |
| `betweenness_n_samples` | int | Pivot nodes used for Brandes sampling (only if betweenness computed) |

**`per_gene/` datasets** — shape (n_genes,) float32 each:

| Dataset | Always present | Description |
|---------|:-:|-------------|
| `eigenvector_centrality` | ✅ | Power-iteration eigenvector centrality using \|r\| as edge weights; L2-normalised |
| `pagerank` | ✅ | PageRank with damping=0.85, row-stochastic \|r\|-weighted transition matrix |
| `approx_betweenness` | Only with `--compute-betweenness` | Approximate Brandes betweenness (normalised); error O(1/√n_samples) |

> **Betweenness note:** Exact Brandes is O(V×E) ≈ 9×10¹⁰ ops for 18K genes / 5M edges and is not feasible. The approximate version using `nSamples=500` takes ~10–30 min. Use `nSamples=1000+` for publication-quality scores.

---

## Pipeline Output Map

```
network_metrics/
├── corr_significant.h5          ← Stage 1: Spearman r + BH-FDR, significant edges only
├── network_edges.h5             ← Stage 2: CSR adjacency + global topology
├── per_gene/
│   ├── 0000_FBgn*.h5            ← Stage 3: L1/L2 metrics + partner indices (per gene)
│   ├── 0001_FBgn*.h5
│   └── ...
├── network_metrics_summary.h5   ← Stage 3b: All per-gene metrics in one array table
├── network_metrics.tsv          ← Stage 3b: Same content as flat TSV (sorted by degree)
└── global_topology.h5           ← Stage 4 (optional): Eigenvector centrality + PageRank
```

---

## Key Differences vs. Differential Pipeline Outputs

| Feature | Network Metrics pipeline | Differential pipeline |
|---------|:---:|:---:|
| Condition split (HIGH / LOW) | ✗ | ✅ |
| Rewiring metrics (Δr, sign_change, disappear…) | ✗ | ✅ |
| Edge-level delta statistics (CI, p-values) | ✗ | ✅ |
| Single Spearman r per edge | ✅ | ✗ |
| `positive_fraction` metric | ✅ | ✗ |
| Eigenvector centrality / PageRank | ✅ (Stage 4) | ✗ |
| Approximate betweenness | ✅ (Stage 4, optional) | ✗ |
| Global topology (clustering, components, scale-free) | ✅ (`network_edges.h5`) | ✅ (`topology/` in per-gene files) |
| CSR adjacency index | ✅ (`network_edges.h5`) | ✅ (`matrices/` in per-gene files) |

---

**Last Updated:** 2026-04-17  
**Status:** ✅ Active reference  
**Related Reading:**
- [RESULT-01-H5-Output-File-Structure.md](RESULT-01-H5-Output-File-Structure.md) — Differential pipeline output structures (summary + per-gene network files with HIGH/LOW split)
- [GUIDE-02-Network-Metrics.md](GUIDE-02-Network-Metrics.md) — Interpretation of L1/L2 topology metrics
- [REFERENCE-01-Statistical-Methods.md](REFERENCE-01-Statistical-Methods.md) — Spearman correlation, BH-FDR, and bootstrap methodology
- [OPTIMIZATION-02-Storage.md](OPTIMIZATION-02-Storage.md) — Storage modes and compression settings
