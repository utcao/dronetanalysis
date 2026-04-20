# H5 Output File Structure Reference

## Overview

This document describes the structure and contents of the two HDF5 output types produced by the dronetanalysis pipeline: the genome-wide **differential network summary** file and the per-gene **network detail files**. Understanding these structures is essential for downstream analysis and result interpretation.

---

## File Types

### 1. Differential Network Summary (`differential_network_summary_*.h5`)

A single genome-wide file. For every gene in the expression matrix it stores ~60 pre-computed metrics derived from the differential co-expression network. Designed for fast cross-gene ranking and comparison without loading individual network files.

**Top-level structure:**

| Path | Shape | Content |
|------|-------|---------|
| `gene_names` | (n_genes,) | FBgn gene IDs (bytes) |
| `per_gene/<metric>` | (n_genes,) | One array per metric (see tables below) |
| `meta/` | — | Run parameters stored as HDF5 group attributes |

**`per_gene/` metrics — Focus gene degree:**

| Field | Type | Description |
|-------|------|-------------|
| `focus_deg_high` | int32 | # direct edges in HIGH condition |
| `focus_deg_low` | int32 | # direct edges in LOW condition |

**`per_gene/` metrics — L1 (direct / 1-hop neighborhood):**

| Field | Type | Description |
|-------|------|-------------|
| `L1_n_edges_high` | int32 | Edge count in HIGH network |
| `L1_n_edges_low` | int32 | Edge count in LOW network |
| `L1_n_nodes` | int32 | Unique L1 neighbor count |
| `L1_n_sign_chg` | int32 | Edges that flipped sign |
| `L1_n_disappear` | int32 | Edges lost in HIGH vs LOW |
| `L1_n_new` | int32 | Edges gained in HIGH vs LOW |
| `L1_n_strengthen` | int32 | Edges that strengthened |
| `L1_n_weaken` | int32 | Edges that weakened |
| `L1_rewire` | int32 | Total rewired edges |
| `L1_frac_rewire` | float32 | Fraction of edges rewired |
| `L1_frac_disappear` | float32 | Fraction of edges lost |
| `L1_frac_new` | float32 | Fraction of edges gained |
| `L1_conn_high` | float32 | Sum of absolute correlations (HIGH) |
| `L1_conn_low` | float32 | Sum of absolute correlations (LOW) |
| `L1_conn_diff` | float32 | Absolute difference in connectivity |
| `L1_conn_mean_high` | float32 | Mean absolute correlation (HIGH) |
| `L1_conn_mean_low` | float32 | Mean absolute correlation (LOW) |
| `L1_mean_abs_dr` | float32 | Mean \|Δr\| across L1 edges |
| `L1_mean_delta` | float32 | Mean signed Δr across L1 edges |
| `L1_clique_density` | float32 | Clique density within L1 neighborhood |

**`per_gene/` metrics — L2 (2-hop neighborhood):**

Mirrors the L1 set with `L2_` prefix, plus:

| Field | Type | Description |
|-------|------|-------------|
| `L2_n_edges` | int32 | Total edges in L2 subnetwork |
| `L2_n_nodes` | int32 | Unique L2 neighbor count |

**`per_gene/` metrics — Ratios and global:**

| Field | Type | Description |
|-------|------|-------------|
| `HL_conn_L1` | float32 | HIGH/LOW connectivity ratio at L1 |
| `HL_conn_L2` | float32 | HIGH/LOW connectivity ratio at L2 |
| `L2L1_conn` | float32 | L2/L1 connectivity ratio |
| `L2L1_deg` | float32 | L2/L1 degree ratio |
| `L2L1_rewire` | float32 | L2/L1 rewiring ratio |
| `max_abs_delta` | float32 | Maximum \|Δr\| in this gene's network |
| `mean_abs_delta` | float32 | Mean \|Δr\| across all network edges |
| `n_sig_edges_diff` | int32 | Total significant differential edges |
| `n_sign_change` | int32 | Sign-change edges (network-wide) |
| `n_disappear` | int32 | Disappearing edges (network-wide) |
| `n_new` | int32 | New edges (network-wide) |
| `n_strengthen` | int32 | Strengthened edges (network-wide) |
| `n_weaken` | int32 | Weakened edges (network-wide) |
| `n_l1_to_l1_edges` | int32 | Edges between L1 neighbors |
| `n_l1_to_l2_edges` | int32 | Edges from L1 to L2 nodes |
| `full_*` | float32/int32 | Full-network connectivity metrics (when enabled) |
| `gene_index` | int32 | Index of this gene in `gene_names` |
| `gene_name` | object | FBgn ID of this gene |

---

### 2. Per-gene Network Files (`networks/<index>_<FBgn>.h5`)

One file per focus gene. Named `<gene_index>_<FBgn>.h5` (e.g. `7246_FBgn0030331.h5`). Contains the same metrics as the summary file **plus** full edge-level detail, partner gene identities, a sparse delta-r matrix, and global topology.

**Top-level groups:**

| Group | Content |
|-------|---------|
| `edges/` | All significant differential edges as flat arrays |
| `focus_gene/` | Focus-gene-centric view: partner lists + metric attrs |
| `matrices/` | Sparse CSR delta-r matrix for all gene pairs |
| `topology/` | Degree distributions for high / low / diff networks |
| `gene_names` | (n_genes,) — shared gene ID lookup |
| `meta/` | Run parameters as HDF5 group attributes |

**`edges/` datasets** — one entry per significant edge (e.g. 68,260 edges):

| Dataset | Type | Description |
|---------|------|-------------|
| `gene_i` | int32 | Index of first gene in pair |
| `gene_j` | int32 | Index of second gene in pair |
| `r_high` | float32 | Pearson r in HIGH group |
| `r_low` | float32 | Pearson r in LOW group |
| `delta_base` | float32 | Base Δr (r_high − r_low) |
| `delta_boot_mean` | float32 | Bootstrap mean Δr |
| `delta_boot_std` | float32 | Bootstrap std of Δr |
| `bias` | float32 | Bootstrap bias |
| `ci_low` | float32 | Lower CI bound (alpha=0.05) |
| `ci_high` | float32 | Upper CI bound |
| `pval_boot` | float32 | Bootstrap p-value |
| `pval_diff` | float32 | Differential p-value |
| `qval_diff` | float32 | FDR-corrected q-value |
| `qual_label` | object | Quality label: `sign_change`, `strengthen`, `weaken`, `new`, `disappear` |
| `qual_score` | int8 | Integer quality score |

**`edges/qual_summary` attributes** — aggregate counts across all edges:

| Attribute | Example value | Description |
|-----------|--------------|-------------|
| `n_sign_change` | 53049 | Edges with correlation sign flip |
| `n_disappear` | 7591 | Edges absent in HIGH |
| `n_new` | 7230 | Edges absent in LOW |
| `n_strengthen` | 205 | Edges with larger \|r\| in HIGH |
| `n_weaken` | 185 | Edges with smaller \|r\| in HIGH |
| `n_unchanged` | 0 | Edges with no significant change |

**`focus_gene/` contents:**

| Path | Type | Description |
|------|------|-------------|
| `direct_partners` | int32 (n_L1,) | Gene indices of L1 direct neighbors |
| `indirect_partners` | int32 (n_L2,) | Gene indices of L2 indirect neighbors |
| `direct_stats` (attrs) | — | L1 edge summary: `n_edges`, `n_sign_change`, `conn_sum_high/low`, `mean_abs_delta`, etc. |
| `two_layer_stats` (attrs) | — | L2 edge summary (same fields) |
| `metrics` (attrs) | — | Full metric set (identical to `per_gene/` row in summary file) |

**`matrices/` datasets** — sparse CSR format of all Δr values:

| Dataset | Shape | Description |
|---------|-------|-------------|
| `delta_data` | (2 × n_edges,) | Non-zero Δr values (symmetric) |
| `delta_indices` | (2 × n_edges,) | Column indices |
| `delta_indptr` | (n_genes + 1,) | Row pointers |

**`topology/` datasets:**

| Path | Shape | Description |
|------|-------|-------------|
| `global_high/degrees` | (n_genes,) | Degree per gene in HIGH network |
| `global_low/degrees` | (n_genes,) | Degree per gene in LOW network |
| `global_diff/degrees` | (n_genes,) | Degree per gene in differential network |

**`meta/` attributes** (example values from `run_voomct`):

| Attribute | Value | Description |
|-----------|-------|-------------|
| `k_high` | 187 | Sample size of HIGH group |
| `k_low` | 187 | Sample size of LOW group |
| `n_genes` | 8763 | Total genes in expression matrix |
| `n_significant` | 68260 | Significant differential edges retained |
| `n_tests` | 38,390,703 | Total gene-pair tests performed |
| `fdr_alpha` | 0.7 | FDR threshold for edge selection |
| `ci_alpha` | 0.05 | CI confidence level |
| `corr_threshold` | 1e-4 | Minimum correlation threshold |
| `edge_selection` | `sig_differential` | Edge selection strategy |
| `require_ci_exclude_zero` | True | Whether CI must exclude zero |
| `min_effect` | 0.0 | Minimum effect size filter |

---

## Summary vs. Network Files: What Each Adds

| Information | Summary `.h5` | Network `.h5` |
|-------------|:---:|:---:|
| Per-gene metrics (L1/L2/ratios) | ✅ all 8763 genes | ✅ focus gene only (in attrs) |
| Individual edge statistics (r, Δr, p, q, CI) | ✗ | ✅ |
| Edge quality labels per edge | ✗ | ✅ |
| Partner gene identities (FBgn IDs) | ✗ | ✅ via `direct_partners` + `gene_names` |
| Sparse full delta-r matrix | ✗ | ✅ |
| Global degree distributions | ✗ | ✅ |
| Run parameters | ✗ | ✅ in `meta/` attrs |
| Edge quality aggregate counts | ✗ | ✅ in `edges/qual_summary` attrs |

**Use the summary file** for genome-wide ranking, filtering, and cross-gene comparisons.  
**Use the network files** when you need to inspect specific edges, retrieve partner gene lists, or reconstruct the network for a focus gene.

---

**Last Updated:** 2026-04-17  
**Status:** ✅ Active reference  
**Related Reading:**
- [GUIDE-04-Qualitative-Change-Metrics.md](GUIDE-04-Qualitative-Change-Metrics.md) — Quality labels (sign_change / disappear / new / strengthen / weaken) explained
- [GUIDE-02-Network-Metrics.md](GUIDE-02-Network-Metrics.md) — L1/L2 network topology metrics interpretation
- [REFERENCE-01-Statistical-Methods.md](REFERENCE-01-Statistical-Methods.md) — Bootstrap delta-r testing and FDR methodology
- [OPTIMIZATION-02-Storage.md](OPTIMIZATION-02-Storage.md) — Storage modes that affect what is saved in network files
