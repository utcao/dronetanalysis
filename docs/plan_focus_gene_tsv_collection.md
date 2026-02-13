# Plan: Collect focus_gene analysis into a TSV across all per-gene h5 files

## Context

In `03_reconstruct_diff_network.py`, the **per-gene directory mode** (`collect_per_gene_networks`) calls `analyze_focus_gene()` for each reference gene and collects a comprehensive TSV with edge-level summaries AND focus_gene neighborhood metrics (degree, rewiring, strength at L1 and L1+L2 layers, plus ratios).

## File modified

**[03_reconstruct_diff_network.py](src/scripts/10spearman_corr/03_reconstruct_diff_network.py)**

No changes to `02b_bootstrap_significant_edges.py` — all metrics are derivable from existing stage 2a+2b outputs.

## Metric Definitions

**"Present"**: `present_low(edge) = sig_low AND |r_low| >= threshold` (same for high). This is consistent with `classify_qualitative_change()`.

**Strength**: `str_low = sum(|r_low|)` for present_low edges only; `str_diff = sum(|Δr|)` for all diff edges.

**Rewiring**: `n_disappear + n_new + n_sign_change` (strengthen/weaken NOT included in rewiring).

**Strengthen/Weaken**: Computed **independently** of sign_change. Among all "both present" edges (including sign_change), edges where `|r_high| > |r_low|` count as strengthen, `|r_high| < |r_low|` as weaken. This means a sign_change edge is also counted as strengthen or weaken.

**L1** = edges directly touching focus gene in diff network.
**L2** = full two-layer edge set: focus→L1 + L1↔L1 + L1→L2 edges (no duplicate edges).

## TSV Column Design

One row per reference gene with `n_sig_total > 0`, sorted by `L2L1_deg` descending.

### Basic info
| Column | Meaning |
|---|---|
| `gene_idx` | Reference gene index |
| `gene_id` | Reference gene ID |
| `n_sig_total` | Total significant differential edges in this gene's whole network |

### L1 degree
| Column | Meaning |
|---|---|
| `L1_deg_diff` | Focus gene partners in diff network (= sum of all non-overlapping qual categories from qual_score) |

### L1 qualitative counts
| Column | Meaning |
|---|---|
| `L1_n_disappear` | present_low AND NOT present_high at L1 |
| `L1_n_new` | NOT present_low AND present_high at L1 |
| `L1_n_sign_chg` | Both present, sign flipped at L1 |
| `L1_n_strengthen` | Both present, \|r_high\| > \|r_low\| at L1 (independent of sign_change) |
| `L1_n_weaken` | Both present, \|r_high\| < \|r_low\| at L1 (independent of sign_change) |
| `L1_rewire` | n_disappear + n_new + n_sign_chg at L1 |

### L1 strength
| Column | Meaning |
|---|---|
| `L1_str_low` | sum\|r_low\| for present_low edges at L1 |
| `L1_str_high` | sum\|r_high\| for present_high edges at L1 |
| `L1_str_diff` | sum\|Δr\| at L1 |
| `L1_mean_abs_dr` | mean\|Δr\| at L1 |

### L2 metrics (full two-layer neighborhood in diff network)
| Column | Meaning |
|---|---|
| `L2_deg_diff` | L1+L2 neighborhood size in diff network |
| `L2_rewire` | Rewiring score over full L2 edge set |
| `L2_str_low` | sum\|r_low\| for present_low edges at L2 |
| `L2_str_high` | sum\|r_high\| for present_high edges at L2 |
| `L2_str_diff` | sum\|Δr\| over full L2 edge set |

### L2/L1 expansion ratios (diff network)
| Column | Meaning |
|---|---|
| `L2L1_deg` | L2_deg_diff / L1_deg_diff |
| `L2L1_rewire` | L2_rewire / L1_rewire |
| `L2L1_str` | L2_str_diff / L1_str_diff |

### High/Low condition ratios
| Column | Meaning |
|---|---|
| `HL_str_L1` | L1_str_high / L1_str_low |
| `HL_str_L2` | L2_str_high / L2_str_low |

## Bug Fixes Applied

### Bug 1: adj_low/adj_high inconsistency
**Before**: `adj_low` used only `|r_low| >= threshold`, while `classify_qualitative_change` used `sig_low AND |r_low| >= threshold`.
**Fix**: `analyze_focus_gene()` now takes `sig_low`/`sig_high` parameters and builds `adj_low`/`adj_high` using the same "present" definition: `sig AND |r| >= threshold`.

### Bug 2: strengthen/weaken excluded sign_change edges
**Before**: `classify_qualitative_change` assigns `QUAL_SIGN_CHANGE` exclusively, so `n_strengthen`/`n_weaken` from `qual_score` would be 0 for sign_change edges (5 sign changes → 0 strengthen, 0 weaken).
**Fix**: `compute_edge_stats()` now counts strengthen/weaken independently by comparing `|r_high|` vs `|r_low|` for ALL "both present" edges (including sign_change). Note: `qual_score` itself is unchanged — the independence is only in the aggregate counts.

### Cleanup: Dropped redundant columns
- Removed `L1_deg_low`, `L1_deg_high`, `L2_deg_low`, `L2_deg_high` (derivable from qualitative counts).
- Renamed "connectivity" → "strength" (proper network science term).
- Renamed `ratio_*` → `L2L1_*` and added `HL_str_*` for High/Low condition ratios.

## Verification

```bash
python src/scripts/10spearman_corr/03_reconstruct_diff_network.py \
    --base-dir results_test_hs/base_correlations \
    --boot-dir results_test_hs/bootstrap_significant \
    --out-h5 /tmp/test_summary.h5 \
    --out-focus-tsv /tmp/test_focus.tsv
```
Checks (all verified):
- `test_focus.tsv` has one header + 95 data rows (one per active gene)
- Rows sorted by `L2L1_deg` descending
- `L1_rewire` = `L1_n_disappear + L1_n_new + L1_n_sign_chg` for all rows
- Where `L1_n_sign_chg > 0`, `L1_n_strengthen + L1_n_weaken > 0` (bug fix confirmed)
- `HL_str_L1` = `L1_str_high / L1_str_low` (within tolerance)
- `L2L1_deg` = `L2_deg_diff / L1_deg_diff` (within tolerance)
- No duplicate edges in L2 computation
- `test_summary.h5` has `per_gene/` group with all arrays
