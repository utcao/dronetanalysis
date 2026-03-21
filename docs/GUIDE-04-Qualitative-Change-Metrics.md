# Qualitative Change Metrics Guide

## Overview

Stage 3 classifies every significant differential edge into one of **6 mutually exclusive categories** based on whether the edge is "present" in the low and/or high condition network. These categories, together with the L1/L2 focus-gene neighbourhood metrics, are the primary outputs for interpreting *how* co-expression rewires between conditions.

This guide documents:
1. The 6 qualitative categories and their exact definitions
2. What "present" means and why it matters
3. The partition identity (sanity check)
4. L1/L2 focus-gene neighbourhood metrics
5. `L1_n_edges_low` / `L1_n_edges_high` and their relationship to `focus_deg_*`
6. How to sanity-check output tables

---

## 1. What "Present" Means

An edge is **present** in a condition when it passes **both** criteria:

```
present_low  = sig_low  AND |r_low|  >= corr_threshold
present_high = sig_high AND |r_high| >= corr_threshold
```

- `sig_low` / `sig_high`: individually significant in the low / high condition (from bootstrap)
- `corr_threshold` (default `0.0001`): minimum absolute correlation to be considered real

An edge can be in the **differential network** (significant delta between conditions) without being "present" in either individual condition — this produces UNCHANGED edges (see section 3).

---

## 2. The 6 Qualitative Categories

| Code | Name | Condition | Biological meaning |
|------|------|-----------|-------------------|
| 0 | **UNCHANGED** | `~present_low AND ~present_high` | Delta is significant but both correlations are too weak to be called "present" in either condition. Rare with `sig_edges` selection; possible with `sig_differential`. |
| 1 | **DISAPPEAR** | `present_low AND ~present_high` | Correlation exists in low expression, lost in high expression. |
| 2 | **NEW** | `~present_low AND present_high` | Correlation absent in low, gains significance in high expression. |
| 3 | **SIGN_CHANGE** | `present_low AND present_high AND sign(r_low) ≠ sign(r_high)` | Correlation flips direction (positive ↔ negative) between conditions. Strong rewiring signal. |
| 4 | **STRENGTHEN** | `present_low AND present_high AND same sign AND \|r_high\| > \|r_low\|` | Correlation intensifies in high expression, same direction. |
| 5 | **WEAKEN** | `present_low AND present_high AND same sign AND \|r_high\| < \|r_low\|` | Correlation weakens in high expression, same direction. |

**Note on same-sign ties** (`|r_high| == |r_low|`): these edges remain UNCHANGED (code 0). Practically impossible with continuous float values.

### Conceptual grouping

```
                 present_low   present_high
DISAPPEAR            ✓             ✗        ← unidirectional rewiring
NEW                  ✗             ✓        ← unidirectional rewiring
SIGN_CHANGE          ✓             ✓        ← bidirectional rewiring (sign flip)
STRENGTHEN           ✓             ✓        ← stable, gaining strength
WEAKEN               ✓             ✓        ← stable, losing strength
UNCHANGED            ✗             ✗        ← in diff network, absent from both
```

**Rewiring edges**: DISAPPEAR + NEW + SIGN_CHANGE (structural change in network topology)
**Magnitude-change edges**: STRENGTHEN + WEAKEN (topology preserved, strength changes)

---

## 3. Partition Identity (Global Sanity Check)

The 6 categories partition **all** significant differential edges exactly:

```
n_unchanged + n_disappear + n_new + n_sign_change + n_strengthen + n_weaken = n_sig_edges_diff
```

In the per-gene summary output (`rewiring_hubs.tsv`, `differential_network_summary.h5`), only the 5 non-UNCHANGED categories are stored. If UNCHANGED edges exist:

```
n_disappear + n_new + n_sign_change + n_strengthen + n_weaken < n_sig_edges_diff
```

The gap equals `n_unchanged`. With `edge_selection = sig_edges` (default = `sig_differential`), UNCHANGED edges are rare or absent. With `sig_differential`, they can exist for edges with significant delta but weak correlations in both individual conditions.

---

## 4. STRENGTHEN and WEAKEN: Pure Definition

**STRENGTHEN and WEAKEN are pure qual_score codes** (codes 4 and 5). They are assigned only when:
- edge is present in **both** conditions
- signs are the **same**
- the magnitudes differ

SIGN_CHANGE edges are **not** included in STRENGTHEN or WEAKEN counts, even if `|r_high| > |r_low|` for a sign-changing edge. This ensures the clean partition holds and avoids double-counting.

This applies consistently everywhere: global `qual_summary`, per-gene summary stats in `03b`, and L1/L2 focus gene stats.

---

## 5. L1 Focus-Gene Neighbourhood Metrics

### Layer definitions

For a given focus gene, all metrics are computed on the **differential network** (significant delta edges only):

| Layer | Contents |
|-------|---------|
| **L1 (direct)** | Edges between focus gene and its direct differential partners |
| **L2 outer** | Edges among L1 partners (L1↔L1) and from L1 partners to 2nd-hop nodes (L1→L2) |
| **Full 2-layer** | L1 + L2 outer combined |

### Key identity: `L1_n_nodes == n_direct_edges`

Each direct partner contributes exactly **one** edge to the focus gene. Therefore:

```
L1_n_nodes = number of direct differential partners = number of L1 edges
```

### L1 partition identity

```
L1_n_disappear + L1_n_new + L1_n_sign_chg + L1_n_strengthen + L1_n_weaken + n_unchanged_L1
  = L1_n_nodes
```

When UNCHANGED ≈ 0 (typical):

```
L1_n_disappear + L1_n_new + L1_n_sign_chg + L1_n_strengthen + L1_n_weaken ≈ L1_n_nodes
```

### L1 rewiring score

```
L1_rewire = L1_n_disappear + L1_n_new + L1_n_sign_chg
```

This counts structurally rewired edges only (existence or sign changed). STRENGTHEN and WEAKEN are stable edges that do **not** contribute to rewiring.

---

## 6. Condition-Specific Degree: `focus_deg_low/high` and `L1_n_edges_low/high`

### What these metrics measure

`focus_deg_low` and `L1_n_edges_low` measure the **same quantity**: among the L1 differential edges of the focus gene, how many are also "present" in the low condition network.

```
L1_n_edges_low  = focus_deg_low  = |{ L1 edges : present_low  = True }|
L1_n_edges_high = focus_deg_high = |{ L1 edges : present_high = True }|
```

`focus_deg_low` is computed by `len(adj_low[focus_gene])` where `adj_low` is built only from the differential edge set — so it is **not** the full degree in the low network, only among differential partners.

### Derivation from categories

```
L1_n_edges_low  = L1_n_disappear + L1_n_sign_chg + L1_n_strengthen + L1_n_weaken
                = L1_n_nodes - L1_n_new - n_unchanged_L1

L1_n_edges_high = L1_n_new      + L1_n_sign_chg + L1_n_strengthen + L1_n_weaken
                = L1_n_nodes - L1_n_disappear - n_unchanged_L1
```

Equivalently:
```
L1_n_edges_low  = L1_n_nodes - L1_n_new      (approximately, when UNCHANGED ≈ 0)
L1_n_edges_high = L1_n_nodes - L1_n_disappear (approximately, when UNCHANGED ≈ 0)
```

**Interpretation:**
- A high `focus_deg_low` relative to `L1_n_nodes` means most differential partners were already present in the low network (many disappear/stabilise).
- A high `focus_deg_high` means most partners are present in the high network (many new/stabilise).
- `focus_deg_low + focus_deg_high > L1_n_nodes` is **expected** because SIGN_CHANGE, STRENGTHEN, and WEAKEN edges contribute to both counts (they are present in both conditions).

### L2 equivalents

`L2_n_edges_low` / `L2_n_edges_high` apply the same logic to the L2 outer-layer edges (L1↔L1 and L1→L2).

---

## 7. Complete Sanity Check Table

| Identity | Scope | Exactness |
|----------|-------|-----------|
| `D + N + SC + ST + WK + U = n_sig_edges_diff` | Global | Exact |
| `D + N + SC + ST + WK + U = L1_n_nodes` | L1 | Exact |
| `L1_rewire = D + N + SC` | L1 | Exact (definition) |
| `L1_n_edges_low = D + SC + ST + WK` | L1 | Exact |
| `L1_n_edges_high = N + SC + ST + WK` | L1 | Exact |
| `focus_deg_low == L1_n_edges_low` | L1 | Always true |
| `L1_n_edges_low + L1_n_edges_high = L1_n_nodes + SC + ST + WK` | L1 | Exact (overlap = SC+ST+WK) |
| `L1_n_edges_low ≈ L1_n_nodes - L1_n_new` | L1 | Approx (exact when U=0) |
| `L1_n_edges_high ≈ L1_n_nodes - L1_n_disappear` | L1 | Approx (exact when U=0) |

*Abbreviations: D=disappear, N=new, SC=sign_change, ST=strengthen, WK=weaken, U=unchanged.*

---

## 8. Output Columns in `rewiring_hubs.tsv`

| Column | Type | Description |
|--------|------|-------------|
| `n_sig_edges_diff` | int | Total significant differential edges for this focus gene |
| `n_disappear` | int | Global: edges present in low only |
| `n_new` | int | Global: edges present in high only |
| `n_sign_change` | int | Global: edges with sign flip |
| `n_strengthen` | int | Global: edges that strengthen (pure, excludes sign_change) |
| `n_weaken` | int | Global: edges that weaken (pure, excludes sign_change) |
| `focus_deg_low` | int | L1: how many differential partners present in low network (= `L1_n_edges_low`) |
| `focus_deg_high` | int | L1: how many differential partners present in high network (= `L1_n_edges_high`) |
| `L1_n_nodes` | int | L1: number of direct differential partners (= number of L1 edges) |
| `L1_n_disappear` | int | L1: disappear edges among direct partners |
| `L1_n_new` | int | L1: new edges among direct partners |
| `L1_n_sign_chg` | int | L1: sign-change edges among direct partners |
| `L1_n_strengthen` | int | L1: strengthen edges among direct partners (pure, excludes sign_change) |
| `L1_n_weaken` | int | L1: weaken edges among direct partners (pure, excludes sign_change) |
| `L1_rewire` | int | L1: `L1_n_disappear + L1_n_new + L1_n_sign_chg` |
| `L1_n_edges_low` | int | L1: edges present in low network (= `focus_deg_low`) |
| `L1_n_edges_high` | int | L1: edges present in high network (= `focus_deg_high`) |
| `L1_frac_rewire` | float | `L1_rewire / L1_n_nodes` |
| `L1_clique_density` | float | Fraction of possible L1↔L1 edges present in diff network |
| `L2_n_nodes` | int | L2: number of pure 2nd-hop nodes (not focus gene or L1 partners) |
| `L2_n_edges` | int | L2: total outer-layer edges (L1↔L1 + L1→L2) |
| `L2_n_edges_low` | int | L2: outer-layer edges present in low network |
| `L2_n_edges_high` | int | L2: outer-layer edges present in high network |
| `L2_rewire` | int | L2: `L2_n_disappear + L2_n_new + L2_n_sign_chg` |
| `L2L1_deg` | float | `L2_n_nodes / L1_n_nodes` — expansion ratio of neighbourhood |
| `L2L1_rewire` | float | `L2_rewire / L1_rewire` — relative rewiring in outer vs inner layer |
| `HL_conn_L1` | float | `L1_conn_high / L1_conn_low` — connectivity ratio high/low at L1 |

---

## 9. Edge Selection and Its Effect on UNCHANGED

The pipeline supports two `edge_selection` modes:

| Mode | Edges selected | UNCHANGED edges |
|------|---------------|-----------------|
| `sig_differential` (default) | Edges where delta is statistically significant (CI excludes 0) | Possible — edge with significant delta but weak individual correlations |
| `sig_edges` | Edges significant in at least one condition (`sig_low OR sig_high`) | Very rare — would require an edge to be sig but below `corr_threshold` |

For most analyses with the default mode, UNCHANGED edges are rare and `n_disappear + n_new + n_sign_change + n_strengthen + n_weaken ≈ n_sig_edges_diff`.

---

## Related Reading

- [GUIDE-02-Network-Metrics.md](GUIDE-02-Network-Metrics.md) — Global network topology metrics (degree distribution, clustering, scale-free properties)
- [REFERENCE-01-Statistical-Methods.md](REFERENCE-01-Statistical-Methods.md) — Bootstrap significance testing and the `sig_low`/`sig_high`/`sig_differential` flags
- [GUIDE-01-Complete-Workflow.md](GUIDE-01-Complete-Workflow.md) — Full pipeline and how Stage 3/3b produce these metrics
