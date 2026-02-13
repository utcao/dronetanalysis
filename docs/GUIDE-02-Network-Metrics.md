# Network Topology Metrics Guide

## Understanding Your Current Metrics

### Current Metrics Explained

Your current `compute_global_topology()` function returns:

```python
{
    "n_nodes_total": 20000,      # Total genes
    "n_nodes_active": 5000,      # Genes with at least 1 connection
    "n_edges": 10000,            # Number of connections
    "density": 0.00005,          # How "connected" the network is
    "avg_degree": 1.0,           # Average connections per gene
    "max_degree": 150,           # Most connected gene
    "median_degree": 2.0,        # Median connections (active genes only)
    "degrees": array([...])      # Connections for each gene
}
```

### What These Metrics Tell You

| Metric | What It Means | Good Value | Interpretation |
|--------|---------------|------------|----------------|
| **n_nodes_active** | Genes participating in network | Depends on FDR | More = denser network |
| **n_edges** | Number of significant co-expressions | Thousands | More = more signal |
| **density** | % of possible edges present | 0.001-0.01 | Higher = more connected |
| **avg_degree** | Avg connections per gene | 5-50 | Higher = hub-rich |
| **max_degree** | Largest hub size | 50-500 | Identifies super-hubs |
| **median_degree** | Typical gene connectivity | 2-10 | More robust than mean |

---

## Additional Important Metrics You Should Add

### 1. **Scale-Free Properties** ⭐ VERY IMPORTANT

**What:** Do you have few super-hubs and many low-connected genes (power-law)?

**Why:** Most biological networks are scale-free (like social networks, internet)

```python
def compute_scale_free_metrics(degrees: np.ndarray) -> dict:
    """
    Test if network follows scale-free topology.
    
    Biological networks typically have γ = 2-3
    """
    from scipy import stats
    
    active_degrees = degrees[degrees > 0]
    if len(active_degrees) < 10:
        return {"power_law_exponent": np.nan, "is_scale_free": False}
    
    # Degree distribution
    unique_degrees, counts = np.unique(active_degrees, return_counts=True)
    mask = counts >= 2
    
    if mask.sum() < 5:
        return {"power_law_exponent": np.nan, "is_scale_free": False}
    
    # Fit power law: log(P(k)) = -γ * log(k)
    log_k = np.log10(unique_degrees[mask])
    log_p = np.log10(counts[mask])
    
    slope, intercept, r_value, _, _ = stats.linregress(log_k, log_p)
    gamma = -slope
    
    # Scale-free if γ ∈ [2,3.5] and good fit
    is_scale_free = (2.0 <= gamma <= 3.5) and (r_value**2 > 0.8)
    
    return {
        "power_law_exponent": gamma,
        "power_law_r_squared": r_value**2,
        "is_scale_free": is_scale_free,
    }
```

**Interpretation:**
- **γ = 2-3**: Scale-free (typical for biological networks)
- **γ < 2**: Super-hubs dominate
- **γ > 3**: More uniform distribution
- **R² > 0.8**: Good fit

---

### 2. **Clustering Coefficient** ⭐ IMPORTANT

**What:** Do your neighbors connect to each other (triangles)?

**Why:** High clustering = functional modules/pathways

```python
def compute_clustering_coefficient(gene_i, gene_j, n_genes) -> float:
    """
    Measure tendency to form triangles (modules).
    
    High clustering = genes in same pathway cluster together
    """
    from collections import defaultdict
    
    adj = defaultdict(set)
    for i, j in zip(gene_i, gene_j):
        adj[i].add(j)
        adj[j].add(i)
    
    local_clustering = []
    for node in range(n_genes):
        neighbors = list(adj[node])
        k = len(neighbors)
        
        if k < 2:
            continue
        
        # Count triangles
        triangles = sum(
            1 for i, n1 in enumerate(neighbors)
            for n2 in neighbors[i+1:]
            if n2 in adj[n1]
        )
        
        # Local clustering
        possible = k * (k - 1) / 2
        local_clustering.append(triangles / possible)
    
    return np.mean(local_clustering) if local_clustering else 0.0
```

**Interpretation:**
- **C > 0.3**: Highly modular (genes form tight groups)
- **C = 0.1-0.3**: Moderate modularity (typical)
- **C < 0.1**: Sparse, non-modular

---

### 3. **Connected Components** ⭐ IMPORTANT

**What:** How many separate "islands" in your network?

**Why:** Fragmented network = poor connectivity, multiple independent modules

```python
def compute_connected_components(gene_i, gene_j, n_genes) -> dict:
    """Find disconnected subnetworks."""
    from collections import defaultdict, deque
    
    adj = defaultdict(set)
    for i, j in zip(gene_i, gene_j):
        adj[i].add(j)
        adj[j].add(i)
    
    visited = set()
    components = []
    
    for start in range(n_genes):
        if start in visited or len(adj[start]) == 0:
            continue
        
        # BFS
        component = set()
        queue = deque([start])
        
        while queue:
            node = queue.popleft()
            if node in visited:
                continue
            visited.add(node)
            component.add(node)
            queue.extend(n for n in adj[node] if n not in visited)
        
        components.append(len(component))
    
    components = sorted(components, reverse=True)
    
    return {
        "n_components": len(components),
        "largest_component_size": components[0] if components else 0,
        "largest_component_fraction": components[0] / n_genes if components else 0,
    }
```

**Interpretation:**
- **n_components = 1**: Fully connected (ideal)
- **n_components = 2-10**: Few large modules (good)
- **n_components > 100**: Very fragmented (bad)
- **Largest fraction > 0.8**: Most genes in one network (good)

---

### 4. **Assortativity** (Hub Connectivity)

**What:** Do hubs connect to other hubs, or to peripheral nodes?

```python
def compute_assortativity(gene_i, gene_j, degrees) -> float:
    """
    Pearson correlation of degrees at edge endpoints.
    
    r > 0: Assortative (hubs ↔ hubs)
    r < 0: Disassortative (hubs ↔ periphery)
    """
    from scipy.stats import pearsonr
    
    if len(gene_i) == 0:
        return 0.0
    
    deg_i = degrees[gene_i]
    deg_j = degrees[gene_j]
    
    r, _ = pearsonr(deg_i, deg_j)
    return r
```

**Interpretation:**
- **r > 0.3**: Hubs connect to hubs (hierarchical core)
- **r ≈ 0**: Random mixing
- **r < -0.3**: Hubs connect periphery (star-like)

---

## Comprehensive Metrics Function

Here's the **improved version** with all important metrics:

```python
def compute_comprehensive_topology(
    gene_i: np.ndarray,
    gene_j: np.ndarray,
    values: np.ndarray,
    n_genes: int,
    corr_threshold: float = 0.1,
    gene_names: Optional[list] = None,
) -> dict:
    """
    Compute comprehensive network topology metrics.
    
    Returns complete characterization of network structure.
    """
    from collections import defaultdict
    from scipy import stats
    
    # Filter edges
    mask = np.abs(values) >= corr_threshold
    gi, gj = gene_i[mask], gene_j[mask]
    n_edges = len(gi)
    
    if n_edges == 0:
        return _empty_topology(n_genes)
    
    # Build adjacency
    adj = defaultdict(set)
    for i, j in zip(gi, gj):
        adj[i].add(j)
        adj[j].add(i)
    
    # === BASIC METRICS ===
    nodes_active = set(gi) | set(gj)
    n_nodes_active = len(nodes_active)
    degrees = np.array([len(adj[g]) for g in range(n_genes)], dtype=np.int32)
    
    density = 2 * n_edges / (n_genes * (n_genes - 1)) if n_genes > 1 else 0
    avg_degree = 2 * n_edges / n_genes if n_genes > 0 else 0
    
    active_degrees = degrees[degrees > 0]
    max_degree = int(degrees.max())
    median_degree = float(np.median(active_degrees))
    
    # === SCALE-FREE PROPERTIES ===
    scale_free = compute_scale_free_metrics(degrees)
    
    # === CLUSTERING ===
    clustering = compute_clustering_coefficient(gi, gj, n_genes)
    
    # === COMPONENTS ===
    components = compute_connected_components(gi, gj, n_genes)
    
    # === ASSORTATIVITY ===
    assortativity = compute_assortativity(gi, gj, degrees)
    
    # === TOP HUBS ===
    top_indices = np.argsort(degrees)[::-1][:20]
    if gene_names is not None:
        top_hubs = [
            {"gene": gene_names[i], "degree": int(degrees[i])}
            for i in top_indices if degrees[i] > 0
        ]
    else:
        top_hubs = [
            {"gene": f"gene_{i}", "degree": int(degrees[i])}
            for i in top_indices if degrees[i] > 0
        ]
    
    # === DEGREE PERCENTILES ===
    percentiles = {
        f"degree_p{p}": float(np.percentile(active_degrees, p))
        for p in [25, 50, 75, 90, 95, 99]
    }
    
    return {
        # Basic
        "n_nodes_total": n_genes,
        "n_nodes_active": n_nodes_active,
        "n_nodes_inactive": n_genes - n_nodes_active,
        "fraction_active": n_nodes_active / n_genes,
        "n_edges": n_edges,
        "density": density,
        
        # Degree distribution
        "avg_degree": avg_degree,
        "median_degree": median_degree,
        "max_degree": max_degree,
        "min_degree": int(active_degrees.min()) if len(active_degrees) > 0 else 0,
        "std_degree": float(active_degrees.std()) if len(active_degrees) > 0 else 0,
        **percentiles,
        
        # Scale-free
        "power_law_exponent": scale_free["power_law_exponent"],
        "power_law_r_squared": scale_free.get("power_law_r_squared", np.nan),
        "is_scale_free": scale_free["is_scale_free"],
        
        # Clustering
        "global_clustering": clustering,
        
        # Components
        "n_components": components["n_components"],
        "largest_component_size": components["largest_component_size"],
        "largest_component_fraction": components["largest_component_fraction"],
        
        # Assortativity
        "assortativity": assortativity,
        
        # Top genes
        "top_hubs": top_hubs[:10],
        
        # Raw
        "degrees": degrees,
    }
```

---

## Better Reporting

### 1. **Summary Table (Human-Readable)**

```python
def write_topology_summary(topo_low, topo_high, topo_diff, output_file):
    """
    Create a comparison table across networks.
    
    Output: topology_summary.tsv
    """
    metrics = [
        ("Active Nodes", "n_nodes_active"),
        ("Edges", "n_edges"),
        ("Density", "density"),
        ("Avg Degree", "avg_degree"),
        ("Median Degree", "median_degree"),
        ("Max Degree", "max_degree"),
        ("Clustering", "global_clustering"),
        ("Components", "n_components"),
        ("Largest Component (%)", "largest_component_fraction"),
        ("Power Law γ", "power_law_exponent"),
        ("Scale-Free?", "is_scale_free"),
        ("Assortativity", "assortativity"),
    ]
    
    with open(output_file, 'w') as f:
        f.write("Metric\tLow_Network\tHigh_Network\tDifferential\n")
        for label, key in metrics:
            low_val = format_value(topo_low.get(key))
            high_val = format_value(topo_high.get(key))
            diff_val = format_value(topo_diff.get(key))
            f.write(f"{label}\t{low_val}\t{high_val}\t{diff_val}\n")
    
    print(f"✓ Saved: {output_file}")


def format_value(val):
    """Format for display."""
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return "N/A"
    elif isinstance(val, bool):
        return "Yes" if val else "No"
    elif isinstance(val, float):
        if abs(val) < 0.01:
            return f"{val:.2e}"
        else:
            return f"{val:.3f}"
    elif isinstance(val, int):
        return f"{val:,}"
    return str(val)
```

### 2. **Visual Summary (PDF)**

```python
import matplotlib.pyplot as plt
import numpy as np

def create_topology_comparison_plot(topo_low, topo_high, topo_diff, output_file):
    """
    Create 3-panel comparison figure.
    
    Panel 1: Degree distributions (log-log)
    Panel 2: Metric comparison (bar chart)
    Panel 3: Top hubs comparison
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # Panel 1: Degree distributions
    for topo, label, color in [
        (topo_low, "Low", "blue"),
        (topo_high, "High", "red"),
        (topo_diff, "Diff", "green"),
    ]:
        degrees = topo["degrees"][topo["degrees"] > 0]
        if len(degrees) > 0:
            axes[0].hist(degrees, bins=50, alpha=0.5, label=label, 
                        color=color, edgecolor='black')
    
    axes[0].set_xlabel("Degree (log scale)")
    axes[0].set_ylabel("Frequency (log scale)")
    axes[0].set_title("Degree Distribution Comparison")
    axes[0].set_xscale('log')
    axes[0].set_yscale('log')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # Panel 2: Metric comparison
    metrics = ['density', 'avg_degree', 'global_clustering', 'assortativity']
    labels = ['Density', 'Avg Degree', 'Clustering', 'Assortativity']
    
    x = np.arange(len(labels))
    width = 0.25
    
    low_vals = [topo_low.get(m, 0) for m in metrics]
    high_vals = [topo_high.get(m, 0) for m in metrics]
    diff_vals = [topo_diff.get(m, 0) for m in metrics]
    
    # Normalize for visibility
    max_vals = [max(abs(l), abs(h), abs(d)) for l, h, d in zip(low_vals, high_vals, diff_vals)]
    low_norm = [l/m if m > 0 else 0 for l, m in zip(low_vals, max_vals)]
    high_norm = [h/m if m > 0 else 0 for h, m in zip(high_vals, max_vals)]
    diff_norm = [d/m if m > 0 else 0 for d, m in zip(diff_vals, max_vals)]
    
    axes[1].bar(x - width, low_norm, width, label='Low', color='blue', alpha=0.7)
    axes[1].bar(x, high_norm, width, label='High', color='red', alpha=0.7)
    axes[1].bar(x + width, diff_norm, width, label='Diff', color='green', alpha=0.7)
    
    axes[1].set_ylabel('Normalized Value')
    axes[1].set_title('Network Metrics Comparison')
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=45, ha='right')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3, axis='y')
    
    # Panel 3: Top hubs
    top_n = 10
    low_hubs = topo_low.get("top_hubs", [])[:top_n]
    high_hubs = topo_high.get("top_hubs", [])[:top_n]
    diff_hubs = topo_diff.get("top_hubs", [])[:top_n]
    
    y_pos = np.arange(top_n)
    
    # Show all three side by side
    axes[2].barh(y_pos, [h["degree"] for h in low_hubs] + [0]*(top_n-len(low_hubs)), 
                alpha=0.3, color='blue', label='Low')
    axes[2].barh(y_pos, [h["degree"] for h in high_hubs] + [0]*(top_n-len(high_hubs)), 
                alpha=0.3, color='red', label='High')
    axes[2].barh(y_pos, [h["degree"] for h in diff_hubs] + [0]*(top_n-len(diff_hubs)), 
                alpha=0.3, color='green', label='Diff')
    
    axes[2].set_yticks(y_pos)
    axes[2].set_yticklabels([h.get("gene", "") for h in diff_hubs] + [""]*(top_n-len(diff_hubs)))
    axes[2].invert_yaxis()
    axes[2].set_xlabel('Degree')
    axes[2].set_title('Top Hub Genes (Differential)')
    axes[2].legend()
    axes[2].grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Saved: {output_file}")
```

---

## What to Report in Your Paper

### **Essential Metrics (Always)**
1. ✅ **n_nodes_active / n_edges** - Network size
2. ✅ **avg_degree / median_degree** - Connectivity level
3. ✅ **max_degree** - Largest hub
4. ✅ **top_hubs** (top 10 genes) - Biological interpretation

### **Important (Recommended)**
5. ✅ **global_clustering** - Modularity
6. ✅ **n_components** - Connectivity
7. ✅ **power_law_exponent** - Scale-free property
8. ✅ **is_scale_free** (Yes/No)

### **Nice to Have**
9. ⚠️ **assortativity** - Hub organization
10. ⚠️ **density** - Overall connectivity

---

## Example Report Output

```
NETWORK TOPOLOGY SUMMARY
========================

Low Expression Network:
  Nodes:          5,234 / 20,000 (26.2% active)
  Edges:          12,345
  Density:        0.00045
  Avg degree:     4.72 (median: 3.0)
  Max degree:     234 (BRCA1)
  Clustering:     0.34 (modular)
  Components:     3 (98.5% in largest)
  Scale-free:     Yes (γ=2.4, R²=0.92)
  Assortativity:  0.12 (slight hub clustering)

High Expression Network:
  Nodes:          4,892 / 20,000 (24.5% active)
  Edges:          10,234
  Density:        0.00043
  Avg degree:     4.19 (median: 2.0)
  Max degree:     198 (TP53)
  Clustering:     0.31 (modular)
  Components:     5 (95.2% in largest)
  Scale-free:     Yes (γ=2.6, R²=0.89)
  Assortativity:  0.09 (slight hub clustering)

Differential Network (Rewired Edges):
  Edges:          2,345 (rewired)
  Avg degree:     1.36
  Clustering:     0.12 (less modular)
  Components:     12 (fragmented)
  Scale-free:     No (γ=1.8, R²=0.65)
  
Top Rewiring Hubs:
  1. FOXM1    (degree: 87)
  2. MYC      (degree: 65)
  3. CCNB1    (degree: 54)
  ...
```

Would you like me to integrate these into your Stage 3 script?

---

## Related Reading

- [GUIDE-01-Complete-Workflow.md](GUIDE-01-Complete-Workflow.md) - Full pipeline workflow and usage
- [FIX-01-Critical-Issues-Summary.md](FIX-01-Critical-Issues-Summary.md) - Summary of enhanced metrics and fixes
- [REFERENCE-01-Statistical-Methods.md](REFERENCE-01-Statistical-Methods.md) - Statistical methodology reference