#!/usr/bin/env python3
"""
Gene expression boxplot/violin plot with individual selection and dot highlighting.

Input CSVs: rows = genes, columns = individuals, values = logCPM (or similar).

Individual selection (--individual_sel):
  1.0      → use all individuals from each file
  0 < x < 1 → select top x% AND bottom x% of individuals per file, ranked by
               expression of --gene (or mean across all genes if not specified)

X-axis:
  Multiple files → x-axis shows group names (one violin per file)
  Single file + individual_sel < 1 → x-axis shows "top X%" and "bottom X%"

Dot display (all selected dots shown by default):
  Non-highlighted top-selected  : salmon
  Non-highlighted bottom-selected: lightblue
  Non-highlighted (when all selected) : lightgray

Dot highlighting (--highlight_sel):
  0 or 1   → no highlighting
  0 < x < 1 → within selected individuals, highlight top x% (red) and bottom x% (blue)
               This is computed AFTER individual selection.

Statistics:
  Pairwise Wilcox (Mann-Whitney U, two-sided) between x-axis categories.
  BH-FDR correction. Exact p-values displayed.

Output:
  results/boxplot_violinplot/{prefix}.pdf / .svg / .png
  results/boxplot_violinplot/{prefix}_stats.csv

Example:
  python src/scripts/plot_expr_violin.py \\
      --files dataset/raw/logCPM_Ctrl_Dros.csv dataset/raw/logCPM_HS_Dros.csv \\
      --groups CT HS \\
      --individual-sel 0.1 \\
      --highlight-sel 0.05 \\
      --output-prefix CT_vs_HS_top10pct
"""

import argparse
import os
import sys
import warnings

import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns

warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# Pure-numpy stats (replaces scipy + statsmodels to avoid DLL issues)
# ---------------------------------------------------------------------------

def _rankdata(a):
    """Average-rank with tie handling (mirrors scipy.stats.rankdata)."""
    a = np.asarray(a, dtype=float)
    n = len(a)
    order = np.argsort(a, kind="mergesort")
    ranks = np.empty(n, dtype=float)
    ranks[order] = np.arange(1, n + 1, dtype=float)
    # resolve ties
    a_sorted = a[order]
    i = 0
    while i < n:
        j = i
        while j < n - 1 and a_sorted[j] == a_sorted[j + 1]:
            j += 1
        if j > i:
            avg = (i + j + 2) / 2.0   # average of 1-based ranks i+1 .. j+1
            ranks[order[i : j + 1]] = avg
        i = j + 1
    return ranks


def _norm_sf(z):
    """Survival function of standard normal (P(Z > z))."""
    return 0.5 * math.erfc(z / math.sqrt(2))


def mannwhitneyu(x, y):
    """Two-sided Mann-Whitney U test, normal approximation (no continuity correction)."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    n1, n2 = len(x), len(y)
    combined = np.concatenate([x, y])
    ranks = _rankdata(combined)
    R1 = ranks[:n1].sum()
    U1 = R1 - n1 * (n1 + 1) / 2.0
    U2 = n1 * n2 - U1
    U = min(U1, U2)
    mu = n1 * n2 / 2.0
    sigma = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
    if sigma == 0:
        return U, 1.0
    z = abs(U - mu) / sigma
    p = 2.0 * _norm_sf(z)
    return U, min(p, 1.0)


def bh_correction(pvals):
    """Benjamini-Hochberg FDR correction. Returns adjusted p-values."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    if n == 0:
        return pvals.copy()
    order = np.argsort(pvals)
    adj = np.empty(n, dtype=float)
    adj[order] = pvals[order] * n / (np.arange(n, dtype=float) + 1)
    # enforce monotonicity: step through sorted order from right
    for i in range(n - 2, -1, -1):
        adj[order[i]] = min(adj[order[i]], adj[order[i + 1]])
    return np.minimum(adj, 1.0)


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Gene expression boxplot/violin plot",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--files", nargs="+", required=True,
                   help="Input CSV files (rows=genes, cols=individuals)")
    p.add_argument("--groups", nargs="+", required=True,
                   help="Group name for each input file (same order as --files)")
    p.add_argument("--individual-sel", type=float, default=1.0,
                   help=(
                       "Individual selection fraction. "
                       "1 = use all samples; "
                       "0 < x < 1 = select top x%% AND bottom x%% per file "
                       "(ranked by expression). Default: 1.0"
                   ))
    p.add_argument("--highlight-sel", type=float, default=0.0,
                   help=(
                       "Dot-highlighting fraction (applied after individual selection). "
                       "0 or 1 = no highlighting; "
                       "0 < x < 1 = highlight top x%% (red) and bottom x%% (blue). "
                       "Default: 0 (no highlighting)"
                   ))
    p.add_argument("--gene", type=str, default=None,
                   help=(
                       "Gene ID to plot (must be a row index in the CSV). "
                       "If omitted, the mean expression across all genes is used."
                   ))
    p.add_argument("--gene-symbol", type=str, default=None,
                   help=(
                       "Human-readable gene symbol (e.g. Hsp83). Matches --gene-symbol in "
                       "plot_gene_mad_variability.R and plot_sample_variability.R. "
                       "When provided with --gene, '{gene_id}_{gene_symbol}' is appended "
                       "to output file names."
                   ))
    p.add_argument("--output-prefix", type=str, default="expr_violin",
                   help="Output file name prefix (no extension). Default: expr_violin")
    p.add_argument("--output-dir", type=str, default=None,
                   help=(
                       "Output directory. Default: results/expr_violin "
                       "(relative to cwd, or absolute path)."
                   ))
    p.add_argument("--figsize", nargs=2, type=float, default=None,
                   metavar=("WIDTH", "HEIGHT"),
                   help="Figure size in inches. Auto-scaled if omitted.")
    p.add_argument("--dpi", type=int, default=150,
                   help="DPI for PNG output. Default: 150")
    p.add_argument("--dot-alpha", type=float, default=0.35,
                   help="Transparency of non-highlighted dots. Default: 0.35")
    p.add_argument("--dot-size", type=float, default=8,
                   help="Marker size of non-highlighted dots. Default: 8")
    p.add_argument("--violin-alpha", type=float, default=0.55,
                   help="Transparency of violin fill. Default: 0.55")
    p.add_argument("--no-violin", action="store_true",
                   help="Skip violin layer (show box plot only)")
    p.add_argument("--no-box", action="store_true",
                   help="Skip box plot layer (show violin only)")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Data loading & individual selection
# ---------------------------------------------------------------------------

def load_file(filepath, group_name, gene, individual_sel, highlight_sel):
    """
    Load one expression CSV and return a tidy DataFrame with columns:
      expression, group, sel_type, highlight_type
    """
    df = pd.read_csv(filepath, index_col=0)
    # df shape: (genes, individuals)

    # Per-individual expression (used for ranking)
    if gene is not None:
        if gene not in df.index:
            sys.exit(f"ERROR: gene '{gene}' not found in {filepath}")
        expr = df.loc[gene].astype(float)
    else:
        expr = df.mean(axis=0).astype(float)

    n = len(expr)

    # ---- Individual selection ----
    if individual_sel == 1.0:
        sel_index = expr.index
        sel_type_map = {idx: "all" for idx in sel_index}
    elif 0.0 < individual_sel < 1.0:
        k = max(1, int(np.floor(n * individual_sel)))
        top_idx = set(expr.nlargest(k).index)
        bot_idx = set(expr.nsmallest(k).index)
        # In rare ties both sets may overlap; top takes priority
        sel_index = pd.Index(sorted(top_idx | bot_idx, key=lambda i: expr[i], reverse=True))
        sel_type_map = {}
        for idx in sel_index:
            if idx in top_idx:
                sel_type_map[idx] = "top"
            else:
                sel_type_map[idx] = "bottom"
    else:
        sys.exit(f"ERROR: --individual_sel must be 1 or strictly in (0, 1), got {individual_sel}")

    selected_expr = expr.loc[sel_index]

    # ---- Highlight selection (within selected) ----
    if highlight_sel == 0.0 or highlight_sel == 1.0:
        highlight_map = {idx: "normal" for idx in sel_index}
    elif 0.0 < highlight_sel < 1.0:
        k_h = max(1, int(np.floor(len(selected_expr) * highlight_sel)))
        hi_top_idx = set(selected_expr.nlargest(k_h).index)
        hi_bot_idx = set(selected_expr.nsmallest(k_h).index)
        highlight_map = {}
        for idx in sel_index:
            if idx in hi_top_idx:
                highlight_map[idx] = "top_highlight"
            elif idx in hi_bot_idx:
                highlight_map[idx] = "bottom_highlight"
            else:
                highlight_map[idx] = "normal"
    else:
        sys.exit(f"ERROR: --highlight_sel must be 0, 1, or strictly in (0, 1), got {highlight_sel}")

    records = pd.DataFrame({
        "individual": sel_index,
        "expression": selected_expr.values,
        "group": group_name,
        "sel_type": [sel_type_map[i] for i in sel_index],
        "highlight_type": [highlight_map[i] for i in sel_index],
    })
    return records


# ---------------------------------------------------------------------------
# Statistical tests
# ---------------------------------------------------------------------------

def run_wilcox(data, x_col, y_col="expression"):
    """Pairwise Mann-Whitney U (two-sided) + BH correction."""
    categories = list(data[x_col].unique())
    pairs = [(a, b) for i, a in enumerate(categories) for b in categories[i + 1:]]

    if not pairs:
        return []

    raw_pvals = []
    for g1, g2 in pairs:
        v1 = data.loc[data[x_col] == g1, y_col].dropna().values
        v2 = data.loc[data[x_col] == g2, y_col].dropna().values
        if len(v1) < 3 or len(v2) < 3:
            raw_pvals.append(np.nan)
        else:
            _, p = mannwhitneyu(v1, v2)
            raw_pvals.append(p)

    valid_mask = [not np.isnan(p) for p in raw_pvals]
    valid_pvals = [p for p, m in zip(raw_pvals, valid_mask) if m]

    if valid_pvals:
        adj_valid = bh_correction(np.array(valid_pvals))
        adj_iter = iter(adj_valid)
        adj_pvals = [float(next(adj_iter)) if m else np.nan
                     for m in valid_mask]
    else:
        adj_pvals = list(raw_pvals)

    return list(zip(pairs, raw_pvals, adj_pvals))


def fmt_p(p):
    if np.isnan(p):
        return "NA"
    if p < 0.0001:
        return f"{p:.2e}"
    return f"{p:.4f}"


# ---------------------------------------------------------------------------
# Significance annotation bars
# ---------------------------------------------------------------------------

def add_significance_bars(ax, test_results, x_order, y_data, y_col="expression"):
    """Draw bracketed p-value annotations above the plot."""
    if not test_results:
        return

    y_max = y_data[y_col].max()
    y_range = y_data[y_col].max() - y_data[y_col].min()
    step = y_range * 0.09
    y_pos = y_max + step * 0.5

    for (g1, g2), raw_p, adj_p in test_results:
        if np.isnan(adj_p):
            continue
        x1 = x_order.index(g1)
        x2 = x_order.index(g2)
        bar_h = step * 0.25
        ax.plot(
            [x1, x1, x2, x2],
            [y_pos, y_pos + bar_h, y_pos + bar_h, y_pos],
            lw=1.2, c="black",
        )
        label = f"p={fmt_p(adj_p)} (BH-FDR)"
        ax.text(
            (x1 + x2) / 2, y_pos + bar_h + step * 0.05,
            label, ha="center", va="bottom", fontsize=8,
        )
        y_pos += step * 1.3

    # Expand y-axis so bars fit
    ax.set_ylim(top=y_pos + step * 0.5)


# ---------------------------------------------------------------------------
# Dot appearance helper
# ---------------------------------------------------------------------------

DOT_STYLES = {
    # (sel_type, highlight_type): (color, alpha_override, size_override, zorder)
    ("all", "normal"):          ("lightgray",  None,  None, 3),
    ("all", "top_highlight"):   ("red",        0.85,  20,   6),
    ("all", "bottom_highlight"):("royalblue",  0.85,  20,   6),
    ("top", "normal"):          ("#f4a582",    None,  None, 3),   # light salmon
    ("top", "top_highlight"):   ("red",        0.90,  20,   6),
    ("top", "bottom_highlight"):("royalblue",  0.90,  20,   6),  # shouldn't happen, fallback
    ("bottom", "normal"):       ("#92c5de",    None,  None, 3),   # light steel blue
    ("bottom", "top_highlight"):("red",        0.90,  20,   6),  # fallback
    ("bottom", "bottom_highlight"):("royalblue",0.90, 20,   6),
}


def scatter_dots(ax, x_pos, y_vals, sel_types, highlight_types,
                 dot_alpha, dot_size, rng):
    """Draw dots for one x-axis category with appropriate styling."""
    df = pd.DataFrame({
        "y": y_vals,
        "sel_type": sel_types,
        "highlight_type": highlight_types,
    })

    for (sel_t, hi_t), (color, a_ov, s_ov, zo) in DOT_STYLES.items():
        mask = (df["sel_type"] == sel_t) & (df["highlight_type"] == hi_t)
        sub = df[mask]
        if sub.empty:
            continue
        y = sub["y"].values
        jitter = rng.uniform(-0.18, 0.18, size=len(y))
        ax.scatter(
            x_pos + jitter, y,
            color=color,
            alpha=a_ov if a_ov is not None else dot_alpha,
            s=s_ov if s_ov is not None else dot_size,
            zorder=zo,
            linewidths=0,
        )


# ---------------------------------------------------------------------------
# Main plot
# ---------------------------------------------------------------------------

def build_legend(individual_sel, highlight_sel):
    handles = []

    if individual_sel == 1.0:
        handles.append(mpatches.Patch(color="lightgray", label="Selected (all)"))
    else:
        pct = f"{individual_sel * 100:.0f}%"
        handles.append(mpatches.Patch(color="#f4a582", label=f"Top {pct} selected"))
        handles.append(mpatches.Patch(color="#92c5de", label=f"Bottom {pct} selected"))

    if 0.0 < highlight_sel < 1.0:
        h_pct = f"{highlight_sel * 100:.0f}%"
        handles.append(
            Line2D([0], [0], marker="o", color="w", markerfacecolor="red",
                   markersize=8, label=f"Top {h_pct} highlighted")
        )
        handles.append(
            Line2D([0], [0], marker="o", color="w", markerfacecolor="royalblue",
                   markersize=8, label=f"Bottom {h_pct} highlighted")
        )

    return handles


def main():
    args = parse_args()

    # Validate
    if len(args.files) != len(args.groups):
        sys.exit("ERROR: --files and --groups must have the same number of entries")

    # ------------------------------------------------------------------ #
    # 1. Load & process all files
    # ------------------------------------------------------------------ #
    all_frames = []
    for filepath, group in zip(args.files, args.groups):
        if not os.path.isfile(filepath):
            sys.exit(f"ERROR: file not found: {filepath}")
        print(f"  Loading {filepath} → group '{group}'")
        frame = load_file(filepath, group, args.gene,
                          args.individual_sel, args.highlight_sel)
        all_frames.append(frame)
        print(f"    {len(frame)} individuals selected")

    data = pd.concat(all_frames, ignore_index=True)

    # ------------------------------------------------------------------ #
    # 2. Determine x-axis column and ordering
    # ------------------------------------------------------------------ #
    single_file = len(args.files) == 1
    use_subgroup_x = single_file and (0.0 < args.individual_sel < 1.0)

    if use_subgroup_x:
        x_col = "sel_type"
        pct = f"{args.individual_sel * 100:.0f}%"
        x_labels = {
            "bottom": f"Bottom {pct}",
            "top": f"Top {pct}",
        }
        x_order_internal = ["bottom", "top"]   # low → high for intuitive layout
        x_order_display  = [x_labels[k] for k in x_order_internal]
        # remap sel_type values for display
        data["x_cat"] = data["sel_type"].map(x_labels)
        x_col = "x_cat"
        x_order = x_order_display
    else:
        x_col = "group"
        x_order = args.groups  # preserve user-specified order

    # ------------------------------------------------------------------ #
    # 3. Statistical tests
    # ------------------------------------------------------------------ #
    test_results = run_wilcox(data, x_col)

    print("\nWilcox tests (BH-FDR corrected):")
    if not test_results:
        print("  No pairwise comparisons (only one category).")
    for (g1, g2), raw_p, adj_p in test_results:
        print(f"  {g1} vs {g2}:  raw p = {fmt_p(raw_p)},  adj p (BH) = {fmt_p(adj_p)}")

    # ------------------------------------------------------------------ #
    # 4. Figure layout
    # ------------------------------------------------------------------ #
    n_x = len(x_order)
    if args.figsize:
        figw, figh = args.figsize
    else:
        figw = max(5, n_x * 2.8 + 1.5)
        figh = 6.5

    palette = sns.color_palette("Set2", n_x)
    color_map = {g: palette[i] for i, g in enumerate(x_order)}

    fig, ax = plt.subplots(figsize=(figw, figh))
    sns.set_style("whitegrid")

    # ---- Violin ----
    if not args.no_violin:
        sns.violinplot(
            data=data, x=x_col, y="expression", order=x_order,
            palette=color_map, inner=None, alpha=args.violin_alpha,
            linewidth=0.8, ax=ax, cut=0,
        )

    # ---- Box ----
    if not args.no_box:
        sns.boxplot(
            data=data, x=x_col, y="expression", order=x_order,
            palette=color_map, width=0.18, fliersize=0,
            boxprops=dict(alpha=0.85, zorder=4),
            whiskerprops=dict(zorder=4),
            medianprops=dict(color="black", linewidth=1.5, zorder=5),
            ax=ax,
        )

    # ---- Dots ----
    rng = np.random.default_rng(42)
    for i, x_val in enumerate(x_order):
        sub = data[data[x_col] == x_val]
        scatter_dots(
            ax, i,
            sub["expression"].values,
            sub["sel_type"].values,
            sub["highlight_type"].values,
            args.dot_alpha, args.dot_size, rng,
        )

    # ---- Significance bars ----
    add_significance_bars(ax, test_results, x_order, data)

    # ---- Labels & aesthetics ----
    if args.gene:
        gene_label = f"{args.gene} ({args.gene_symbol})" if args.gene_symbol else args.gene
        y_label = f"{gene_label} expression (logCPM)"
        title = f"Expression of {gene_label}"
    else:
        y_label = "Mean expression (logCPM, across all genes)"
        title = "Per-individual mean gene expression"

    if args.individual_sel < 1.0:
        title += f" | top/bottom {args.individual_sel * 100:.0f}% individuals"

    ax.set_ylabel(y_label, fontsize=11)
    ax.set_xlabel("", fontsize=11)
    ax.set_title(title, fontsize=12, pad=10)
    ax.set_xticks(range(n_x))
    ax.set_xticklabels(x_order, fontsize=11)

    # n per category as subtitle on x-axis ticks
    tick_labels = []
    for x_val in x_order:
        n_ind = (data[x_col] == x_val).sum()
        tick_labels.append(f"{x_val}\n(n={n_ind})")
    ax.set_xticklabels(tick_labels, fontsize=10)

    legend_handles = build_legend(args.individual_sel, args.highlight_sel)
    if legend_handles:
        ax.legend(handles=legend_handles, fontsize=9,
                  loc="upper right", framealpha=0.85)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()

    # ------------------------------------------------------------------ #
    # 5. Save outputs
    # ------------------------------------------------------------------ #
    if args.output_dir:
        outdir = args.output_dir
    else:
        outdir = os.path.join("results", "expr_violin")

    os.makedirs(outdir, exist_ok=True)

    # Auto-append gene ID + readable name to file names when --gene_name is given
    file_prefix = args.output_prefix
    if args.gene_symbol:
        gene_tag = f"{args.gene}_{args.gene_symbol}" if args.gene else args.gene_symbol
        # sanitise: replace spaces/slashes with underscores
        gene_tag = gene_tag.replace(" ", "_").replace("/", "_")
        file_prefix = f"{file_prefix}_{gene_tag}"

    prefix = os.path.join(outdir, file_prefix)

    fig.savefig(prefix + ".pdf", bbox_inches="tight")
    fig.savefig(prefix + ".svg", bbox_inches="tight")
    fig.savefig(prefix + ".png", bbox_inches="tight", dpi=args.dpi)
    plt.close(fig)

    # Stats CSV
    if test_results:
        stats_rows = []
        for (g1, g2), raw_p, adj_p in test_results:
            n1 = (data[x_col] == g1).sum()
            n2 = (data[x_col] == g2).sum()
            stats_rows.append({
                "group1": g1, "n1": n1,
                "group2": g2, "n2": n2,
                "raw_p": raw_p,
                "adj_p_BH": adj_p,
            })
        stats_df = pd.DataFrame(stats_rows)
        stats_df.to_csv(prefix + "_stats.csv", index=False)
        print(f"\nStats saved: {prefix}_stats.csv")

    print(f"\nPlots saved:")
    print(f"  {prefix}.pdf")
    print(f"  {prefix}.svg")
    print(f"  {prefix}.png")


if __name__ == "__main__":
    main()
