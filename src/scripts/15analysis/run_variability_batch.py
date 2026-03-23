#!/usr/bin/env python3
# ==============================================================================
# run_variability_batch.py
#
# Orchestrates plot_gene_mad_variability.R and plot_sample_variability.R
# across multiple conditions, genes, and sample metrics. Aggregates the
# resulting CSV outputs into a single Excel summary table.
#
# Usage:
#   # Dry run — print all Rscript commands without executing:
#   python run_variability_batch.py --config config/variability_batch.yaml --dry-run
#
#   # Full run (4 parallel workers):
#   python run_variability_batch.py --config config/variability_batch.yaml --jobs 4
#
#   # Skip existing PDFs, re-run only missing ones:
#   python run_variability_batch.py --config config/variability_batch.yaml --skip-existing
#
#   # Re-generate Excel from existing CSVs (skip all Rscript calls):
#   python run_variability_batch.py --config config/variability_batch.yaml --excel-only
#
# Config YAML keys (see config/variability_batch.yaml):
#   output_dir, low_frac, high_frac, sample_metrics, save_csv, jobs,
#   conditions (name → {expr_file, label}), genes (list of {id, symbol}),
#   rscript_path
#
# Output layout:
#   {output_dir}/{cond}/
#       {gene_id}_mad_variability.pdf / .csv
#       {gene_id}_sample_variability_{metric}.pdf / .csv
#   {output_dir}/variability_summary.xlsx
#       Sheet "gene_mad"   : one row per (focus_gene × condition)
#       Sheet "sample_itv" : one row per (focus_gene × condition × metric)
# ==============================================================================

import argparse
import logging
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
import yaml
from scipy.stats import mannwhitneyu

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Script paths (relative to CWD — matches Snakemake working directory)
# ---------------------------------------------------------------------------
_MAD_SCRIPT    = "code/dronetanalysis/src/scripts/15analysis/plot_gene_mad_variability.R"
_SAMPLE_SCRIPT = "code/dronetanalysis/src/scripts/15analysis/plot_sample_variability.R"


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------
@dataclass
class Job:
    kind: str           # "gene_mad" | "sample_itv"
    cond_name: str
    expr_file: str
    cond_label: str
    gene_id: str
    gene_symbol: Optional[str]
    output_dir: str     # per-condition sub-dir
    low_frac: float
    high_frac: float
    metric: Optional[str] = None   # sample_itv only
    rscript: str = "Rscript"
    save_csv: bool = True

    def pdf_path(self) -> Path:
        if self.kind == "gene_mad":
            return Path(self.output_dir) / f"{self.gene_id}_mad_variability.pdf"
        return Path(self.output_dir) / f"{self.gene_id}_sample_variability_{self.metric}.pdf"

    def csv_path(self) -> Path:
        return self.pdf_path().with_suffix(".csv")

    def build_cmd(self) -> List[str]:
        script = _MAD_SCRIPT if self.kind == "gene_mad" else _SAMPLE_SCRIPT
        cmd = [
            self.rscript, script,
            "--expr-file",       self.expr_file,
            "--focus-gene",      self.gene_id,
            "--output-dir",      self.output_dir,
            "--low-frac",        str(self.low_frac),
            "--high-frac",       str(self.high_frac),
            "--condition-label", self.cond_label,
        ]
        if self.gene_symbol:
            cmd += ["--gene-symbol", self.gene_symbol]
        if self.kind == "sample_itv":
            cmd += ["--summary-metric", self.metric]
        if self.save_csv:
            cmd.append("--save-csv")
        return cmd


# ---------------------------------------------------------------------------
# Config loading and validation
# ---------------------------------------------------------------------------
def load_config(path: str) -> dict:
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    return cfg


def validate_config(cfg: dict, excel_only: bool) -> None:
    for cond, info in cfg["conditions"].items():
        if not os.path.exists(info["expr_file"]):
            raise FileNotFoundError(
                f"Condition '{cond}' expr_file not found: {info['expr_file']}"
            )
    if not cfg.get("genes"):
        raise ValueError("Config 'genes' list is empty.")
    if not excel_only and not cfg.get("save_csv", True):
        raise ValueError(
            "save_csv must be true for Excel aggregation. "
            "Set save_csv: true in config or use --excel-only."
        )


# ---------------------------------------------------------------------------
# Build job list
# ---------------------------------------------------------------------------
def build_jobs(cfg: dict, output_dir: str) -> List[Job]:
    jobs: List[Job] = []
    rscript   = cfg.get("rscript_path", "Rscript")
    low_frac  = cfg.get("low_frac", 0.2)
    high_frac = cfg.get("high_frac", 0.2)
    save_csv  = cfg.get("save_csv", True)
    metrics   = cfg.get("sample_metrics", ["median", "mean", "sum"])

    for cond_name, cond_info in cfg["conditions"].items():
        cond_dir = os.path.join(output_dir, cond_name)
        for gene in cfg["genes"]:
            gene_id     = gene["id"]
            gene_symbol = gene.get("symbol") or None  # treat empty string as None

            # Gene-level MAD job
            jobs.append(Job(
                kind        = "gene_mad",
                cond_name   = cond_name,
                expr_file   = cond_info["expr_file"],
                cond_label  = cond_info.get("label", cond_name),
                gene_id     = gene_id,
                gene_symbol = gene_symbol,
                output_dir  = cond_dir,
                low_frac    = low_frac,
                high_frac   = high_frac,
                rscript     = rscript,
                save_csv    = save_csv,
            ))

            # Sample-level ITV jobs (one per metric)
            for metric in metrics:
                jobs.append(Job(
                    kind        = "sample_itv",
                    cond_name   = cond_name,
                    expr_file   = cond_info["expr_file"],
                    cond_label  = cond_info.get("label", cond_name),
                    gene_id     = gene_id,
                    gene_symbol = gene_symbol,
                    output_dir  = cond_dir,
                    low_frac    = low_frac,
                    high_frac   = high_frac,
                    metric      = metric,
                    rscript     = rscript,
                    save_csv    = save_csv,
                ))
    return jobs


# ---------------------------------------------------------------------------
# Run a single job
# ---------------------------------------------------------------------------
def run_job(job: Job, skip_existing: bool, dry_run: bool) -> str:
    """Execute one R script job. Returns 'run', 'skipped', or 'failed'."""
    pdf = job.pdf_path()

    if skip_existing and pdf.exists():
        log.debug("SKIP (exists): %s", pdf)
        return "skipped"

    os.makedirs(job.output_dir, exist_ok=True)
    cmd = job.build_cmd()
    cmd_str = " ".join(cmd)

    if dry_run:
        print(cmd_str)
        return "dry_run"

    log.info("RUN [%s | %s | %s]: %s",
             job.cond_name, job.gene_id,
             job.metric or "mad", pdf.name)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.error("FAILED: %s\n  stdout: %s\n  stderr: %s",
                  cmd_str,
                  result.stdout[-500:] if result.stdout else "",
                  result.stderr[-500:] if result.stderr else "")
        return "failed"
    return "run"


# ---------------------------------------------------------------------------
# Run all jobs (optionally parallel)
# ---------------------------------------------------------------------------
def run_all_jobs(
    jobs: List[Job],
    skip_existing: bool,
    dry_run: bool,
    n_workers: int,
) -> Dict[str, int]:
    counts = {"run": 0, "skipped": 0, "failed": 0, "dry_run": 0}

    if n_workers <= 1 or dry_run:
        for job in jobs:
            status = run_job(job, skip_existing, dry_run)
            counts[status] += 1
    else:
        with ThreadPoolExecutor(max_workers=n_workers) as pool:
            futures = {
                pool.submit(run_job, job, skip_existing, dry_run): job
                for job in jobs
            }
            for future in as_completed(futures):
                status = future.result()
                counts[status] += 1

    return counts


# ---------------------------------------------------------------------------
# Compute per-group summary stats from a long-format CSV
# ---------------------------------------------------------------------------
def _group_stats(df: pd.DataFrame, value_col: str) -> dict:
    """Compute mean, median, Wilcoxon p from a LOW/HIGH long-format frame."""
    low  = df.loc[df["group"] == "LOW",  value_col].dropna()
    high = df.loc[df["group"] == "HIGH", value_col].dropna()

    if len(low) < 2 or len(high) < 2:
        p_val = float("nan")
    else:
        _, p_val = mannwhitneyu(low, high, alternative="two-sided")

    return {
        "n_low":        len(low),
        "n_high":       len(high),
        "mean_low":     float(low.mean()),
        "mean_high":    float(high.mean()),
        "median_low":   float(low.median()),
        "median_high":  float(high.median()),
        "wilcoxon_p":   p_val,
    }


# ---------------------------------------------------------------------------
# Aggregate CSVs → Excel
# ---------------------------------------------------------------------------
def aggregate_to_excel(jobs: List[Job], excel_path: str) -> None:
    mad_rows   = []
    itv_rows   = []

    for job in jobs:
        csv = job.csv_path()
        if not csv.exists():
            log.warning("CSV not found (skipping): %s", csv)
            continue

        try:
            df = pd.read_csv(csv)
        except Exception as exc:
            log.warning("Cannot read CSV %s: %s", csv, exc)
            continue

        symbol = job.gene_symbol or job.gene_id

        if job.kind == "gene_mad":
            stats = _group_stats(df, "mad_value")
            mad_rows.append({
                "focus_gene":    job.gene_id,
                "symbol":        symbol,
                "condition":     job.cond_name,
                "n_genes":       stats["n_low"],   # n_low == n_high for gene MAD
                **{k: stats[k] for k in
                   ("mean_low", "mean_high", "median_low", "median_high", "wilcoxon_p")},
            })

        else:  # sample_itv
            stats = _group_stats(df, "itv")
            itv_rows.append({
                "focus_gene":    job.gene_id,
                "symbol":        symbol,
                "condition":     job.cond_name,
                "metric":        job.metric,
                **stats,
            })

    if not mad_rows and not itv_rows:
        log.error("No CSVs were readable — Excel not written.")
        return

    df_mad = pd.DataFrame(mad_rows)
    df_itv = pd.DataFrame(itv_rows)

    # Sort for readability
    if not df_mad.empty:
        df_mad.sort_values(["condition", "focus_gene"], inplace=True)
    if not df_itv.empty:
        df_itv.sort_values(["condition", "focus_gene", "metric"], inplace=True)

    os.makedirs(os.path.dirname(excel_path) or ".", exist_ok=True)
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        df_mad.to_excel(writer, sheet_name="gene_mad",   index=False)
        df_itv.to_excel(writer, sheet_name="sample_itv", index=False)

    log.info("Excel written: %s  (gene_mad: %d rows, sample_itv: %d rows)",
             excel_path, len(df_mad), len(df_itv))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Batch variability analysis across conditions and genes.\n"
            "Calls plot_gene_mad_variability.R and plot_sample_variability.R\n"
            "for all configured combinations, then aggregates to Excel."
        )
    )
    p.add_argument("--config",       required=True, help="Path to variability_batch.yaml")
    p.add_argument("--output-dir",   default=None,  help="Override config output_dir")
    p.add_argument("--dry-run",      action="store_true",
                   help="Print Rscript commands without executing")
    p.add_argument("--skip-existing", action="store_true",
                   help="Skip if output PDF already exists")
    p.add_argument("--excel-only",   action="store_true",
                   help="Skip R scripts; aggregate existing CSVs → Excel")
    p.add_argument("--jobs",         type=int, default=None,
                   help="Parallel Rscript workers (overrides config)")
    return p.parse_args()


def main():
    args = parse_args()

    cfg = load_config(args.config)
    validate_config(cfg, excel_only=args.excel_only)

    output_dir = args.output_dir or cfg.get("output_dir", "results/variability_combined")
    n_workers  = args.jobs if args.jobs is not None else cfg.get("jobs", 1)

    jobs = build_jobs(cfg, output_dir)

    mad_jobs = [j for j in jobs if j.kind == "gene_mad"]
    itv_jobs = [j for j in jobs if j.kind == "sample_itv"]
    log.info("Jobs to run: %d gene-MAD  +  %d sample-ITV  =  %d total",
             len(mad_jobs), len(itv_jobs), len(jobs))

    # --- Run R scripts ---
    if not args.excel_only:
        counts = run_all_jobs(jobs, args.skip_existing, args.dry_run, n_workers)
        if not args.dry_run:
            log.info("Run summary — ran: %d | skipped: %d | failed: %d",
                     counts["run"], counts["skipped"], counts["failed"])
            if counts["failed"] > 0:
                log.warning("%d job(s) failed; Excel may be incomplete.", counts["failed"])

    # --- Aggregate to Excel ---
    if not args.dry_run:
        excel_path = os.path.join(output_dir, "variability_summary.xlsx")
        aggregate_to_excel(jobs, excel_path)

    log.info("Done.")


if __name__ == "__main__":
    main()
