#!/usr/bin/env python3
"""
generate_pseudo_phenotype.py
----------------------------
Generate simulated phenotype data that shares sample identifiers with an
existing expression matrix.

Only the header line of each input file is read, so even very large matrices
(hundreds of thousands of rows) impose no memory cost.

Usage
-----
    python generate_pseudo_phenotype.py --config config.yaml
    python generate_pseudo_phenotype.py --config config.yaml --dataset-tag ctrl_voom
    python generate_pseudo_phenotype.py --config config.yaml --output-dir /tmp/out --seed 0
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import yaml

# ---------------------------------------------------------------------------
# Distribution registry
# Add new distributions here without touching the rest of the code.
# Signature: (n_samples: int, params: dict, rng: np.random.Generator) -> np.ndarray
# ---------------------------------------------------------------------------
DIST_GENERATORS = {
    "gaussian": lambda n, params, rng: rng.normal(
        loc=params["mean"], scale=params["std"], size=n
    ),
}


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def load_config(config_path: str) -> dict:
    with open(config_path) as fh:
        return yaml.safe_load(fh)


def read_sample_names(
    input_file: str,
    sep: str = "\t",
    skip_first_col: bool = True,
) -> list:
    """Read only the first line of *input_file* and return sample name list.

    Parameters
    ----------
    input_file:
        Path to the expression matrix.
    sep:
        Field delimiter.
    skip_first_col:
        When True (default) the first token on the header line is a feature-ID
        column (e.g. gene IDs) and is dropped — sample names begin at index 1.
    """
    with open(input_file) as fh:
        header = fh.readline().rstrip("\n").split(sep)
    if skip_first_col:
        header = header[1:]
    if not header:
        raise ValueError(f"No sample names found in header of {input_file!r}")
    return header


def write_output(
    rows: list,          # list of (tag: str, values: np.ndarray)
    sample_names: list,
    out_path: Path,
) -> None:
    df = pd.DataFrame(
        {tag: values for tag, values in rows},
        index=sample_names,
    ).T
    df.index.name = "phenotype_id"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t")
    logging.info("  Written: %s  (%d phenotypes × %d samples)", out_path, len(rows), len(sample_names))


# ---------------------------------------------------------------------------
# Core generation
# ---------------------------------------------------------------------------

def generate_phenotype(
    pheno_cfg: dict,
    n_samples: int,
    global_seed: int,
) -> tuple:
    """Draw samples for one phenotype entry.

    Parameters
    ----------
    pheno_cfg:
        Single phenotype block from config (tag, distribution, params, optional seed).
    n_samples:
        Number of samples to draw (= number of columns in the expression matrix).
    global_seed:
        Fallback seed when the phenotype block does not specify its own.

    Returns
    -------
    (tag, values) where values is a 1-D np.ndarray of length n_samples.
    """
    tag = pheno_cfg["tag"]
    distribution = pheno_cfg.get("distribution", "gaussian")
    params = pheno_cfg.get("params", {})
    seed = pheno_cfg.get("seed", global_seed)

    if distribution not in DIST_GENERATORS:
        raise ValueError(
            f"Unknown distribution {distribution!r}. "
            f"Supported: {list(DIST_GENERATORS)}"
        )

    rng = np.random.default_rng(seed)
    values = DIST_GENERATORS[distribution](n_samples, params, rng)
    return tag, values


def process_dataset(
    dataset_cfg: dict,
    global_seed: int,
    output_dir: Path,
    config_dir: Path,
) -> None:
    """Process one dataset entry: read sample names, generate all phenotypes, write TSV.

    Parameters
    ----------
    dataset_cfg:
        One element of cfg["datasets"].
    global_seed:
        Seed resolved from CLI / global config block.
    output_dir:
        Resolved output directory (may override config value).
    config_dir:
        Directory of the config file — used to resolve relative input_file paths.
    """
    raw_path = dataset_cfg["input_file"]
    # Resolve relative paths with respect to the config file location
    input_file = (config_dir / raw_path).resolve()
    dataset_tag = dataset_cfg["dataset_tag"]
    sep = dataset_cfg.get("sep", "\t")
    skip_first_col = dataset_cfg.get("skip_first_col", True)
    phenotypes = dataset_cfg.get("phenotypes", [])

    logging.info("Dataset: %s  →  tag=%s", input_file, dataset_tag)

    sample_names = read_sample_names(str(input_file), sep=sep, skip_first_col=skip_first_col)
    n_samples = len(sample_names)
    logging.info("  Samples found: %d", n_samples)

    rows = []
    for pheno_cfg in phenotypes:
        tag, values = generate_phenotype(pheno_cfg, n_samples, global_seed)
        rows.append((tag, values))
        logging.info("  Phenotype generated: %s", tag)

    out_path = output_dir / f"{dataset_tag}_pseudo_phenotypes.tsv"
    write_output(rows, sample_names, out_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: Optional[list] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate pseudo phenotype data matching sample names of expression matrices.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--config", required=True,
        help="Path to YAML config file.",
    )
    parser.add_argument(
        "--dataset-tag",
        help="Run only the dataset with this dataset_tag (runs all if omitted).",
    )
    parser.add_argument(
        "--output-dir",
        help="Override global.output_dir from config.",
    )
    parser.add_argument(
        "--seed", type=int,
        help="Override global.seed from config.",
    )
    parser.add_argument(
        "--log-level", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity (default: INFO).",
    )
    return parser.parse_args(argv)


def main(argv: Optional[list] = None) -> None:
    args = parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(levelname)s: %(message)s",
    )

    cfg = load_config(args.config)
    config_dir = Path(args.config).resolve().parent

    # Resolve global settings (CLI overrides config)
    global_cfg = cfg.get("global", {})
    output_dir = Path(args.output_dir or global_cfg.get("output_dir", "results/pseudo_phenotype"))
    if not output_dir.is_absolute():
        output_dir = (config_dir / output_dir).resolve()

    global_seed = args.seed if args.seed is not None else global_cfg.get("seed", 42)

    datasets = cfg.get("datasets", [])
    if not datasets:
        logging.error("No datasets defined in config.")
        sys.exit(1)

    # Optional filter by dataset_tag
    if args.dataset_tag:
        datasets = [d for d in datasets if d.get("dataset_tag") == args.dataset_tag]
        if not datasets:
            logging.error("No dataset with dataset_tag=%r found in config.", args.dataset_tag)
            sys.exit(1)

    for dataset_cfg in datasets:
        process_dataset(dataset_cfg, global_seed, output_dir, config_dir)

    logging.info("Done.")


if __name__ == "__main__":
    main()
