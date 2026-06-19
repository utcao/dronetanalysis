#!/bin/sh
# ==============================================================================
# submit_limma_voom_q1q5.sh — qsub master job for the limma-voom Q1-vs-Q5 workflow
#
# Submits this script as a lightweight master job. The master runs Snakemake,
# which in turn submits one qsub child job per focal gene x condition pair using
# the cluster-generic executor plugin (Snakemake >=8 requirement).
#
# Architecture:
#   master job (this script, 1 slot, base conda env)
#     |-- child qsub: FBgn0027560 CT  (1 slot, 10 GB, 2 h, ganlss Rscript)
#     |-- child qsub: FBgn0027560 HS  (1 slot, 10 GB, 2 h, ganlss Rscript)
#     |-- child qsub: FBgn0001234 CT  ...
#     ...
#
# Child job resources and qsub flags are defined in:
#   src/q1q5_degdvg/q1q5_limma_voom/profiles/sge/config.yaml
#
# Focal gene list: src/q1q5_degdvg/shared/focal_genes.csv  (gene_id + annotation columns)
#
# SUBMIT FROM: code/dronetanalysis/
#   qsub src/q1q5_degdvg/q1q5_limma_voom/submit_limma_voom_q1q5.sh
#
# DRY RUN (interactive, no jobs submitted):
#   /tmp/global2/caoyt/miniforge3/bin/snakemake \
#     -s src/q1q5_degdvg/q1q5_limma_voom/compute_deg_limma_voom_q1q5.snakemake \
#     --profile src/q1q5_degdvg/q1q5_limma_voom/profiles/sge -n -p
# ==============================================================================

# --- master job resources ---
# The master only runs Snakemake; all heavy work is in child jobs.
# h_rt must exceed the longest expected child job (2 h) plus queue wait time.
#$ -A limma_voom_q1q5
#$ -pe parallel 1
#$ -l h_vmem=4G
#$ -l h_rt=24:0:0

#$ -N limma_voom_master
#$ -j y
#$ -m beas
#$ -M yutao.cao@tuebingen.mpg.de

#$ -o ../../logs/limma_voom/$JOB_NAME.$JOB_ID.out
#$ -S /bin/bash
#$ -cwd

set -e
set -x

SNAKEMAKE=/tmp/global2/caoyt/miniforge3/bin/snakemake

mkdir -p ../../results/limma_voom
mkdir -p ../../logs/limma_voom

$SNAKEMAKE \
  -s src/q1q5_degdvg/q1q5_limma_voom/compute_deg_limma_voom_q1q5.snakemake \
  --profile src/q1q5_degdvg/q1q5_limma_voom/profiles/sge
