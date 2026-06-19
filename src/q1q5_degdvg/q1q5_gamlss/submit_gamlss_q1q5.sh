#!/bin/sh
# ==============================================================================
# submit_gamlss_q1q5.sh — qsub master job for the GAMLSS Q1-vs-Q5 workflow
#
# Submits this script as a lightweight master job. The master runs Snakemake,
# which in turn submits one qsub child job per focal gene x condition pair using
# the cluster-generic executor plugin (Snakemake >=8 requirement).
#
# Architecture:
#   master job (this script, 1 slot, base conda env)
#     |-- child qsub: FBgn0027560 CT  (4 slots, 25 GB, 48 h, ganlss Rscript)
#     |-- child qsub: FBgn0027560 HS  (4 slots, 25 GB, 48 h, ganlss Rscript)
#     |-- child qsub: FBgn0001234 CT  ...
#     ...
#
# Child job resources and qsub flags are defined in:
#   src/q1q5_degdvg/q1q5_gamlss/profiles/sge/config.yaml
#
# Focal gene list: src/q1q5_degdvg/shared/focal_genes.csv  (gene_id + annotation columns)
#
# SUBMIT FROM: code/dronetanalysis/
#   qsub src/q1q5_degdvg/q1q5_gamlss/submit_gamlss_q1q5.sh
#
# DRY RUN (interactive, no jobs submitted):
#   /tmp/global2/caoyt/miniforge3/bin/snakemake \
#     -s src/q1q5_degdvg/q1q5_gamlss/compute_degdvg_gamlss_q1q5.snakemake \
#     --profile src/q1q5_degdvg/q1q5_gamlss/profiles/sge -n -p
# ==============================================================================

# --- master job resources ---
# The master only runs Snakemake; all heavy work is in child jobs.
# h_rt must exceed the longest expected child job (48 h) plus queue wait time.
#$ -A gamlss_q1q5
#$ -pe parallel 1
#$ -l h_vmem=4G
#$ -l h_rt=560:0:0

#$ -N gamlss_master
#$ -j y
#$ -m beas
#$ -M yutao.cao@tuebingen.mpg.de

#$ -o ../../logs/gamlss_0609/$JOB_NAME.$JOB_ID.out
#$ -S /bin/bash
#$ -cwd

set -e
set -x

# Snakemake runs from the base conda env, which has the cluster-generic executor
# plugin (snakemake-executor-plugin-cluster-generic) installed.
# The ganlss env is NOT activated here; child jobs call Rscript via its full
# path (/tmp/global2/caoyt/miniforge3/envs/ganlss/bin/Rscript).
SNAKEMAKE=/tmp/global2/caoyt/miniforge3/bin/snakemake

# --- create output directories before any jobs start ---
mkdir -p ../../results/gamlss_0609
mkdir -p ../../logs/gamlss_0609

# --- launch Snakemake in cluster mode ---
# All cluster settings (qsub flags, job count, latency-wait, etc.) are in
# the profile config; only the snakemake file path is given here.
$SNAKEMAKE \
  -s src/q1q5_degdvg/q1q5_gamlss/compute_degdvg_gamlss_q1q5.snakemake \
  --profile src/q1q5_degdvg/q1q5_gamlss/profiles/sge
