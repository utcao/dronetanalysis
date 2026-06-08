#!/bin/sh
# ==============================================================================
# submit_gamlss_q1q5.sh — Snakemake cluster-mode master job
#
# This script is a lightweight master job. It runs Snakemake, which submits
# each focal-gene × condition pair as its own qsub child job. The master
# only monitors job status; the heavy computation happens in the child jobs.
#
# Architecture:
#   master job (this script, 1 slot)
#     └─ qsub → gamlss_quintile: FBgn0027560 CT   (4 slots, 25G, 48h)
#     └─ qsub → gamlss_quintile: FBgn0027560 HS   (4 slots, 25G, 48h)
#     └─ qsub → gamlss_quintile: FBgn0001234 CT   (4 slots, 25G, 48h)
#     ...
#
# SUBMIT FROM: code/dronetanalysis/
#   qsub tmp/submit_gamlss_q1q5.sh
#
# DRY RUN (see which jobs would be submitted, no actual submission):
#   snakemake -s tmp/compute_degdvg_gamlss_q1q5.snakemake -n -p
# ==============================================================================

# --- master job resources: small — it only runs Snakemake itself ---
#$ -A gamlss_q1q5
#$ -pe parallel 1
#$ -l h_vmem=4G
#$ -l h_rt=560:0:0    # must stay alive until the last child job finishes

#$ -N gamlss_master
#$ -j y
#$ -m beas
#$ -M yutao.cao@tuebingen.mpg.de

#$ -o logs/$JOB_NAME.$JOB_ID.out
#$ -S /bin/bash
#$ -cwd

set -e
set -x

# --- environment ---
source /tmp/global2/caoyt/miniforge3/etc/profile.d/conda.sh
conda activate ganlss

# --- create directories the child jobs will write into ---
mkdir -p logs
mkdir -p ../../results/gamlss_0608_2026
mkdir -p ../../logs/gamlss_0608_2026

# --- child job qsub template ---
# Snakemake expands {threads} (from the rule's threads: directive) and
# {rule}, {wildcards.gene}, {wildcards.cond} before passing to qsub.
# Each gamlss_quintile rule declares threads: 4, so child jobs get 4 slots.
CHILD_QSUB="qsub \
  -A child_gamlss_q1q5 \
  -pe parallel {threads} \
  -l h_vmem=25G \
  -l h_rt=48:0:0 \
  -N gamlss_{wildcards.gene}_{wildcards.cond} \
  -j y \
  -cwd \
  -o ../../logs/gamlss_0608_2026/{wildcards.gene}_{wildcards.cond}.out"

# --- run Snakemake in cluster mode ---
# -j 20          : max 20 child jobs in the queue at once; raise if you have
#                  many genes and want them all submitted simultaneously
# --latency-wait : wait up to 120 s for output files to appear on NFS after a
#                  child job finishes (avoids false "missing output" errors)
# --keep-going   : if one gene fails, continue with the rest
# --rerun-incomplete : re-run any partially written output from a prior crash
snakemake \
  -s tmp/compute_degdvg_gamlss_q1q5.snakemake \
  -j 50 \
  --cluster "$CHILD_QSUB" \
  --latency-wait 120 \
  --rerun-incomplete \
  --keep-going \
  -p
