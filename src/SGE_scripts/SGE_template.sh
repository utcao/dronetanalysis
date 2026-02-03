#!/bin/sh
# Add another # at the beginning of the line to comment out a line
# NOTE: adding a switch to the command line will override values in this file.
# These options are MANDATORY in LCRC; Your qsub will fail if you don't provide them.
# Adjust the number of node, processors, memory and maximum mem as needed!
#$ -A randeff_remove
#$ -pe parallel 1
#$ -l h_vmem=50G
#$ -l h_rt=10:0:0
# Job queue specification (optional but recommended)
#$ -q all.q
# Highly recommended
# The first 15 characters of the job name are displayed in the qstat output:
#$ -N randeff_remove_fam
# Job arrays support (uncomment and modify as needed)
# #$ -t 1-1000
# If you want to merge stdout and stderr, use the -j option
# y=yes, n=don't merge
#$ -j y
# Controlling email notifications
# When to send email b=job begin, e=job end, a=job abort, j=subjobs (job arrays), n=no mail
#$ -m beas
#$ -M yutao.cao@tuebingen.mpg.de
#$ -o logs/$JOB_NAME.$JOB_ID.out
#$ -e logs/$JOB_NAME.$JOB_ID.err
# Use /bin/bash to execute this script
#$ -S /bin/bash
# Run job from current working directory
#$ -cwd
# Set up environment
# Load system base, then conda for specific packages)
# module load R/4.3.0
# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate bioinformatics_env

# Print job information for debugging
echo "Job started at: $(date)"
echo "Job ID: $JOB_ID"
echo "Job name: $JOB_NAME"
echo "Node: $(hostname)"
echo "Working directory: $(pwd)"
if [ ! -z "$SGE_TASK_ID" ]; then
    echo "Array task ID: $SGE_TASK_ID"
fi
echo "=========================================="

# Your commands go here
# Example: Rscript your_script.R

# Print completion info
echo "=========================================="
echo "Job completed at: $(date)"