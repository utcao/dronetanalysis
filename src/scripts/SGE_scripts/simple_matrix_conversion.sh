#!/bin/sh
# ==============================================================================
# Simple SGE Script for Matrix Conversion
# 
# This shows how to use SGE with command line arguments
# ==============================================================================

# SGE directives (the #$ lines tell SGE what resources you need)
#$ -A dronetanalysis              # Your account/project
#$ -pe parallel 1                # Number of CPU cores (1 is fine for this)
#$ -l h_vmem=4G                  # Memory per core (4GB should be plenty)
#$ -l h_rt=1:0:0                 # Time limit (1 hour)
#$ -q all.q                      # Queue to use
#$ -N matrix_convert              # Job name (shows up in qstat)
#$ -j y                          # Merge stdout and stderr into one file
#$ -m eas                        # Email when job ends, aborts, or starts
#$ -M gabriel.thornes@tuebingen.mpg.de # Your email address
#$ -o logs/matrix_convert_$JOB_ID.out  # Output file
#$ -S /bin/bash                  # Use bash shell
#$ -cwd                          # Run from current working directory

# ==============================================================================
# Print job information for debugging
# ==============================================================================

echo "=== SGE Job Information ==="
echo "Job ID: $JOB_ID"
echo "Job name: $JOB_NAME"
echo "Node: $(hostname)"
echo "Working directory: $(pwd)"
echo "Started at: $(date)"
echo "=========================="

# ==============================================================================
# Load required modules (modify based on your cluster)
# ==============================================================================

# Uncomment and modify these lines based on your cluster setup:
# module load R/4.3.0
# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate bioinformatics_env

# ==============================================================================
# Run the R script with command line arguments
# ==============================================================================

# Basic example - just specify input file
# Rscript src/scripts/examples/simple_matrix_conversion.R \
#     --input results/spearman_correlation/permutation_test/sig_edges_coexpr_net.csv

# More complete example with all options
Rscript src/scripts/examples/simple_matrix_conversion.R \
    --input results/spearman_correlation/permutation_test/sig_edges_coexpr_net.csv \
    --output-dir results/matrices \
    --gene-pairs-col gene_pairs \
    --value-col rho \
    --verbose

# ==============================================================================
# Print completion information
# ==============================================================================

echo "=========================="
echo "Job completed at: $(date)"
echo "Exit status: $?"
echo "=========================="

# ==============================================================================
# HOW TO USE THIS SGE SCRIPT:
# ==============================================================================
#
# 1. Save this as a file (e.g., run_matrix_conversion.sh)
# 2. Make it executable: chmod +x run_matrix_conversion.sh
# 3. Submit to SGE: qsub run_matrix_conversion.sh
# 4. Check status: qstat
# 5. View output: cat logs/matrix_convert_JOBID.out
#
# To modify for different input files, just change the --input argument above
# ==============================================================================