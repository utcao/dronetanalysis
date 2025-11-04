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
#$ -l h_rt=0:1:0                # Time limit
#$ -q standard.q                 # Queue to use (available: standard.q, long.q, test.q, cryo-em.q)
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

# Create logs directory if it doesn't exist
mkdir -p logs

# ==============================================================================
# Load required modules (modify based on your cluster)
# ==============================================================================

# Environment should be inherited from submission environment
echo "Using inherited conda environment"
echo "R location: $(which R)"
echo "R version: $(R --version | head -1)"

# ==============================================================================
# Run the R script with command line arguments
# ==============================================================================

# Change to the project directory so relative paths work
cd /tmp/global2/gthornes/dronetanalysis
echo "Changed to directory: $(pwd)"
echo "About to run R script..."

# Use your existing script with minimal command line args
# Rscript src/scripts/perm_matrix_convert/perm_matrix_convert.R \
#     --input results/spearman_correlation/perm/permutation_test/sig_edges_coexpr_net.csv \
#     --output results/matrices/correlation_matrix_wide.csv \
#     --gene-pairs-col gene_pairs \
#     --value-col rho

# Alternative: use defaults (no arguments needed)
echo "Starting Rscript with full conda path..."
/tmp/global2/gthornes/miniforge3/envs/dronetanalysis/bin/Rscript src/scripts/perm_matrix_convert/perm_matrix_convert.R

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