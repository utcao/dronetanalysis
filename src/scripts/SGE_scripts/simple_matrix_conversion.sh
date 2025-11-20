#!/bin/sh
# ==============================================================================
# Simple SGE Script for Matrix Conversion
# 
# This shows how to use SGE with command line arguments
# ==============================================================================

# SGE directives (the #$ lines tell SGE what resources you need)
#$ -A dronetanalysis              # Your account/project
#$ -pe parallel 1                # Number of CPU cores (1 is fine for this)
#$ -l h_vmem=16G                  # Memory per core (4GB should be plenty)
#$ -l h_rt=0:2:0                # Time limit
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
# Activate conda environment
# ==============================================================================

echo "Activating conda environment..."

# Option 1: Source the corrected .bashrc to get conda initialization
source ~/.bashrc

# Activate the dronetanalysis environment
conda activate dronetanalysis

# Verify environment activation
echo "Current conda environment: $CONDA_DEFAULT_ENV"
echo "R location: $(which R)"

# ==============================================================================
# Run the R script with command line arguments
# ==============================================================================

# Change to the project directory so relative paths work
cd /tmp/global2/gthornes/dronetanalysis
echo "Changed to directory: $(pwd)"

# Use your existing script with minimal command line args
# Rscript src/scripts/perm_matrix_convert/perm_matrix_convert.R \
#     --input results/spearman_correlation/perm/permutation_test/sig_edges_coexpr_net.csv \
#     --output results/matrices/correlation_matrix_wide.csv \
#     --gene-pairs-col gene_pairs \
#     --value-col rho

echo "Now trying the matrix conversion script..."
timeout 120 Rscript src/scripts/perm_matrix_convert/perm_matrix_convert.R
echo "Matrix conversion completed with exit code: $?"

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