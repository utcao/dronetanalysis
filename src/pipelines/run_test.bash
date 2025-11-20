#!/bin/bash
# Script to test the Snakemake workflow
echo ""
echo "Running Snakemake dry-run..."
snakemake -np --cores 1 -s src/pipelines/Snakemake_template

echo ""
echo "Dry-run completed. Check if all rules are properly scheduled."