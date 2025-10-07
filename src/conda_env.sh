conda install -p /tmp/global2/gthornes/miniforge3/envs/dronetanalysis/ r-base python -y

# Install R packages
conda install -p /tmp/global2/gthornes/miniforge3/envs/dronetanalysis/ \
    r-data.table r-ggplot2 r-dplyr r-tidyr r-readr r-purrr r-stringr \
    r-reshape2 r-glue r-yaml -y

# Install DESeq2 using R BiocManager (more reliable than conda)
# After running this script, start R and run:
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("DESeq2")

# Install Python packages
conda install -p /tmp/global2/gthornes/miniforge3/envs/dronetanalysis/ \
    numpy pandas scikit-learn matplotlib jupyter biopython -y