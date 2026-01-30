# conda env packages
# conda activate dronetanalysis
conda_env_path=/tmp/global2/caoyt/miniforge3/envs/dronetanalysis
# conda create -p $conda_env_path r-base=4.3 python=3.12 -y

# Install R packages
micromamba install -p $conda_env_path \
    r-data.table r-reshape2 r-ggplot2 r-dplyr r-tidyr r-readr r-purrr r-stringr \
    r-dofuture r-foreach \
    r-futile.logger r-glue r-yaml r-argparse r-igraph bioconda::bioconductor-limma bioconda::bioconductor-edger bioconda::r-wgcna -y

# Install Python packages
micromamba install -p $conda_env_path \
    numpy pandas scikit-learn matplotlib jupyter biopython statsmodels -y