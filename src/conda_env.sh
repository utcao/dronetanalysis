conda install -p /tmp/global2/gthornes/miniforge3/envs/dronetanalysis/ r-base python -y

# Install R packages
conda install -p /tmp/global2/gthornes/miniforge3/envs/dronetanalysis/ \
    r-data.table r-ggplot2 r-dplyr r-tidyr r-readr r-purrr r-stringr \
    r-reshape2 r-glue -y

# Install Python packages
conda install -p /tmp/global2/gthornes/miniforge3/envs/dronetanalysis/ \
    numpy pandas scikit-learn matplotlib jupyter biopython -y