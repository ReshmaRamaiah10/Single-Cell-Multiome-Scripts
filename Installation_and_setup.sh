module load gcc/10.2.0
conda create -n multiome_winner -c conda-forge python=3.8 R=4.2
mamba install -c conda-forge get_version
mamba install -c conda-forge rpy2
mamba install -c conda-forge jupyter nb_conda ipykernel nbconvert
mamba install -c conda-forge scanpy python-igraph leidenalg
mamba install -c bioconda anndata2ri
mamba install -c conda-forge r-devtools geos r-xml

pip install fa2
mamba install -c conda-forge pyarrow pytabix pyranges bedtools
 
R
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
devtools::install_github("jokergoo/ComplexHeatmap")
BiocManager::install(c("motifmatchr", "chromVAR", "TFBSTools", "BSgenome", "rtracklayer"))
devtools::install_github("dpeerlab/ArchR", ref="master", repos = BiocManager::repositories())
library(ArchR)
ArchR::installExtraPackages()
 
ipython kernel install --name "multiome_winner" --user
pip install MACS2
pip install PhenoGraph
pip install doubletdetection
mamba install -c anaconda openpyxl
mamba install -c bioconda gseapy
pip install cmake
pip install SEACells
 
R
devtools::install_github("GreenleafLab/chromVARmotifs")
install.packages(tidyr)
BiocManager::install("MAST")
BiocManager::install("fgsea")
