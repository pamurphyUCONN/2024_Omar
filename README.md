This repository contains notebooks for analysis and visualization based on the study by Omar et al. (DOI: 10.1101/2023.12.13.571178v1).

1. Demultiplex_add_metadata.ipynb
Sequencing data were processed with Cell Ranger version 7.0.1 on Xanadu high performance computing cluster (UConn Health Center). Reads from FASTQ files were aligned to mouse mm10 or human GRCh38 reference genome as described by 10x genomics. Paths to feature library (antibody capture and gene expression) were specified during the processing, so we ended with an h5 object that included both gene expression and antibody capture information.  h5 Data for processing can be found at https://murphylabvasculardata.cam.uchc.edu/download as the individual batch files (e.g. Batch2022Dec22)
Notebook uses these h5 objects and describes data demultiplexing by hashtags, merging of batches and annotation. Hashtagged nuclei were demultiplexed using DemuxEM using Pegasus implementation in Python with parameters min_signal=15.0, alpha=0.0, alpha_noise=150.0, and deMULTIplex2. Nuclei with ambiguous hashtag assignment were discarded. Metadata was added from sample annotations, and additional information on other assays (e.g. pTau levels from Western blot analysis). Output h5ad file was used as input for inCITE normalization and scrublet removal of duplicates and clustering.

2. InCITE_normalization_and_Scrublet.ipynb
Nuclei were filtered based on gene expression, mitochondrial content, nuclear pore hashtag counts, and oligonucleotide-conjugated antibody levels. Gene counts were normalized and log-scaled. For cell type clustering, 6,830 variable genes were selected, followed by log count scaling and UMI adjustment. Dimensionality reduction was performed using PCA, with batch effects corrected via Harmony. A nearest-neighbor graph was constructed, and clustering was done using the Leiden algorithm (resolution=0.8), followed by UMAP visualization. Suspected doublets were removed using Scrublet, and PCA and Harmony corrections were recalculated. Cell types were annotated using Pegasus and manual labeling, resulting in a final count of 139,342 nuclei.  This h5ad file can be found at https://murphylabvasculardata.cam.uchc.edu/download as the file “All Data”

3. Diffmarkers_bulk_cell_types.ipynb
Marker genes for cell types were visualized in a dot plot. Endothelial cells (ECs) were identified by ERG and FLT1, astrocytes by GFAP and AQP4, neuronal cells by RBFOX3, GABAergic neurons by GAD1 and GAD2, and glutamatergic neurons by SLC17A7. Microglial cells were marked by CX3CR1, ITGAM, and P2RY12. Pericytes expressed PDGFRB and PTN, smooth muscle cells were identified by ACTA2, and oligodendrocytes by MOG, MAG, and PLP1. OPCs were marked by PDGFRA and CSPG4, and macrophages by CD163, MIRC1, and LYVE1. Further analysis subdivided ECs into artery, capillary, vein, and reactive capillary (REV) using specific markers such as FBLN5, ITIH5, CACNA1C, and IFITM2.

4. Correlation_plots.ipynb
Analysis of compositional changes of disease states across age with correlation plots. Scatter plot showing correlation (Spearmann) between donor age and percent contribution in clusters. This was used to generate figures 3G and 4E for microglia and capillaries, respectively. H5ad files for microglia and capillaries can be found at https://murphylabvasculardata.cam.uchc.edu/download as the files “Microglia” and “Capillaries”.

5. Density_plots.ipynb 
Normalized protein data were visualized using density plots based on predefined quintiles of protein expression. Each cell was labeled with its corresponding quintile. The embedding_density function from Scanpy was used to calculate and visualize the density of protein expression within the embedding space.

6. NFkBvTDP_among_DiseaseStates.ipynb
Correlation analysis of NF-kB and TDP43 protein levels across disease states in figures 6 and SI figure 13. The smoothed line plot comparing NF-kB protein levels (x-axis) to TDP-43 (y-axis) across cell types. 

7. Comparative_DEG_and_pathway.ipynb
To explore gene set overlap, we used GSEApy with pre-ranked lists of differentially expressed genes (based on Log2 fold-change) from healthy versus disease capillary clusters, as well as from experimental models involving TDP-43 in human and mouse brain ECs. Differential gene expression was analyzed using DESeq2, with ranked lists filtered for BaseMean expression >5 (provided as supplemental information in parallel manuscript, https://www.biorxiv.org/content/10.1101/2023.12.13.571184v1 ). Gene set enrichment analysis was performed using MSigDB Hallmark and KEGG pathway gene sets, along with custom gene sets derived from our data (provided as SI Table 8 in the manuscript). Overlap in pathways was visualized using NetworkX, and key enriched pathways were plotted in GSEApy.


Processing requirements:
Analysis and clustering of the entire dataset was performed on a high-performance computational cluster with graphics processing unit (GPU) for iterative functions.  Memory requirements are 200GB.  Smaller portions of the data (e.g. Capillaries or Microglia) could be analyzed on most desktop computers.

Software versions:

Python version: 3.10.10 (main, Mar 21 2023, 18:45:11) [GCC 11.2.0]
TensorFlow version: 2.13.0
SciPy version: 1.10.1
Statsmodels version: 0.14.0
NumPy version: 1.24.3
Pandas version: 2.0.2
Scanpy version: 1.9.3
anndata version: 0.8.0
Matplotlib version: 3.7.1
Seaborn version: 0.12.2
Pegasus version: 1.7.1
rpy2 version: 3.5.13
scVI version: 1.0.2
pybiomart version: 0.2.0

Full Conda environment is provided as a YAML file : inCITE_conda_env.yml
