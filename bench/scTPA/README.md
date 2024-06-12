## What is [scTPA](http://sctpa.bio-data.cn/sctpa)
### Introduction
scTPA is a web tool for single-cell transcriptome analysis and annotation based on biological pathway activation in human and mice. We collected a large number of biological pathways with different functional and taxonomic classifications, which facilitates the identification of key pathway signatures for cell type annotation and interpretation.

### What can scTPA do
* Calculating pathway activity score of single cell
* Dimension reduction
* Clustering of cell population by different methods
* Identifying significantly activated pathways of cell clusterings
* Comparison analysis of the associated gene expression profiles of pathways

## Usage
### Install
* **step1 Download scTPA**
scTPA can be download directly  from `Download ZIP` button. Alternatively, scTPA can be installed through github: enter the directory where you would like to install scTPA and run
```
git clone https://github.com/yupenghe/methylpy.git
cd scTPA/
```
* **step2 Install dependent R packages**
For using scTPA, user must install following packages,:
>Seurat
>foreach
>bigstatsr
>data.table
dplyr
scales
ggplot2
cowplot
pheatmap

To install this packages, start "R" and enter:
```
if (!requireNamespace(c("Seurat","bigstatsr","data.table","foreach","dplyr","scales","ggplot2","cowplot","pheatmap"), quietly = TRUE))
    install.packages(c("Seurat","bigstatsr","data.table","foreach","dplyr","scales","ggplot2","cowplot","pheatmap"))
```
* **step3 Install optional R packages**
If user want to use some specialized method in scTPA, the following R packages are required.
1. `scran` for "scran" normalization method
2. `scImpute` for "scImpute" imputation method
3. `SIMLR` for "simlr" clustering method
4. `dbscan` for "dbscan" clustering method 

**scran**
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("scran", quietly = TRUE))
    BiocManager::install("scran")
```
**scImpute**
```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
if (!requireNamespace("scImpute", quietly = TRUE))
    devtools::install_github("Vivianstats/scImpute")
```
**SIMLR**
```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
if (!requireNamespace("scImpute", quietly = TRUE))
    devtools::install_github("Vivianstats/scImpute")
```
**dbscan**
```
if (!requireNamespace("dbscan", quietly = TRUE))
    BiocManager::install("dbscan")
```
### Test scTPA
```
Rscript /path/to/you/scTPA-master/R/scTPA.R -f /path/to/you/scTPA-master/test/expression.csv --cellType /path/to/you/scTPA-master/test/cell_type.csv --work_dir /path/to/you/scTPA-master/ --species homo --pathway_database kegg --para_size 1 -o /path/to/you/scTPA-master/results/
```
Once the program has run successfully, a series of results files and folders will appear in the results folder.
### Command

```
Rscript /path/to/you/scTPA/scTPA.R -h
Options:
    -f FILE, --file=FILE
        gene expression profile, genes X cells
    --cellType=CELLTYPE
        cell type file. First column is cell name (same as the colnames of gene expression profile), second column is cell type. No header names.[default= NULL]
    --work_dir=WORK_DIR
        Workshop direction. [default= ./]
    --normalize=NORMALIZE_METHOD
        methods used for normalization. "log", "CLR", "RC" or "scran"[default= none]
    --min_cells=MIN_CELLS
        genes must be in a minimum number of cells. Used for filtering genes[default= 3]
    --min_features=MIN_FEATURES
        cells must have at least the minimum number of genes. Used for filtering cells[default= 200]
    --species=SPECIES
        species. "homo" or "mus"[default= homo]
    --imputation=IMPUTATION
        Imputation method. "scImpute" or "none"[default= none]
    --data_type=FILE
        data type of gene expression profile，"TPM" or "count"[default= TPM]
    --pathway_database=PATHWAY_DATABASE
        pathway database, detials see https://github.com/sulab-wmu/scTPA[default= kegg]
    --user_pathway=USER_PATHWAY
        user defined pathway file，only for gmt format[default = NULL]
    --pas_method=PAS_METHOD
        method for calculating PAS. "gsva", "ssgsea", "zscore" or "plage"[default= ssgsea]
    --para_size=PARA_SIZE
        number of kernels used for parallel[default= 4]
    --cluster_method=CLUSTER_METHOD
        clustering method. "seurat", "hclust", "simlr", "kmedoids", "kmeans" or "dbscan"[default= seurat]
    --seurat_dims=SEURAT_DIMS
        dimensions used in Seurat clustering[default= 8]
    --seurat_resolution=SEURAT_RESOLUTION
        resolution used for Seurat clustering[default= 0.5]
    --k_cluster=K_CLUSTER
        number of clusters, useless if clustering method is Seurat or dbscan[default= 5]
    --min_pts=MIN_PTS
        parameter in DBSCAN[default= 3]
    --dims=DIMS
        number of PCA dimensions used for TSNE or UMAP[default= 20]
    --marker_method=FIND_MAKER_METHOD
        method of finding siginificant markers[default= wilcox]
    --logFC_thre=THRESHOLD_LOGFC
        threshold of logFC (Detail see Seurat)[default= 0.25]
    --min_pct=MIN_PCT
        only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations.[default= 0.1]
    --pic_type=PIC_TYPE
        type of picture, png or pdf [default= png]
    -o OUT_DIR, --out_dir=OUT_DIR
        output folder[default= NULL]
    -h, --help
        Show this help message and exit
```
#### Details
**`--normalize`:**
***log:*** Log transform. Feature counts output for each cell is divided by the total counts for that cell and multiplied by 1e4. This is then natural-log transformed.
***CLR:*** Centered log ratio. A commonly used Compositional Data Analysis (CoDA) transformation method.
***RC:*** Relative counts. Feature counts output for each cell are is divided by the total counts for that cell and multiplied by 1e4 (for TPM/CPM/FPKM/RPKM this value is 1e6).
***scran:*** The normalization strategy for scRNA-seq is implemented based on the deconvolutional size factor using the scran R package. Detials see [scran](https://github.com/MarioniLab/scran)
***none***: Do not implement normalization

**`--imputation`:**
***scImpute***: Imputing missing value of data matrix following filtering and normalization steps and this function is performed using scImpute R package
***none***: Do not implement imputation.

**`--data_type`:**
***count***: Discrete Data.
***TPM***: Continuous data.

**`--pathway_database`:**

when "--species" is "homo", "--pathway_database" can be select as follow:
***kegg***: An encyclopaedia for genes reaction and regulation. [KEGG](https://www.genome.jp/kegg/). 
***reactome***: A curated database for biomolecular pathways. [Reactome](https://reactome.org/). 
***biocarta***: A pathway database for gene regulation. [BioCarta](https://www.liebertpub.com/doi/pdf/10.1089/152791601750294344). 
***smpdb***: A small molecules pathway database. [SMPDB](https://smpdb.ca/). 
***humancyc***: A curated pathway database of human metabonomics. [HumanCyc](https://humancyc.org/). 
***panther***: A curated pathway database for protein annotation through evolutionary relationship. [PANTHER](http://www.pantherdb.org/). 
***pharmgkb***: A curated pathway database for pharmacogenomics. [pharmGKB](https://www.pharmgkb.org/). 
***acsn2***: A web-based resource depicting signalling and regulatory molecular processes in cancer cell and tumor microenvironment. [ACSN v2.0](https://acsn.curie.fr/ACSN2/ACSN2.html). 
***rb***: A curated map of molecular interactions about retinoblastoma protein (RB/RB1). [RB-Pathways](http://bioinfo-out.curie.fr/projects/rbpathway/). 
***h.all***: Hallmark gene sets. [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp). 
***c2.cgp***: Chemical and genetic perturbations. [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp). 
***c2.cp***: 
***c4.cgn***: Cancer gene neighborhoods. [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp). 
***c4.cm***: Cancer modules. [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp). 
***c5.bp***: GO biological process. [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp). 
***c5.mf***: GO cellular component. [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp). 
***c5.cc***: GO molecular function. [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp). 
***c6.all***: Oncogenic signatures. [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp). 
***c7.all***: Immunologic signatures. [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp). 

when "--species" is mus, "--pathway_database" can be select as follow:
***kegg***: An encyclopaedia for genes reaction and regulation. [KEGG](https://www.genome.jp/kegg/). 
***reactome***: A curated database for biomolecular pathways. [Reactome](https://reactome.org/). 
***smpdb***: A small molecules pathway database. [SMPDB](https://smpdb.ca/). 
***c5.bp***: GO biological process. [GSKB](http://ge-lab.org/gskb/). 
***c5.mf***: GO cellular component. [GSKB](http://ge-lab.org/gskb/). 
***c5.cc***: GO molecular function. [GSKB](http://ge-lab.org/gskb/). 
***other***: Including "Location", "HPO", "STITCH", "MPO", "T3DB", "PID", "MethyCancer" and "MethCancerDB*, details see [table](http://ge-lab.org/gskb/Table%201-sources.pdf). [GSKB](http://ge-lab.org/gskb/). 


