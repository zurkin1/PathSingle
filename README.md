# PathSingle A Biochemical Pathway Analysis Tool
![Pathway analysis](https://norbis.w.uib.no/files/2016/05/F1.large_-768x623.jpg)

PathSingle is a gene regulatory network analysis tool designed specifically for single-cell data. It can handle the substantial data volumes associated with single-cell experiments, as well as addressing sparse and potentially erroneous data points. As we developed ![PathWeigh](https://github.com/zurkin1/Pathweigh), we recognized the unique challenges posed by single-cell RNA sequencing (scRNA-seq) data, such as dropouts and large sparse matrices. These challenges necessitated a distinct approach and development process, prompting the creation of PathSingle for single cell, a specialized pathway analysis tool designed to handle scRNA-seq data.
Based on Python and integrated with the popular Scanpy package [1], PathSingle accepts a gene expression matrix of cells and produces the activity level of 581 pathways within these cells.

## Running
In order to run PathSingle please refer to the notebook folder where the notebook is availabe, demonstrating single-cell RNAseq data.

### [List of supported pathways.](data/pathnames.txt)

### [Guide for adding a new pathway.](data/guide.md)

Support: zurkin at gmail dot com.

[1] Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. Genome Biol. 2018;19:15.