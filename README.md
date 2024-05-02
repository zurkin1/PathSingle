# PathSingle - A Biochemical Pathway Analysis Tool
![PathSingle](data/TLR%20Signaling%20Cascade%20Pathway.png)
Created with BioRender.com

## Single Cell Data
Single cell RNA sequencing allows researchers to quantify the gene expression profiles of individual cells. Using multiple approaches, cells are isolated, their RNA content is captured, reversed transcribed into cDNA, amplified to create libraries of cDNA fragments, and then sequenced. There are still several challenges with this technique. For example, low RNA capture efficiency from a single cell. This results in many zero values in the data, also called “dropouts”, and a high level of technical noise. Technical variability can make it difficult to distinguish true biological variation from technical noise.
To address this issue of missing data, multiple imputation tools have been developed and are now part of the standard pipeline for single-cell analysis. These tools perform data imputation to estimate missing gene expression values. After evaluating several options, we selected MAGIC [1] for imputation in our analysis.

Single-cell RNA sequencing, as opposed to bulk RNA sequencing, generates large and sparse data matrices. Processing these matrices on a standard machine is extremely time-consuming. We therefore implemented algorithmic improvements to enhance performance. Specifically, PathSingle incorporates two key optimizations. First, the pathway topology is encoded into an adjacency matrix to enable rapid access during propagation. Second, calculations are parallelized across available CPU cores to exploit concurrent processing. The sparsity and scale of single-cell expression matrices poses computational challenges for efficient analysis. By caching topological information and distributing calculations across CPUs, PathSingle achieves substantial performance gains over naive implementations. These optimizations are critical to facilitating rapid whole-transcriptome pathway analysis on conventional hardware. Further refinement of algorithms and data structures tailored to single-cell omics will help unlock deeper biological insights from these emerging high-dimensional datasets.

As we developed [PathWeigh](https://github.com/zurkin1/Pathweigh), we recognized the unique challenges posed by single-cell RNA sequencing (scRNA-seq) data, such as dropouts and large sparse matrices. These challenges necessitated a distinct approach and development process, prompting the creation of PathSingle for single cell, a specialized pathway analysis tool designed to handle scRNA-seq data.
Based on Python and integrated with the popular Scanpy package [1], PathSingle accepts a gene expression matrix of cells and produces the activity level of 581 pathways within these cells.

## Running
In order to run PathSingle please refer to the notebook folder where the notebook is availabe, demonstrating single-cell RNAseq data.


### [List of supported pathways.](data/pathnames.txt)

### [Guide for adding a new pathway.](data/guide.md)

Support: zurkin at gmail dot com.

[1] Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. Genome Biol. 2018;19:15.