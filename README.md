# PathSingle: A Biochemical Pathway Analysis Tool for Single-Cell Data
![Pathway analysis](code/data/Pathsingle.png)

## Overview
PathSingle is a Python-based pathway analysis tool tailored for single-cell data analysis. It employs a unique graph-based algorithm to enable the analysis of diverse cellular states, such as T cell subtypes. Designed to be open-source, extensible, and computationally efficient, PathSingle provides researchers with a versatile framework for uncovering biologically meaningful insights from high-dimensional single-cell transcriptomics data.

## Key Features
- Tailored for single-cell RNA-seq data analysis
- Unique graph-based algorithm for pathway analysis
- Efficient classification of diverse cellular states
- Open-source and extensible
- Computationally efficient

## Installation
`pip install pathsingle`

Or simply clone this repository using git clone command.

## Quick Start
```
anndata = sc.read('./data/activity_df.csv', delimiter=',', cache=False)
calc_activity(anndata)
```

## Usage
For detailed usage instructions, please refer to the scripts and notebooks in the `pathsingle` folder.

[Single Cells Benchmark Readme](code/benchmark.md)

[Feature Reduction Notebook](code/feature_reduction.ipynb)

The [pathway_db](code/pathway_db) subfolder contains scripts for building and customizing the pathways database.
- [kegg](code/pathway_db/kegg): Scripts for downloading pathways from KEGG PATHWAY database and building the existing PathSingle database pathway_relations.csv.
- [reactome](code/pathway_db/reactome): Scripts for downloading and building a new PathSingle pathway_relations.csv database using Reactome Pathway database.
- [pathwaycommon](code/pathway_db/pathwaycommon): Scripts for downloading and building a new PathSingle pathway_relations.csv database using Pathway Commons integrated database.

## Supported Pathways
PathSingle currently supports 357 curated pathways. Click the link to view the full list. [List of supported pathways.](code/data/pathway_relations.csv)

## Contributing
We welcome contributions! Please see our Contributing Guidelines for more information on how to get involved.

## License
PathSingle is available under the MIT license. See the LICENSE file for more details.

## Support
For questions, issues, or feature requests, please open an issue on our GitHub repository.

For additional support, contact: zurkin at yahoo dot com.

## Citation
If you use PathSingle in your research, please cite our paper:
Livne, D., Efroni, S. Pathway metrics accurately stratify T cells to their cells states. BioData Mining 17, 60 (2024). https://doi.org/10.1186/s13040-024-00416-7

## Acknowledgments
We thank the reviewers of this work for their valuable feedback and contributions.