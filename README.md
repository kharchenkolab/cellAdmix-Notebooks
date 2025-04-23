# cellAdmix-Notebooks

To install the dependencies, run these commands from the project root:

```bash
BiocManager::install("sva")
```

```bash
Rscript -e "devtools::install_deps()"
```

To create the necessary data folders, run:

```bash
Rscript -e "dataorganizer::CreateFolders()"
```

## Notebooks

The notebooks are organized by dataset in the [`Analyses`](./Analyses) folder. Within each dataset folder, numbered notebooks should be run first in the order of the numbers. Non-numbered notebooks can be run in any order **after** the numbered ones.

## Data

*The notebooks assume that the data is stored in `./data/spatial_data`. See `R/data_loading.R` for details.*

This section explains how to download and prepare the data required for the notebooks.

To download the data, use the following links:
- [**Human Breast Cancer**](https://www.nature.com/articles/s41467-023-43458-x) (BC)
  - Xenium: [10x Genomics Preview](https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast). Only need 'In Situ Sample 1, Replicate 1', 'Xenium Output Bundle', which should be unzipped into the root dataset folder.
  - scRNA-seq: [GEO GSM7782698](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7782698). Need `GSM7782698_count_raw_feature_bc_matrix.h5` (from GEO) and `umap_annotation.csv` from Source Data -> "Fig. 3E".
- [**Mouse Hypothalamus**](https://doi.org/10.1126/science.aau5324) (Brain)
  - MERFISH: [Dryad](https://doi.org/10.5061/dryad.8t8s248)
  - scRNA-seq: [GEO GSE113576](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113576). Only need `regional_sampling_UMIcounts.txt` and `mouse_int_meta.csv` files.
  - scRNA-seq metadata from [the paper](https://www.science.org/doi/10.1126/science.aau5324) supplementary Table 1 should be saved in `mouse_hypothalamus_rna/scRNA_metadata.xlsx`.
  - There is also `merfish_barcodes.csv` with the molecule information, which wasn't published in the original study.
- **Mouse Ileum** (Gut)
  - MERFISH: [Dryad](https://doi.org/10.5061/dryad.jm63xsjb2)
  - scRNA-seq: [GEO GSE92332](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92332)
- **Human NSCLC** (NSCLC)
  - CosMx: [Nanostring FFPE](https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/).
  - scRNA-seq: [GEO GSE127465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127465)
  - DIALOGUE results: [GitHub](https://github.com/livnatje/DIALOGUE/wiki/Reproducing-the-lung-cancer-spatial-results-and-figures). Put `DIALOGUE1_LungCancer.SMI.rds` in `NSCLC/CosmX`.
  - iTALK LR database: [GitHub](https://github.com/Coolgenome/iTALK/tree/master/data). Put `LR_database.rda` in `NSCLC/CosmX`.
- **Human Ovarian Cancer** (OC)
  - Xenium: [10x Genomics](https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-ovarian-cancer)
  - scRNA-seq: [10x Genomics scFFPE](https://www.10xgenomics.com/datasets/17k-human-ovarian-cancer-scFFPE)
- **Human Pancreatic Cancer** (Pancreas)
  - Xenium: [10x Genomics](https://www.10xgenomics.com/datasets/ffpe-human-pancreas-with-xenium-multimodal-cell-segmentation-1-standard)
  - snRNA-seq: [Single Cell Portal](http://singlecell.charite.de/cellbrowser/pancreas/)

The code assumes that the data is stored in `../data` with a separate subfolder per dataset and per modality. See [data_mapping.yml](./data_mapping.yml) for details. There, you can also change the paths to match your local setup.


## Individual dataset notes

### Mouse Ileum

Run [GUT_preprocessing.ipynb](./Analyses/Gut_analyses/GUT_preprocessing.ipynb) to preprocess the data.