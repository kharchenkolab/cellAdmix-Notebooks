# Requires [OC_benchmark_segmentations.ipynb](./OC_benchmark_segmentations.ipynb) and script 5.

suppressPackageStartupMessages({
  library(cellAdmix)
  library(cowplot)
  library(dataorganizer)
  library(ggplot2)
  library(magrittr)
  library(Matrix)
  library(readr)
  library(Seurat)
  library(gtools)
  library(dplyr)

  devtools::load_all()
})

theme_set(theme_bw())

# Load cell type annotations ---------------------------------------------
cell_types <- CachePath('oc_seg_annot_transfer5.rds') %>% readRDS()
cell_types_methods <- sapply(names(cell_types), function(x) {
  strsplit(x, split='_')[[1]][[1]]
})


# Functions

#' Load segmentation data for a specific method
#'
#' @param method Character string: "xenium", "baysor", "bidcell", or "proseg"
#' @param min_count Minimum number of molecules per cell to keep
#' @return List containing dataframe of molecules and cell annotations
load_df_seg <- function(method, min_count=10) {
  if (method == 'xenium') {
    data <- prepare_OC_sc_spatial()
    df <- data[[1]]
    cell_annot <- data[[2]]
    return(list(df, cell_annot))
  }

  if (method == 'baysor') {
    df <- DatasetPath(
      'human_ovarian_cancer', 'seg_method_results', 'baysor', 'segmentation.parquet'
    ) %>%  arrow::read_parquet()
    df <- as.data.frame(df)
    df <- df[df$cell != "", ]
    df <- df[!df$is_noise, ]
    ct_annot <- cell_types[which(cell_types_methods == 'Baysor')]
    ct_annot_bcodes <- sapply(names(ct_annot), function(x) {
      strsplit(x, split='_')[[1]][[5]]
    })
    # df$cell <- paste0('cell', df$cell)
    ndx_match <- match(df$cell, ct_annot_bcodes)
    df <- df[!is.na(ndx_match), ]
    ndx_match <- match(df$cell, ct_annot_bcodes)
    df$celltype <- ct_annot[ndx_match]

  } else if (method == 'bidcell') {
    df <- DatasetPath(
      'human_ovarian_cancer', 'seg_method_results', 'bidcell',
      'transcripts_processed.csv'
    ) %>% read_csv()
    df <- as.data.frame(df)
    bidcell_ids <- CachePath('oc_df_small_bidcell_ids.csv') %>% read_csv()
    bidcell_ids <- as.data.frame(bidcell_ids)
    df$cell <- bidcell_ids$cell_id
    df <- df[df$cell != 0, ]

    ct_annot <- cell_types[which(cell_types_methods == 'BIDCell')]
    ct_annot_bcodes <- sapply(names(ct_annot), function(x) {
      strsplit(x, split='_')[[1]][[2]]
    })
    ndx_match <- match(df$cell, ct_annot_bcodes)
    df <- df[!is.na(ndx_match), ]
    ndx_match <- match(df$cell, ct_annot_bcodes)
    df$celltype <- ct_annot[ndx_match]

    # Modify column names
    ndx_change <- which(colnames(df) == 'feature_name')
    colnames(df)[ndx_change] <- 'gene'

    ndx_change <- which(colnames(df) == 'x_location')
    colnames(df)[ndx_change] <- 'x'
    ndx_change <- which(colnames(df) == 'y_location')
    colnames(df)[ndx_change] <- 'y'
    ndx_change <- which(colnames(df) == 'z_location')
    colnames(df)[ndx_change] <- 'z'

    df$cell <- paste0('cell', as.character(df$cell))

  } else if (method == 'proseg') {
    df <- DatasetPath(
      'human_ovarian_cancer', 'seg_method_results', 'proseg',
      'transcript-metadata.csv.gz'
    ) %>% read_csv()
    df <- as.data.frame(df)
    ndx_change <- which(colnames(df) == 'assignment')
    colnames(df)[ndx_change] <- 'cell'
    df <- df[df$cell != 4294967295, ]

    ct_annot <- cell_types[which(cell_types_methods == 'ProSeg')]
    ct_annot_bcodes <- sapply(names(ct_annot), function(x) {
      strsplit(x, split='_')[[1]][[2]]
    })
    ndx_match <- match(df$cell, ct_annot_bcodes)
    df <- df[!is.na(ndx_match), ]
    ndx_match <- match(df$cell, ct_annot_bcodes)
    df$celltype <- ct_annot[ndx_match]

    df$cell <- as.character(df$cell)
    df$cell <- paste0('cell', df$cell)
  }

  # Create cell annotation dataframe
  cell_annot <- df %>%
    group_by(cell) %>%
    summarise(
      x = mean(x, na.rm = TRUE),
      y = mean(y, na.rm = TRUE),
      celltype = unique(celltype)
    )
  cell_annot <- as.data.frame(cell_annot)
  cell_annot$cell <- as.character(cell_annot$cell)
  rownames(cell_annot) <- cell_annot$cell

  # Remove cells with very few molecules
  cell_counts <- table(df$cell)
  cells_keep <- names(cell_counts)[cell_counts >= min_count]
  df <- df[df$cell %in% cells_keep, ]
  cell_annot <- cell_annot[rownames(cell_annot) %in% cells_keep, ]

  df$mol_id <- as.character(1:nrow(df))
  rownames(df) <- df$mol_id

  return(list(df, cell_annot))
}

# Execute analysis pipeline ----------------------------------------------

base_dir <- CachePath('OC_scaled_dat6/')
methods_test <- c('baysor', 'bidcell', 'proseg')  # xenium already done
k <- 20

# Run NMF-CRF analysis
run_nmf_crf_all_seg(
  methods_test,
  k=k,
  h_nmf=5,
  h_crf=10,
  base_dir,
  num_cells_samp=2000,
  nmol_dsamp=10000,
  top_g_thresh=.25,
  n_hvgs=NULL,
  mols_test=NULL,
  proj_h=50,
  bridge_ncells_samp=500,
  bridge_knn_k=20,
  n.cores=10
)

# Clean data and save cleaned versions -----------------------------------

for (method in methods_test) {
  print(method)
  df_dat <- load_df_seg(method)
  df <- df_dat[[1]]
  cell_annot <- df_dat[[2]]

  all_ctypes <- unique(cell_annot$celltype)

  # Load bridge results
  if (method == 'xenium') {
    bridge_raw <- readRDS(file=paste0(base_dir, 'bridge_joint_20.rds'))
  } else {
    bridge_raw <- readRDS(file=paste0(base_dir, 'bridge_', method, '.rds'))
  }

  # Parse bridge results
  bridge_res <- extract_bridge_res(
    bridge_raw, all_ctypes, k, nmf.type='joint', p.thresh=.1, adj.pvals=FALSE
  )

  # Load CRF results
  if (method == 'xenium') {
    crf_res <- readRDS(file=paste0(base_dir, 'crf_joint_k', k, '.rds'))
  } else {
    crf_res <- readRDS(file=paste0(base_dir, 'crf_joint_k', k, '_', method, '.rds'))
  }

  # Set up for FP checks
  cell_annot$z <- 1
  crf_all <- list(crf_res)
  names(crf_all) <- method
  bridge_res_all <- list(bridge_res)
  names(bridge_res_all) <- method

  # Run false positive checks
  fp_checks <- check_fp_multi(
    df, cell_annot, crf_all,
    enr_res_all=NULL,
    bridge_res_all=bridge_res_all,
    memb_res_all=NULL,
    do_clean=TRUE,
    knn_k=100,
    median_thresh=.2
  )

  annot_res_all <- fp_checks[[1]]
  orig_nms_all <- fp_checks[[2]]

  # Clean the data and save
  df$factor <- crf_res[, 1]
  df_cln <- df[!(paste(df$factor, df$celltype, sep="_") %in% annot_res_all[[1]]), ]
  saveRDS(df_cln, file=paste0(base_dir, 'df_k', k, '_cln_', method, '.rds'))
}
