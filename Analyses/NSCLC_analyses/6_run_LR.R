
suppressPackageStartupMessages({
  library(cellAdmix)
  library(Seurat)
  library(Giotto)
  library(data.table)
  library(dplyr)
  library(glue)
  library(Matrix)
  library(readr)
  library(dataorganizer)
  library(devtools)

  load_all()
})

suppressMessages({
  load_arial_font()
})

options(future.globals.maxSize = 4000 * 1024^2)

base_dir <- CachePath('nsclc_scaled_dat11/')
so_spatial_orig <- file.path(base_dir, 'so_spatial_orig.rds') %>% readRDS()
so_spatial_cln_fib <- file.path(base_dir, 'so_spatial_cln_15_fib_full.rds') %>% readRDS()
so_spatial_cln_macro <- file.path(base_dir, 'so_spatial_cln_15_macro_full.rds') %>% readRDS()

message('Data loaded')

# loading the already saved, cleaned fibroblasts data
dat <- so_spatial_orig[['RNA']]$counts
meta <- so_spatial_orig@meta.data

dat_cln <- so_spatial_cln_macro[['RNA']]$counts
fib_dat_cln <- so_spatial_cln_fib[['RNA']]$counts

# creating cleaned copy of the original data with the cleaned counts in the right places
dat_cln_full <- dat
cells_replace_macro <- colnames(dat_cln)
cells_replace_fib <- colnames(fib_dat_cln)
dat_cln_full[,cells_replace_macro] <- dat_cln
dat_cln_full[,cells_replace_fib] <- fib_dat_cln

## adding coarse cell type annotation to reduce number of pairwise ct tests run
immune_cell_types <- c(
  'B-cell','NK','T CD4 memory','T CD4 naive','T CD8 memory','T CD8 naive',
  'Treg','plasmablast','mast','mDC','monocyte','pDC','neutrophil'
)
meta$cell_type_coarse <- ifelse(
  meta$cell_type %in% immune_cell_types,
  'immune cell other', meta$cell_type
)

message('Metadata prepared')

# appending region labels to the metadata
cell_annot <- prepare_nsclc_metadata(reps='all')
meta <- cbind.data.frame(meta,cell_annot[rownames(meta),c('niche','Run_Tissue_name','x','y')])

gem_sub_orig <- createGiottoObject(
  dat,
  spatial_locs = meta[,c('x','y','cell')],
  cell_metadata = meta
)

gem_sub_orig <- normalizeGiotto(gobject = gem_sub_orig)

# create network
gem_sub_orig <- createSpatialNetwork(gobject = gem_sub_orig, minimum_k = 2)

message('Spatial network orig created')

gem_sub_cln <- createGiottoObject(
  dat_cln_full,
  spatial_locs = meta[,c('x','y','cell')],
  cell_metadata = meta
)

gem_sub_cln <- normalizeGiotto(gobject = gem_sub_cln)

# create network
gem_sub_cln <- createSpatialNetwork(gobject = gem_sub_cln, minimum_k = 2)

message('Spatial network cln created')

# using iTALK database
load(DatasetPath('nsclc', 'LR_database.rda'))

lr_pairs <- database[,c('Ligand.ApprovedSymbol','Receptor.ApprovedSymbol')]
colnames(lr_pairs) <- c('ligand','receptor')
lr_pairs <- unique(lr_pairs)

lr_pairs <- lr_pairs[lr_pairs$ligand %in% rownames(dat),]
lr_pairs <- lr_pairs[lr_pairs$receptor %in% rownames(dat),]

message('LR pairs prepared')
## get statistical significance of gene pair expression changes upon cell-cell interaction
spatial_all_scores_orig = spatCellCellcom(gem_sub_orig,
                                          spatial_network_name = 'Delaunay_network',
                                          cluster_column = 'cell_type_coarse',
                                          random_iter = 1000,
                                          feat_set_1 = lr_pairs$ligand,
                                          feat_set_2 = lr_pairs$receptor,
                                          adjust_method = 'fdr',
                                          do_parallel = T,
                                          cores = 30,
                                          verbose = 'a lot')

write_rds(spatial_all_scores_orig, file.path(base_dir, 'LR_scores_full_orig.rds'))

message('Spatial scores orig created')

spatial_all_scores_cln = spatCellCellcom(gem_sub_cln,
                                         spatial_network_name = 'Delaunay_network',
                                         cluster_column = 'cell_type_coarse',
                                         random_iter = 1000,
                                         feat_set_1 = lr_pairs$ligand,
                                         feat_set_2 = lr_pairs$receptor,
                                         adjust_method = 'fdr',
                                         do_parallel = T,
                                         cores = 30,
                                         verbose = 'a lot')

write_rds(spatial_all_scores_cln, file.path(base_dir, 'LR_scores_full_cln.rds'))

message('Spatial scores cln created')
message('All done!')