suppressPackageStartupMessages({
  library(dataorganizer)
  library(devtools)
  library(cellAdmix)
  devtools::load_all()
})

suppressMessages({
  load_arial_font()
})

# load data
data <- prepare_brain_tx_and_meta()
cell_annot <- data[[1]]
df <- data[[2]]
markers <- load_markers_brain()

base_dir <- CachePath('brain_scaled_dat4/')

k_joint <- c(5,15,20,30)
k_ct <- c(5)

gc()

# load the precomputed NMF and CRF results (from NSCLC_nmf_crf_benchmark.R)
nmf_crf <- load_nmf_crf_multi(k_joint, k_ct, dir=base_dir)
all_nmf <- nmf_crf[[1]]
crf_all <- nmf_crf[[2]]

### run bridge tests first

# get dataset-wide molecule knn graph to identify neighboring cells and molecules
all_id <- unique(cell_annot$Animal_ID)
all_breg <- unique(cell_annot$Bregma)
all_cells <- unique(cell_annot$cell)
all_ct <- unique(cell_annot$celltype)
ncm_cell <- Matrix(0,nrow = nrow(df),ncol=length(all_cells),sparse=TRUE)
rownames(ncm_cell) <- rownames(df)
colnames(ncm_cell) <- all_cells
all_ncm_ct <- list()
all_ncm_cell <- list()
for (id in all_id) {
  print(id)
  for (breg in all_breg) {
    df_sub <- df[df$Bregma==breg & df$Animal_ID==id,]
    prep_dat <- get_mol_knn(df_sub,knn_k=30)
    ncm_ct_sub <- prep_dat[[1]] # molecules by cell types
    ncm_ct_sub <- ncm_ct_sub[,all_ct]
    ncm_cell_sub <- prep_dat[[2]] # molecules by cells

    all_ncm_ct[[length(all_ncm_ct)+1]] <- ncm_ct_sub
    all_ncm_cell[[length(all_ncm_cell)+1]] <- ncm_cell_sub
  }
}
ncm_ct <- do.call(rbind,all_ncm_ct)

# combine ncm_cell matrices into one
all_row_names <- unique(unlist(lapply(all_ncm_cell, rownames)))
all_col_names <- unique(unlist(lapply(all_ncm_cell, colnames)))
row_map <- setNames(seq_along(all_row_names), all_row_names)
col_map <- setNames(seq_along(all_col_names), all_col_names)
triplet_list <- lapply(all_ncm_cell, function(m) {
  m_t <- as(m, "TsparseMatrix")  # Convert to triplet format (dgTMatrix)
  row_inds <- row_map[rownames(m)][m_t@i + 1]
  col_inds <- col_map[colnames(m)][m_t@j + 1]
  x_vals <- m_t@x  # Values remain the same
  data.frame(i = row_inds, j = col_inds, x = x_vals)
})
all_triplets <- do.call(rbind, triplet_list)
ncm_cell <- sparseMatrix(
  i = all_triplets$i,
  j = all_triplets$j,
  x = all_triplets$x,
  dims = c(length(all_row_names), length(all_col_names)),
  dimnames = list(all_row_names, all_col_names)
)

rm(knn_res)
rm(nmf_crf)
gc()

# result is saved in dir_save/bridge_nmftype_k.rds' for each method variant tested
run_bridge_multi(crf_all,df,cell_annot,ncm_ct,ncm_cell,dir_save=base_dir,
                 ncells_samp=500,n.cores=20)

### now run the marker enrichment tests
# result is saved in dir_save/enrich_results.rds'
enr_res <- run_enr_multi(all_nmf,markers,dir_save=base_dir,adj_pvals=TRUE,pthresh=.1,
                         n_perm=1000,n.cores=10)
