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
data <- prepare_BC_sc_spatial()
df <- data$df_spatial
cell_annot <- data$cell_annot
markers <- load_markers_BC()

base_dir <- CachePath('BC_scaled_dat4/')
dir.create(base_dir, showWarnings=FALSE, recursive=TRUE)

k_joint <- c(5,15,20,30)
k_ct <- c(5)

gc()

# load the precomputed NMF and CRF results (from NSCLC_nmf_crf_benchmark.R)
nmf_crf <- load_nmf_crf_multi(k_joint, k_ct, dir=base_dir)
all_nmf <- nmf_crf[[1]]
crf_all <- nmf_crf[[2]]

### run bridge tests first
knn_res <- get_mol_knn(df,knn_k=10)
ncm_ct <- knn_res[[1]]
ncm_cell <- knn_res[[2]]

rm(knn_res)
rm(nmf_crf)
gc()

# result is saved in dir_save/bridge_nmftype_k.rds' for each method variant tested
run_bridge_multi(crf_all,df,cell_annot,ncm_ct,ncm_cell,dir_save=base_dir,
                 ncells_samp=500,n.cores=15)

### now run the marker enrichment tests
# result is saved in dir_save/enrich_results.rds'
enr_res <- run_enr_multi(all_nmf,markers,dir_save=base_dir,adj_pvals=TRUE,
                         pthresh=.1,n_perm=1000,n.cores=10)
