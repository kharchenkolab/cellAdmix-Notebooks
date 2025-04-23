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
data <- prepare_gut_tx_and_meta()
cell_annot <- data[[1]]
df <- data[[2]]
markers <- load_markers_gut()

base_dir <- CachePath('gut_scaled_dat6/')

k_joint <- c(5,15,20,30)
k_ct <- c(5)

gc()

# load the precomputed NMF and CRF results (from NSCLC_nmf_crf_benchmark.R)
nmf_crf <- load_nmf_crf_multi(k_joint, k_ct, dir=base_dir)
all_nmf <- nmf_crf[[1]]
crf_all <- nmf_crf[[2]]

### run bridge tests first
knn_res <- get_mol_knn(df,knn_k=20)
ncm_ct <- knn_res[[1]]
ncm_cell <- knn_res[[2]]

rm(knn_res)
rm(nmf_crf)
gc()

# result is saved in dir_save/bridge_nmftype_k.rds' for each method variant tested
run_bridge_multi(crf_all,df,cell_annot,ncm_ct,ncm_cell,dir_save=base_dir,
                 ncells_samp=200,knn_k=20,n.cores=40)

## now running membrane test
df$fov <- 1
df$z_index <- dense_rank(df$z)-1 # need z_index for indexing the stain image. Code also assumes it's 0 indexed.
df_for_memb <- df[,c('x','y','z_index','z','cell','celltype','fov','mol_id','gene')]

rm(ncm_ct)
rm(ncm_cell)
gc()

# result is saved in dir_save/membrane_nmftype_k.rds' for each method variant tested
run_memb_multi(df_for_memb,crf_all,im_load_fn=load_images_gut,dir_save=base_dir,
               min_mols=5,max_cells_per_fov=200,
               knn_k=20,n.cores=40)


### now run the marker enrichment tests
# result is saved in dir_save/enrich_results.rds'
enr_res <- run_enr_multi(all_nmf,markers,dir_save=base_dir,adj_pvals=TRUE,pthresh=.1,n_perm=1000,
                         df=df,crf_h=20,n.cores=10)
