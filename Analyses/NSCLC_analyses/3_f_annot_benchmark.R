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
cell_annot <- prepare_nsclc_metadata(reps='one')
df <- prepare_nsclc_transcript_data(cell_annot,reps='one')
markers <- load_markers_nsclc()

base_dir <- CachePath('NSCLC_scaled_dat11')

# k_joint <- c(5,15,20,30)
# k_ct <- c(5)

k_joint <- c(5,15)
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
                 ncells_samp=500,n.cores=30)

### now running membrane test
df_for_memb <- df
# converting back to px from um units
df_for_memb$x <- (df_for_memb$x * 1000) / 180
df_for_memb$y <- (df_for_memb$y * 1000) / 180
df_for_memb$z <- (df_for_memb$z * 1000) / 800

df_for_memb <- df_for_memb[,c('x_local_px','y_local_px','z','cell','fov','celltype','gene')]
colnames(df_for_memb) <- c('x','y','z','cell','fov','celltype','gene')

# # remove z slices that we don't have images for
ndx_keep <- which(df_for_memb$z %in% c(0:8))
df_for_memb <- df_for_memb[ndx_keep,]
# also need to subset the CRF results to these same indices
crf_for_memb <- lapply(crf_all,function(x){
  return(x[ndx_keep,,drop=FALSE])
})

rm(crf_all)
rm(ncm_ct)
rm(ncm_cell)
gc()

# result is saved in dir_save/membrane_nmftype_k.rds' for each method variant tested
run_memb_multi(df_for_memb,crf_for_memb,im_load_fn=load_images_nsclc,dir_save=base_dir,
               min_mols=5,max_cells_per_fov=10,
               knn_k=10,c_dist_expand=.1,n.cores=30)


### now run the marker enrichment tests
# result is saved in dir_save/enrich_results.rds'
enr_res <- run_enr_multi(all_nmf,markers,dir_save=base_dir,adj_pvals=TRUE,pthresh=.1,n_perm=1000,n.cores=30)
