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
data <- prepare_OC_sc_spatial()
df <- data[[1]]
cell_annot <- data[[2]]
markers <- load_markers_OC()

base_dir <- CachePath('OC_scaled_dat6')

k_joint <- c(5,15,20,30)
k_ct <- c(5)


rm(data)
gc()

# load the precomputed NMF and CRF results (from NSCLC_nmf_crf_benchmark.R)
nmf_crf <- load_nmf_crf_multi(k_joint, k_ct, dir=base_dir)
all_nmf <- nmf_crf[[1]]
crf_all <- nmf_crf[[2]]

# ### run bridge tests first
# knn_res <- get_mol_knn_blocks(df,knn_k=10)

# saveRDS(knn_res,file=paste0(base_dir,'knn_res.rds'))
knn_res <- readRDS(file=paste0(base_dir,file='knn_res.rds'))

ncm_ct <- knn_res[[1]]
ncm_cell <- knn_res[[2]]

rm(knn_res)
gc()

run_bridge_multi(crf_all,df,cell_annot,ncm_ct,ncm_cell,dir_save=base_dir,
                 ncells_samp=500,n.cores=9)

### now running membrane test
# convert back to pixels from micrometers
pixel_size = DatasetPath('human_ovarian_cancer', 'experiment.xenium') %>%
  jsonlite::read_json() %$% pixel_size
# for the low res images
pixel_size <- pixel_size * 2
df_for_memb = df %>%
  mutate(
    x=as.integer(x / pixel_size),
    y=as.integer(y / pixel_size),
    z=as.integer(z / pixel_size)
  )
df_for_memb$z_index <- 0 # we don't have z-coordinate information for membrane stains
df_for_memb$fov <- 1 # stains aren't split up per fov in this dataset
df_for_memb <- df_for_memb[,c('x','y','z_index','z','cell','celltype','fov','mol_id')]

rm(df)
rm(knn_res)
gc()


# result is saved in dir_save/membrane_nmftype_k.rds' for each method variant tested
run_memb_multi(df_for_memb,crf_all,im_load_fn=load_images_oc,dir_save=base_dir,
               min_mols=5,max_cells_per_fov=200,
               knn_k=10,c_dist_expand=.05,ncm_ct=ncm_ct,ncm_cell=ncm_cell,n.cores=15)

### now run the marker enrichment tests
enr_res <- run_enr_multi(all_nmf,markers,dir_save=base_dir,adj_pvals=TRUE,pthresh=.1,n_perm=1000,
                         df=df,crf_h=40,n.cores=15)
