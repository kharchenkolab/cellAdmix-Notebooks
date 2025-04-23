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


base_dir <- CachePath('OC_scaled_dat6')

k_joint <- c(5,15,20,30)
k_ct <- c(5)


# run NMF for each k value
run_nmf_multi_k(df,cell_annot,k_joint,k_ct,h=5,dir_save=base_dir,num_cells_samp=2000,
                nmol_dsamp_joint=10000,nmol_dsamp_ct=10000,top_g_thresh=.25,
                n.cores=30)

rm(data)
gc()

# run CRF for each k value
run_crf_multi_k(df,cell_annot,k_joint,k_ct,h=10,dir_save=base_dir,same.label.ratio=5,
                normalize.by='gene',proj_h=50,n.cores=30)
