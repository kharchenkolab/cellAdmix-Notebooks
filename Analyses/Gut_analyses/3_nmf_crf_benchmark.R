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

base_dir <- CachePath('gut_scaled_dat6/')


k_joint <- c(5,15,20,30)
k_ct <- c(5)

# run NMF for each k value
run_nmf_multi_k(df,cell_annot,k_joint,k_ct,h=10,dir_save=base_dir,num_cells_samp=2000,
                nmol_dsamp_joint=20000,nmol_dsamp_ct=4000,n.cores=20)

# run CRF for each k value
run_crf_multi_k(df,cell_annot,k_joint,k_ct,h=10,dir_save=base_dir,same.label.ratio=5,
                normalize.by='gene',n.cores=30)
