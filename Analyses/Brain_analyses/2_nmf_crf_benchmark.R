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

base_dir <- CachePath('brain_scaled_dat4/')

k_joint <- c(5,15,20,30)
k_ct <- c(5)


# run NMF for each k value
# only using a subset of the data here for simplicity
cells_keep <- rownames(cell_annot)[cell_annot$Animal_ID==1 & cell_annot$Bregma==0.21]
df2 <- df[df$cell %in% cells_keep,]
run_nmf_multi_k(df2,cell_annot,k_joint,k_ct,h=10,dir_save=base_dir,num_cells_samp=2000,
                nmol_dsamp_joint=10000,nmol_dsamp_ct=10000,n.cores=30)


# run CRF for each k value
run_crf_multi_k(df,cell_annot,k_joint,k_ct,h=10,dir_save=base_dir,same.label.ratio=5,
                normalize.by='gene',n.cores=15)

