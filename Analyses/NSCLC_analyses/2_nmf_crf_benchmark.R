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
so_rna <- prepare_nsclc_scrna()
markers <- load_markers_nsclc()

base_dir <- CachePath('NSCLC_scaled_dat11')

k_joint <- c(5,15,20,30)
k_ct <- c(5)

# run NMF for each k value
run_nmf_multi_k(df,cell_annot,k_joint,k_ct,h=20,dir_save=base_dir,num_cells_samp=2000,
                nmol_dsamp_joint=10000,nmol_dsamp_ct=4000,n.cores=20)


# run CRF for each k value
run_crf_multi_k(df,cell_annot,k_joint,k_ct,h=10,dir_save=base_dir,same.label.ratio=5,
                normalize.by='gene',n.cores=20)
