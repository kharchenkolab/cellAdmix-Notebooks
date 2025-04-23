suppressPackageStartupMessages({
  library(cellAdmix)
  library(cowplot)
  library(dataorganizer)
  library(ggplot2)
  library(magrittr)
  library(Matrix)
  library(readr)
  library(Seurat)
  library(gtools)
  library(dplyr)

  devtools::load_all()
})

theme_set(theme_bw())


cell_types <- CachePath('gut_cell_types_transferred.csv') %>% read_csv()
cell_types <- as.data.frame(cell_types)

load_df_seg <- function(method,min_count=10) {
  if (method=='Baysor') {
    data <- prepare_gut_tx_and_meta()
    cell_annot <- data[[1]]
    df <- data[[2]]
    return(list(df,cell_annot))
  }

  if (method=='ComSeg') {
    df <- DatasetPath('mouse_gut', 'comseg_format', 'results', '11_d6_h5_min20_s7_r1413', 'segmentation.csv') %>%
      read_csv(show_col_types=FALSE) %>% rename(cell=cell_index_pred)
    df <- as.data.frame(df)
    df$cell <- paste0('ComSeg_',df$cell)
    ct_annot <- cell_types[cell_types$protocol==method,]
    ct_annot_bcodes <- ct_annot$cell
    ndx_match <- match(df$cell,ct_annot_bcodes)
    df <- df[!is.na(ndx_match),]
    ndx_match <- match(df$cell,ct_annot_bcodes)
    df$celltype <- ct_annot[ndx_match,'cell_type']
    df$cell <- as.character(df$cell)
  } else if (method=='ProSeg') {
    df <- DatasetPath('mouse_gut', 'seg_method_results', 'proseg', 'transcript-metadata.csv') %>%
      read_csv(show_col_types=FALSE) %>%
      mutate(cell=ifelse(assignment == 4294967295, 0, assignment + 1))
    df <- as.data.frame(df)
    df$cell <- paste0('ProSeg_',df$cell)
    ct_annot <- cell_types[cell_types$protocol==method,]
    ct_annot_bcodes <- ct_annot$cell
    ndx_match <- match(df$cell,ct_annot_bcodes)
    df <- df[!is.na(ndx_match),]
    ndx_match <- match(df$cell,ct_annot_bcodes)
    df$celltype <- ct_annot[ndx_match,'cell_type']
    df$cell <- as.character(df$cell)
  }

  cell_annot <- df %>%
    group_by(cell) %>%
    summarise(x = mean(x, na.rm = TRUE),
              y = mean(y, na.rm = TRUE),
              celltype=unique(celltype)
    )
  cell_annot <- as.data.frame(cell_annot)
  rownames(cell_annot) <- cell_annot$cell

  # removing cells with very few molecules
  cell_counts <- table(df$cell)
  cells_keep <- names(cell_counts)[cell_counts>=min_count]
  df <- df[df$cell %in% cells_keep,]
  cell_annot <- cell_annot[rownames(cell_annot) %in% cells_keep,]

  df$mol_id <- as.character(1:nrow(df))
  rownames(df) <- df$mol_id

  return(list(df,cell_annot))
}


base_dir <- CachePath('gut_scaled_dat6/')

methods_test <- c('ComSeg', 'ProSeg') # baysor is done already as it's the default

run_nmf_crf_all_seg(methods_test,k=20,h_nmf=10,h_crf=10,base_dir,num_cells_samp=2000,
                    nmol_dsamp=20000,top_g_thresh=NULL,n_hvgs=NULL,mols_test=NULL,proj_h=NULL,bridge_ncells_samp=200,bridge_knn_k=20,
                    n.cores=30)

## now cleaning the data and saving these versions
methods_test <- c('ComSeg', 'ProSeg', 'Baysor')
k <- 20
for (method in methods_test) {
  print(method)
  df_dat <- load_df_seg(method)
  df <- df_dat[[1]]
  cell_annot <- df_dat[[2]]

  all_ctypes <- unique(cell_annot$celltype)

  if (method=='Baysor') {
    bridge_raw <- readRDS(file=paste0(base_dir,'bridge_joint_20.rds'))
  } else {
    bridge_raw <- readRDS(file=paste0(base_dir,'bridge_',method,'.rds'))
  }

  # parsing bridge results
  bridge_res <- extract_bridge_res(bridge_raw,all_ctypes,k,nmf.type='joint',p.thresh=.3,
                                   adj.pvals=FALSE)

  if (method=='Baysor') {
    crf_res <- readRDS(file=paste0(base_dir,'crf_joint_k',k,'.rds'))
  } else {
    crf_res <- readRDS(file=paste0(base_dir,'crf_joint_k',k,'_',method,'.rds'))
  }

  cell_annot$z <- 1
  crf_all <- list(crf_res)
  names(crf_all) <- method
  bridge_res_all <- list(bridge_res)
  names(bridge_res_all) <- method
  fp_checks <- check_fp_multi(df,cell_annot,crf_all,enr_res_all=NULL,bridge_res_all=bridge_res_all,
                              memb_res_all=NULL,do_clean=TRUE,knn_k=100,median_thresh=.2)

  annot_res_all <- fp_checks[[1]]
  orig_nms_all <- fp_checks[[2]]

  df$factor <- crf_res[,1]
  df_cln <- df[!(paste(df$factor, df$celltype, sep = "_") %in% annot_res_all[[1]]), ]
  saveRDS(df_cln,file=paste0(base_dir,'df_k',k,'_cln_',method,'.rds'))
}
