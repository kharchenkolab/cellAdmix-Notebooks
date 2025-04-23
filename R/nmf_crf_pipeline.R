#' @importFrom dataorganizer ReadOrCreate
NULL

### NMF pipeline first

run_nmf_multi_k_joint <- function(df,cell_annot,k_joint,h,dir_save,num_cells_samp=NULL,
                            nmol_dsamp=10000,top_g_thresh=NULL,n_hvgs=NULL,mols_test=NULL,
                            n.cores=1) {
  if (!is.null(num_cells_samp)) {
    df <- samp_ct_equal(df,cell_annot,num_cells_samp)
  }

  if (!is.null(top_g_thresh) | !is.null(n_hvgs)) {
    df <- subset_genes(df,cell_annot,top_g_thresh=top_g_thresh,n_hvgs=n_hvgs)
  }

  ## run NMF for all cell types together aka joint
  for (k in k_joint) {
    if (file.exists(paste0(dir_save,'nmf_joint_k',k,'.rds'))) {
      next
    }

    message(paste0('starting NMF for joint k=',k))
    res <- run_knn_nmf(df,k=k,h=h,nmol_dsamp=nmol_dsamp,mols_test=mols_test,n.cores=n.cores)
    saveRDS(res,file=paste0(dir_save,'nmf_joint_k',k,'.rds'))
  }

  return(NULL)
}


run_nmf_multi_k_ct <- function(df,cell_annot,k_ct,h,dir_save,num_cells_samp=NULL,
                            nmol_dsamp=10000,top_g_thresh=NULL,n_hvgs=NULL,
                            n.cores=1) {

  if (!is.null(top_g_thresh) | !is.null(n_hvgs)) {
    df <- subset_genes(df,cell_annot,top_g_thresh=top_g_thresh,n_hvgs=n_hvgs)
  }

  ## running NMF separately per cell type
  # create the per cell type subdirectory if it doesn't already exist
  if (!dir.exists(paste0(dir_save,'per_ct_nmfs/'))) {
    dir.create(paste0(dir_save,'per_ct_nmfs/'))
  }

  all_ct <- unique(cell_annot$celltype)

  for (k in k_ct) {
    for (ct in all_ct) {
      message(paste0('starting NMF for ct ',which(all_ct==ct),' out of ',length(all_ct)))

      # subset df to just have the cell type
      df_ct <- df[df$celltype==ct,]

      n_cells <- length(unique(df_ct$cell))
      if (n_cells>num_cells_samp) {
        cells_samp <- sample(unique(df_ct$cell),num_cells_samp)
        df_ct <- df_ct[df_ct$cell %in% cells_samp,]
      }

      if (file.exists(paste0(dir_save,'per_ct_nmfs/k',k,'_',ct,'.rds'))) {
        next
      }

      res_ct <- run_knn_nmf(df_ct,k=k,h=h,nmol_dsamp=nmol_dsamp,n.cores=n.cores)
      saveRDS(res_ct,file=paste0(dir_save,'per_ct_nmfs/k',k,'_',ct,'.rds'))

    }
  }
  return(NULL)
}


run_nmf_multi_k <- function(df,cell_annot,k_joint,k_ct,h,dir_save,num_cells_samp=NULL,
                            nmol_dsamp_joint=10000,nmol_dsamp_ct=4000,top_g_thresh=NULL,
                            n_hvgs=NULL,mols_test=NULL,n.cores=1) {

  run_nmf_multi_k_joint(df,cell_annot,k_joint,h=h,dir_save=dir_save,num_cells_samp=num_cells_samp,
                        nmol_dsamp=nmol_dsamp_joint,top_g_thresh=top_g_thresh,n_hvgs=n_hvgs,
                        mols_test=mols_test,n.cores=n.cores)

  run_nmf_multi_k_ct(df,cell_annot,k_ct,h=h,dir_save=dir_save,num_cells_samp=num_cells_samp,
                     nmol_dsamp=nmol_dsamp_ct,top_g_thresh=top_g_thresh,n_hvgs=n_hvgs,n.cores=n.cores)

  return(NULL)
}





## CRF pipeline

run_crf_multi_k_joint <- function(df,cell_annot,k_joint,h,dir_save,same.label.ratio=5,
                                  normalize.by='gene',proj_h=NULL,n.cores=1) {
  ## run NMF for all cell types together aka joint
  for (k in k_joint) {
    # check if the result already exists
    if (file.exists(paste0(dir_save,'crf_joint_k',k,'.rds'))) {
      next
    }

    # require that a corresponding nmf result exists
    if (!file.exists(paste0(dir_save,'nmf_joint_k',k,'.rds'))) {
      stop(paste0('Missing NMF result for k=',k))
    }

    res <- readRDS(file=paste0(dir_save,'nmf_joint_k',k,'.rds'))

    message(paste0('starting CRF for joint k=',k))
    crf_res <- run_crf_all(df,res,num_nn=h,same.label.ratio=same.label.ratio,
                           normalize.by=normalize.by,proj_h=proj_h,n.cores=n.cores)
    saveRDS(crf_res,file=paste0(dir_save,'crf_joint_k',k,'.rds'))
  }

  return(NULL)
}


run_crf_multi_k_ct <- function(df,cell_annot,k_ct,h,dir_save,same.label.ratio=5,
                               normalize.by='gene',proj_h=NULL,n.cores=1) {

  all_ct <- unique(cell_annot$celltype)

  for (k in k_ct) {
    if (file.exists(paste0(dir_save,'crf_per_ct_k',k,'.rds'))) {
      next
    }

    for (ct in all_ct) {
      message(paste0('starting CRF for ct ',which(all_ct==ct),' out of ',length(all_ct)))

      if (file.exists(paste0(dir_save,'per_ct_nmfs/crf_k',k,'_',ct,'.rds'))) {
        next
      }

      # require that a corresponding nmf result exists
      if (!file.exists(paste0(dir_save,'per_ct_nmfs/k',k,'_',ct,'.rds'))) {
        stop(paste0('Missing NMF result for k=',k,' ',ct))
      }

      res <- readRDS(file=paste0(dir_save,'per_ct_nmfs/k',k,'_',ct,'.rds'))

      crf_res <- run_crf_all(df,res,num_nn=h,same.label.ratio=same.label.ratio,
                             normalize.by=normalize.by,proj_h=proj_h,n.cores=n.cores)
      saveRDS(crf_res,file=paste0(dir_save,'per_ct_nmfs/crf_k',k,'_',ct,'.rds'))

      rm(crf_res)
      gc()
    }

    # combine the crfs for the current k value into one df and save it
    files <- list.files(path = paste0(dir_save,'per_ct_nmfs'), pattern = paste0("^crf_k",k), full.names = TRUE)
    fnms <- sapply(files,function(x){
      fsplt <- strsplit(basename(x),split='_')[[1]]
      ct <- paste0(fsplt[3:length(fsplt)],collapse = '_')
      ct <- substr(ct,1,nchar(ct)-4)
      return(ct)
    })
    names(files) <- fnms

    crf_all_ct <- list()
    for (ct in names(files)) {
      f <- files[ct]
      crf_res_ct <- readRDS(f)
      crf_all_ct[[ct]] <- crf_res_ct
    }
    crf_res <- do.call(cbind.data.frame,crf_all_ct)
    colnames(crf_res) <- names(files)
    saveRDS(crf_res,file=paste0(dir_save,'crf_per_ct_k',k,'.rds'))
  }

  return(NULL)
}



run_crf_multi_k <- function(df,cell_annot,k_joint,k_ct,h,dir_save,same.label.ratio=5,
                            normalize.by='gene',proj_h=NULL,n.cores=1) {

  run_crf_multi_k_joint(df,cell_annot,k_joint,h,dir_save,same.label.ratio=same.label.ratio,
                        normalize.by=normalize.by,proj_h=proj_h,n.cores=n.cores)

  run_crf_multi_k_ct(df,cell_annot,k_ct,h,dir_save,same.label.ratio=same.label.ratio,
                     normalize.by=normalize.by,proj_h=proj_h,n.cores=n.cores)

  return(NULL)
}




load_nmf_crf_multi <- function(k_joint, k_ct, dir) {
  all_nmf <- list()
  crf_all <- list()
  for (k in k_joint) {
    res_nm <- paste0('joint_',k)
    res <- readRDS(paste0(dir,'nmf_joint_k',k,'.rds'))
    all_nmf[[res_nm]] <- res

    crf_res <- readRDS(paste0(dir,'crf_joint_k',k,'.rds'))
    crf_all[[res_nm]] <- crf_res
  }
  for (k in k_ct) {
    res_nm <- paste0('ct_',k)
    pattern <- paste0("^k", k, "_.*\\.rds$")
    rds_files <- list.files(paste0(dir,'per_ct_nmfs'), pattern = pattern, full.names = TRUE)
    rds_list <- lapply(rds_files, readRDS)
    names(rds_list) <- sub("^k[0-9]+_(.+)\\.rds$", "\\1", basename(rds_files))
    all_nmf[[res_nm]] <- rds_list

    crf_res <- readRDS(paste0(dir,'crf_per_ct_k',k,'.rds'))
    crf_all[[res_nm]] <- crf_res
  }

  return(list(all_nmf,crf_all))
}


#' Run NMF-CRF analysis for all segmentation methods
#'
#' @param methods_test Vector of method names to analyze
#' @param k Number of factors for NMF
#' @param h_nmf Parameter for NMF
#' @param h_crf Parameter for CRF
#' @param dir_save Directory to save results
#' @param num_cells_samp Number of cells to sample per cell type
#' @param nmol_dsamp Number of molecules to downsample
#' @param top_g_thresh Threshold for top genes
#' @param n_hvgs Number of highly variable genes
#' @param mols_test Molecules to test
#' @param proj_h Projection parameter
#' @param bridge_ncells_samp Number of cells to sample for bridge test
#' @param bridge_knn_k Number of nearest neighbors for bridge test
#' @param n.cores Number of cores to use
#' @return NULL (saves results to files)
run_nmf_crf_all_seg <- function(
    methods_test, k, h_nmf, h_crf, dir_save,
    num_cells_samp=NULL, nmol_dsamp=10000, top_g_thresh=NULL,
    n_hvgs=NULL, mols_test=NULL, proj_h=NULL,
    bridge_ncells_samp=500, bridge_knn_k=20, n.cores=1
  ) {

  ## Run NMF for all cell types together (joint analysis)
  for (method in methods_test) {
    # Skip the method if already run
    if (file.exists(paste0(dir_save, 'bridge_', method, '.rds'))) {
      next
    }

    df_dat <- load_df_seg(method)
    df <- df_dat[[1]]
    cell_annot <- df_dat[[2]]

    if (!is.null(num_cells_samp)) {
      df_sub <- samp_ct_equal(df, cell_annot, num_cells_samp)
    } else {
      df_sub <- df
    }

    if (!is.null(top_g_thresh) | !is.null(n_hvgs)) {
      df_sub <- subset_genes(df_sub, cell_annot, top.g.thresh=top_g_thresh, n.hvgs=n_hvgs)
    }

    message(paste0('Starting NMF for method ', method))
    nmf_file <- paste0(dir_save, 'nmf_joint_k', k, '_', method, '.rds')

    ReadOrCreate(nmf_file, {
      run_knn_nmf(
        df_sub, k=k, h=h_nmf, nmol.dsamp=nmol_dsamp,
        mols.test=mols_test, n.cores=n.cores
      )
    })


    rm(df_sub)
    gc()

    # Run CRF
    crf_file <- paste0(dir_save, 'crf_joint_k', k, '_', method, '.rds')

    ReadOrCreate(crf_file, {
      run_crf_all(
        df, res, num.nn=h_crf, same.label.ratio=5,
        normalize.by='gene', proj.h=proj_h, n.cores=n.cores
      )
    })

    # Run bridge test
    bridge_file <- paste0(dir_save, 'bridge_', method, '.rds')

    ReadOrCreate(bridge_file, {
      knn_res <- get_mol_knn_blocks(df, knn_k=10)
      ncm_ct <- knn_res[[1]]
      ncm_cell <- knn_res[[2]]

      run_bridge_test(
        df, crf_res, cell_annot, ncm_ct, ncm_cell,
        ncells.samp=bridge_ncells_samp,
        knn.k=bridge_knn_k, n.cores=n.cores
      )
    })
  }

  return(NULL)
}