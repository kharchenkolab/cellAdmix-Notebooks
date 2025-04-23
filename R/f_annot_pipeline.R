#' @import Matrix
NULL

run_bridge_multi <- function(crf_all,df,cell_annot,ncm_ct,ncm_cell,dir_save,
                             ncells_samp=200,knn_k=20,allow_target_dup=FALSE,n.cores=1) {
  for (res_nm in names(crf_all)) {
    # don't run if the result already exists
    if (file.exists(paste0(dir_save,'bridge_',res_nm,'.rds'))) {
      next
    }

    print(paste0('Running test ',which(names(crf_all)==res_nm),' out of ', length(crf_all)))
    crf_res <- crf_all[[res_nm]]
    bridge_raw <- run_bridge_test(
      df, crf_res, cell_annot, ncm_ct, ncm_cell,
      ncells.samp=ncells_samp, knn.k=knn_k, n.cores=n.cores
    )

    saveRDS(bridge_raw, paste0(dir_save,'bridge_',res_nm,'.rds'))
  }
}

extract_bridge_multi <- function(all_nmf,dir,all_ctypes,bridge_pthresh=.01,
                                 bridge_adj_pvals=TRUE) {

  all_bridge_res <- list()
  for (res_nm in names(all_nmf)) {
    res_nm_splt <- strsplit(res_nm,split='_')[[1]]
    k <- as.numeric(res_nm_splt[[length(res_nm_splt)]])
    nmf_type <- res_nm_splt[[1]]

    if (!file.exists(paste0(dir,'bridge_',res_nm,'.rds'))) {
      next
    }
    bridge_raw <- readRDS(file=paste0(dir,'bridge_',res_nm,'.rds'))

    # parsing bridge results
    bridge_res <- extract_bridge_res(
      bridge_raw, all_ctypes, n.factors=k, nmf.type=nmf_type,
      p.thresh=bridge_pthresh, adj.pvals=bridge_adj_pvals
    )
    all_bridge_res[[res_nm]] <- bridge_res
  }
  return(all_bridge_res)
}


run_memb_multi <- function(df_for_memb,crf_all,im_load_fn,dir_save,min_mols=5,max_cells_per_fov=10,
                           knn_k=10,use_block_knn=FALSE,all_fov=NULL,c_dist_expand=.025,
                           ncm_ct=NULL,ncm_cell=NULL,n.cores=1) {

  for (res_nm in names(crf_all)) {
    # don't run if the result already exists
    if (file.exists(paste0(dir_save,'membrane_',res_nm,'.rds'))) {
      next
    }

    print(paste0('Running test ',which(names(crf_all)==res_nm),' out of ', length(crf_all)))
    res_nm_splt <- strsplit(res_nm,split='_')[[1]]
    num_factors <- as.numeric(res_nm_splt[[length(res_nm_splt)]])
    crf_res <- crf_all[[res_nm]]

    memb_raw <- run_memb_test(
      df_for_memb, crf_res, im_load_fn, num_factors, min.mols=min_mols,
      max.cells.per.fov=max_cells_per_fov, knn.k=knn_k, all.fov=all_fov,
      use.block.knn=use_block_knn, c.dist.expand=c_dist_expand,
      ncm.ct=ncm_ct, ncm.cell=ncm_cell, n.cores=n.cores
    )

    saveRDS(memb_raw,file=paste0(dir_save,'membrane_',res_nm,'.rds'))
  }
}


extract_memb_multi <- function(all_nmf,all_ctypes,dir,memb_pthresh=.01,
                               memb_adj_pvals=TRUE) {

  all_memb_res <- list()
  for (res_nm in names(all_nmf)) {
    res_nm_splt <- strsplit(res_nm,split='_')[[1]]
    k <- as.numeric(res_nm_splt[[length(res_nm_splt)]])
    nmf_type <- res_nm_splt[[1]]

    if (!file.exists(paste0(dir,'membrane_',res_nm,'.rds'))) {
      next
    }
    memb_raw <- readRDS(file=paste0(dir,'membrane_',res_nm,'.rds'))

    # parsing bridge results
    memb_res <- extract_memb_res(
      memb_raw,all_ctypes,nmf_type,p.thresh=memb_pthresh, adj.pvals=memb_adj_pvals
    )

    all_memb_res[[res_nm]] <- memb_res
  }
  return(all_memb_res)
}



run_enr_multi <- function(all_nmf,markers,dir_save=NULL,adj_pvals=TRUE,pthresh=.1,n_perm=1000,
                        df=NULL,crf_h=NULL,n.cores=1) {
  enr_res_all <- list()
  for (res_nm in names(all_nmf)) {
    print(paste0('Running test ',which(names(all_nmf)==res_nm),' out of ', length(all_nmf)))
    res_nm_splt <- strsplit(res_nm,split='_')[[1]]
    k <- as.numeric(res_nm_splt[[length(res_nm_splt)]])
    nmf_type <- res_nm_splt[[1]]

    nmf_res <- all_nmf[[res_nm]]

    # run enrichment test
    if (nmf_type=='ct') {
      enr_res <- list()
      for (ct in names(nmf_res)) {
        nmf_res_ct <- nmf_res[[ct]]

        enr_res_ct <- get_enr(
          nmf_res_ct, markers, n.perm=n_perm, p.thresh=pthresh,
          adj.pvals=adj_pvals, crf.h=crf_h, df=df, n.cores=n.cores
        )

        enr_res_to <- enr_res_ct[[1]]
        enr_res_to_cln <- lapply(enr_res_to, \(x) {if (ct %in% x) {ct} else {NULL}})
        enr_res_ct[[1]] <- enr_res_to_cln
        enr_res[[ct]] <- enr_res_ct
      }
    } else {
      enr_res <- get_enr(
        nmf_res, markers, n.perm=n_perm, p.thresh=pthresh,
        adj.pvals=adj_pvals, crf.h=crf_h, df=df, n.cores=n.cores
      )
    }

    enr_res_all[[res_nm]] <- enr_res
  }

  if (!is.null(dir_save)) {
    saveRDS(enr_res_all,file=paste0(dir_save,'enrich_results.rds'))
  }

  return(enr_res_all)
}

check_fp_multi <- function(
    df, cell_annot, crf_all, enr_res_all=NULL, bridge_res_all=NULL,
    memb_res_all=NULL, do_clean=TRUE, knn_k=100, median_thresh=.1
  ) {
  ##### first get cell neighbor knn matrix
  meta_sub <- cell_annot[,c('x','y','z')]
  knn_mat <- knn_adjacency_matrix(meta_sub, k=knn_k)
  colnames(knn_mat) <- rownames(meta_sub)
  rownames(knn_mat) <- rownames(meta_sub)

  # now convert this matrix into a matrix of cell type that the molecule belongs to
  j = factor(cell_annot$celltype)
  ncm_ct_full = knn_mat %*% sparseMatrix(i=1:length(j), j=as.integer(j), x=rep(1, length(j)))
  colnames(ncm_ct_full) = levels(j)

  cell_types <- cell_annot %$% setNames(celltype, rownames(.))

  # now compute factor fracs for each ct
  ct_fracs_per_ct_all <- list()
  for (res_nm in names(crf_all)) {
    nm_splt <- strsplit(res_nm,split='_')[[1]]
    nmf_type <- nm_splt[[1]]
    crf_res <- crf_all[[res_nm]]

    if (nmf_type=='ct') {
      ct_fracs_per_ct <- list()
      for (ct in colnames(crf_res)) {
        df$factor <- crf_res[,ct]
        factor_counts <- table(df[,c('cell','factor')])
        cf_fracs <- factor_counts / Matrix::rowSums(factor_counts)
        ct_fracs_per_ct[[ct]] <- cf_fracs
      }
    } else {
      df$factor <- crf_res[,1]
      factor_counts <- table(df[,c('cell','factor')])
      cf_fracs <- factor_counts / Matrix::rowSums(factor_counts)
      ct_fracs_per_ct <- cf_fracs
    }
    ct_fracs_per_ct_all[[res_nm]] <- ct_fracs_per_ct
  }

  check_helper <- function(res_all,do_clean=TRUE,median_thresh=.1) {
    final_all <- lapply(names(res_all), \(res_nm) {
      nmf_type <- strsplit(res_nm,split='_')[[1]][[1]]
      res <- res_all[[res_nm]]
      ct_fracs_per_ct <- ct_fracs_per_ct_all[[res_nm]]
      if (nmf_type=='ct') {
        check_res <- check_f_rm_per_ct(
          cell_types, res, ct_fracs_per_ct, ncm_ct_full,
          do.clean=do_clean, median.thresh=median_thresh
        )
      } else {
        check_res <- check_f_rm(
          cell_types, res[[1]], res[[2]], ct_fracs_per_ct, ncm_ct_full,
          do.clean=do_clean, median.thresh=median_thresh
        )
      }
      return(check_res)
    })
    names(final_all) <- names(res_all)
    return(final_all)
  }

  orig_nms_all <- c()
  annot_res_all <- list()
  if (!is.null(enr_res_all)) {
    enr_final_all <- check_helper(enr_res_all,do_clean=TRUE,median_thresh=median_thresh) # always FP check enr tests
    orig_nms_all <- c(orig_nms_all,names(enr_final_all))
    names(enr_final_all) <- paste0(names(enr_final_all),'_enr')
    annot_res_all <- c(annot_res_all,enr_final_all)
  }

  if (!is.null(bridge_res_all)) {
    bridge_final_all <- check_helper(bridge_res_all,do_clean=do_clean,median_thresh=median_thresh)
    orig_nms_all <- c(orig_nms_all,names(bridge_final_all))
    names(bridge_final_all) <- paste0(names(bridge_final_all),'_bridge')
    annot_res_all <- c(annot_res_all,bridge_final_all)
  }

  if (!is.null(memb_res_all)) {
    memb_final_all <- check_helper(memb_res_all,do_clean=do_clean,median_thresh=median_thresh)
    orig_nms_all <- c(orig_nms_all,names(memb_final_all))
    names(memb_final_all) <- paste0(names(memb_final_all),'_memb')
    annot_res_all <- c(annot_res_all,memb_final_all)
  }

  return(list(annot_res_all,orig_nms_all))
}



