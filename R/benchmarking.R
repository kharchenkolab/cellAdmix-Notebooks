#' @importFrom sccore plapply
#' @import ggplot2
NULL

#' Computes admixture scores from
#'
#' @param score_params list of parameters including   signal_thres, exclude_genes,
#' p.c, min.expr.frac, exclude.cell.types, and max.prob.high
#' @param df Data Frame with transcript coordinates, gene, cell, and cell_type columns
#' @param so_rna Seurat object with scRNA-seq data for comparison
#' @param cell_type_adj_mat Cell-type adjecency matrix
#' @param get_sc_scores Logical whether to compute scores for the scRNA data as well
#'
#' @return a Data Frame of admixture scores for all cells
#' @export
get_dbl_scores <- function(score_params, df, so_rna, cell_type_adj_mat, get_sc_scores=FALSE) {
  so_spatial <- get_counts_meta_seurat(df)
  so_spatial <- Seurat::NormalizeData(so_spatial)

  common_genes <- intersect(rownames(so_rna), rownames(so_spatial))

  so_spatial %<>% subset(features=common_genes)
  so_rna %<>% subset(features=common_genes)

  signal_thres <- score_params[[1]]
  exclude_genes <- score_params[[2]]
  p.c <- score_params[[3]]
  min.expr.frac <- score_params[[4]]
  exclude.cell.types <- score_params[[5]]
  max.prob.high <- score_params[[6]]

  cont_info <- estimate_contamination_scores_seurat(
    so_rna, so_spatial, cell.type.adj.mat=cell_type_adj_mat,
    p.c=p.c, signal.thres=signal_thres, min.expr.frac=min.expr.frac, max.prob.high=max.prob.high,
    exclude.cell.types=exclude.cell.types, exclude.genes=exclude_genes
  )

  doublet_df <- cont_info$doublet.scores %>%
    {data.frame(scores=., cell_type=so_spatial$cell_type[names(.)], cell=names(.), dataset='spatial')}

  if (!get_sc_scores) {
    return(doublet_df)
  }

  cont_info_rna <- estimate_contamination_scores_seurat(
    so_rna, so_rna, cell.type.adj.mat=cell_type_adj_mat,
    p.c=0.25, signal.thres=signal_thres, min.expr.frac=0.05,
    exclude.genes=exclude_genes, exclude.cell.types=exclude.cell.types)

  doublet_df_rna <- cont_info_rna$doublet.scores %>%
    {data.frame(scores=., cell_type=so_rna$cell_type[names(.)], cell=names(.), dataset='scRNAseq')}

  return(list(doublet_df, doublet_df_rna))
}

#' @importFrom dplyr %>%
#' @import magrittr
#'
#' @export
get_scores_multi <- function(annot_res_all,orig_nms_all,score_params,df,so_rna,
                             crf_all,cell_adj_df,n.cores=1) {
  sc_ctypes <- unique(so_rna$cell_type)
  ndx_keep <- which(df$celltype %in% sc_ctypes)
  df <- df[ndx_keep,]

  cell_type_adj_mat <- estimate_cell_type_adjacency(cell_adj_df)

  scores_orig <- get_dbl_scores(score_params,df,so_rna,cell_type_adj_mat,get_sc_scores=TRUE)
  spatial_orig <- scores_orig[[1]]
  sc_orig <- scores_orig[[2]]
  spatial_orig$dataset <- 'spatial_orig'
  sc_orig$dataset <- 'scRNAseq'
  spatial_orig$nmf_type <- 'none'
  sc_orig$nmf_type <- 'none'

  scores_cln_all <- plapply(1:length(annot_res_all),function(res_ndx) {
    print(paste0('Running test ',res_ndx,' out of ',length(annot_res_all)))
    res <- annot_res_all[[res_ndx]]
    res_nm <- names(annot_res_all)[res_ndx]
    res_nm_orig <- orig_nms_all[res_ndx]
    if (!is.list(res)) {
      if (is.na(res[1])) {
        scores_cln <- spatial_orig
        scores_cln$dataset <- res_nm
        scores_cln$nmf_type <- res_nm_orig
        return(scores_cln)
      }
    }

    if (is.list(res)) {
      df_ct_cln <- list()
      for (ct in names(res)) {
        if (!(ct %in% df$celltype)) {
          next
        }
        df$factor <- crf_all[[res_nm_orig]][ndx_keep,ct]
        if (is.na(res[[ct]][1])) {
          df_ct <- df[df$celltype==ct,]
          df_ct$factor <- NULL
          df_ct_cln[[ct]] <- df_ct
        } else {
          f_rm <- sapply(res[[ct]],function(x){
            strsplit(x,split='_')[[1]][[1]]
          })
          f_rm <- as.numeric(f_rm)
          df_ct <- df[df$celltype==ct,]
          df_ct <- df_ct[!(df_ct$factor %in% f_rm), ]
          df_ct$factor <- NULL
          df_ct_cln[[ct]] <- df_ct
        }
      }
      df_cln <- do.call(rbind.data.frame,df_ct_cln)
    } else {
      df$factor <- crf_all[[res_nm_orig]][ndx_keep,1]
      df_cln <- df[!(paste(df$factor, df$celltype, sep = "_") %in% res), ]
    }
    scores_cln <- get_dbl_scores(score_params,df_cln,so_rna,cell_type_adj_mat,get_sc_scores=FALSE)
    scores_cln$dataset <- res_nm
    scores_cln$nmf_type <- res_nm_orig
    return(scores_cln)
  },mc.preschedule=FALSE,n.cores=n.cores,progress=TRUE)
  # })
  scores_final <- c(scores_cln_all,list(spatial_orig,sc_orig))
  scores_final_df <- do.call(rbind.data.frame,scores_final)
  return(scores_final_df)
}


#' @export
theme_legend <- function(position, background.alpha=0.5) {
  theme(
    legend.position=position, legend.justification=position,
    legend.background=element_rect(fill=alpha("white", background.alpha))
  )
}

#' @export
plot_ct_av_norm_scores <- function(scores_final_df,trim_level=NULL,min_mean_thresh=1e-3,dodge.width=NULL) {
  all_ds_nm <- unique(scores_final_df$dataset)
  all_ds_nm <- all_ds_nm[1:(length(all_ds_nm)-2)]

  pattern <- "([^_]+)(?:_([^0-9_]+))?_([0-9]+)_([^_]+)"
  df <- as.data.frame(
    strcapture(pattern, all_ds_nm,
               proto = list(mytype = character(), mygroup2 = character(), myk = integer(), mymethod = character()))
  )
  df$original <- all_ds_nm
  mytype_levels <- c("joint", "ct")
  mymethod_levels <- c("enr","memb","bridge")
  df$mytype <- factor(df$mytype, levels = mytype_levels)
  df$mymethod <- factor(df$mymethod, levels = mymethod_levels)
  df_sorted <- df[order(df$mytype, is.na(df$mygroup2), df$mygroup2, df$myk, df$mymethod), ]

  sorted_vec <- sapply(1:nrow(df_sorted),function(i){
    if (nchar(df_sorted[i,2])>0) {
      return(paste0(df_sorted[i,1],'_',df_sorted[i,2],'_',df_sorted[i,3]))
    } else {
      return(paste0(df_sorted[i,1],'_',df_sorted[i,3]))
    }
  })

  sorted_vec <- unique(c("none",sorted_vec))
  sorted_vec2 <- c("spatial_orig",df_sorted$original,"scRNAseq")
  scores_final_df$dataset <- factor(scores_final_df$dataset,levels=sorted_vec2)

  if (!is.null(trim_level)) {
    # now computing trimmed average scores per cell type
    doublet_df_agg <- scores_final_df %>%
      group_by(nmf_type, dataset, cell_type) %>%
      summarise(scores = mean(scores, trim=trim_level, na.rm = TRUE))
  } else {
    doublet_df_agg <- scores_final_df %>%
      group_by(nmf_type, dataset, cell_type) %>%
      summarise(scores = mean(scores, na.rm = TRUE))
  }

  norm_scores <- doublet_df_agg %>% filter(dataset=='spatial_orig') %$% setNames(scores, cell_type)
  doublet_df_agg <- doublet_df_agg %>%
    filter(norm_scores[cell_type] > min_mean_thresh) %>%
    mutate(scores = scores / norm_scores[cell_type])

  doublet_df_agg$nmf_type <- factor(doublet_df_agg$nmf_type, levels = sorted_vec)

  y_min <- 0
  y_max <- max(doublet_df_agg$scores)
  y_max <- y_max + (.1*y_max) # add some buffer

  p <- ggplot(doublet_df_agg) +
    geom_boxplot(aes(x=dataset, y=scores), outliers=FALSE) +
    ggbeeswarm::geom_beeswarm(aes(x=dataset, y=scores, color=cell_type),dodge.width = dodge.width) +
    facet_wrap(~ nmf_type, scales = "free_x", drop=TRUE, nrow=1) +
    coord_cartesian(ylim = c(y_min,y_max)) +
    theme_legend(position=c(1, 1)) +
    geom_hline(yintercept=1,color='red',linetype='dashed') +
    guides(color=guide_legend("Cell type",keywidth=0.1,
                              keyheight=0.15,
                              default.unit="inch")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(p)
}



localmoran <- function(x, W, conditional = TRUE, mlvar = TRUE) {
  n <- length(x)

  # Center the variable
  xx <- mean(x)
  z <- x - xx

  # Compute the spatial lag of z using the weight matrix
  lz <- as.vector(W %*% z)

  # Compute s2: variance measure (mlvar = TRUE uses division by n)
  if (mlvar) {
    s2 <- sum(z^2) / n
  } else {
    s2 <- sum(z^2) / (n - 1)
  }

  # Local Moran's I statistic for each observation
  Ii <- (z / s2) * lz

  # Compute Wi: sum of weights per observation; and Wi2: sum of squared weights per observation
  Wi <- rowSums(W)
  Wi2 <- rowSums(W^2)

  # Compute m2 (same as s2 with mlvar = TRUE)
  m2 <- sum(z^2) / n

  # Expected value of Ii under the null (conditional version)
  if (conditional) {
    EIi <- - (z^2 * Wi) / ((n - 1) * m2)
  } else {
    EIi <- - Wi / (n - 1)
  }

  # Compute b2, a moment term used in the variance calculation
  if (mlvar) {
    b2 <- (sum(z^4) / n) / (s2^2)
  } else {
    b2 <- (sum(z^4) / (n - 1)) / (s2^2)
  }

  # Compute the variance of local Moran's I (conditional version)
  if (conditional) {
    VarIi <- ((z / m2)^2) * (n / (n - 2)) *
      (Wi2 - (Wi^2 / (n - 1))) * (m2 - (z^2 / (n - 1)))
  } else {
    VarIi <- NA  # Not implemented for unconditional case in this replication
  }

  # Calculate the Z-score for each observation
  ZIi <- (Ii - EIi) / sqrt(VarIi)

  # Two-sided p-value based on the normal approximation
  p.value <- 2 * pnorm(abs(ZIi), lower.tail = FALSE)

  # Organize the results in a matrix similar to spdep's output
  res <- cbind(Ii, EIi, VarIi, ZIi, p.value)
  colnames(res) <- c("Ii", "E.Ii", "Var.Ii", "Z.Ii", "Pr(z != E(Ii))")
  return(res)
}

get_bg_cl_fracs <- function(score_params,df,so_rna,cell_adj_df,n_samp=5000,n.cores=5) {

  signal.thres <- score_params[[1]]
  p.c <- score_params[[3]]
  min.expr.frac <- score_params[[4]]
  exclude.cell.types <- score_params[[5]]
  max.prob.high <- score_params[[6]]

  # first need to compute contamination probabilities per gene for each ctype
  so_spatial <- get_counts_meta_seurat(df)
  so_spatial <- NormalizeData(so_spatial)

  cm.rna=so_rna[['RNA']]$counts
  cm.spatial=so_spatial[['RNA']]$counts

  common.genes <- intersect(rownames(cm.rna), rownames(cm.spatial))
  cm.rna <- cm.rna[common.genes,,drop=FALSE]
  cm.spatial <- cm.spatial[common.genes,,drop=FALSE]

  cell.type.adj.mat <- estimate_cell_type_adjacency(cell_adj_df)

  prob.gene.per.type <- estimate_gene_prob_per_type(cm.rna, cell.annot=as.factor(so_rna$cell_type),
                                                use.counts=FALSE)
  prior.cont.probs <- cell.type.adj.mat %>% {1 - diag(.) / colSums(.)}
  expr.frac.mask <- (matrixStats::colMaxs(prob.gene.per.type) > min.expr.frac)

  cont.probs.per.type <- colnames(cell.type.adj.mat) %>% setNames(., .) %>% setdiff(exclude.cell.types) %>%
    sapply(\(ct) {
      cont.probs <- prob.gene.per.type %>%
        cellAdmix:::estimate_gene_scores(cell.type.adj.mat, cell.type=ct) %>%
        cellAdmix:::estimate_gene_contamination_probabilities(p.c=p.c, prior.cont.prob=prior.cont.probs[ct])

      cont.probs[!((cont.probs > signal.thres) & expr.frac.mask[names(cont.probs)])] <- 0

      return(cont.probs)
    })

  all_ct <- unique(df$celltype)
  sig_frac_full <- list()
  for (ct in all_ct) {
    cont.probs <- cont.probs.per.type[,ct]
    cont.genes <- names(cont.probs)[cont.probs > signal.thres]

    cells_test <- unique(df$cell[df$celltype==ct])
    df_sub <- df[df$cell%in%cells_test,]

    cells_test = df_sub %>% dplyr::count(cell) %>% dplyr::filter(n > 50) %>% dplyr::pull(cell)

    if (length(cells_test)>n_samp) {
      cells_test <- sample(cells_test,n_samp)
    }
    df_sub <- df[df$cell%in%cells_test,]

    cell_sig <- plapply(cells_test,function(cell){
      print(cell)
      df_c <- df_sub[df_sub$cell==cell,]

      if (sum(df_c$gene %in% cont.genes)<5) {
        return(NA)
      }

      df_c$high_prob <- df_c$gene %in% cont.genes
      df_c$high_prob <- as.numeric(df_c$high_prob)

      # Create a matrix of coordinates
      coords <- cbind(df_c$x, df_c$y)

      # ##### using spdep package
      # # Define neighbors using, e.g., 4-nearest neighbors
      # knn <- spdep::knearneigh(coords, k = 20)
      # nb <- spdep::knn2nb(knn)
      #
      # # Convert neighbors into a spatial weights list (row-standardized weights)
      # lw <- spdep::nb2listw(nb, style = "W")
      #
      # # Compute Moran's I
      # moran_result <- spdep::moran.test(df_c$high_prob, lw)
      # m_est <- moran_result$estimate[1]
      # m_pv <- moran_result$p.value
      # #####

      ## without spdep package
      knn_result <- FNN::get.knn(coords, k = 20)

      # Build the weight matrix W (non-symmetric in this example)
      n <- nrow(coords)
      W <- matrix(0, n, n)
      for (i in 1:n) {
        neighbors <- knn_result$nn.index[i, ]
        # Optionally assign weights (here simply 1 for each neighbor)
        W[i, neighbors] <- 1
      }
      moran_result <- ape::Moran.I(df_c$high_prob, W)
      m_est <- moran_result$observed
      m_pv <- moran_result$p.value

      g_unq <- unique(df_c$gene[df_c$gene %in% cont.genes])
      if (m_est>0 & m_pv<.01) {
        # Compute local Moran's I
        # local_i <- spdep::localmoran(df_c$high_prob, lw)
        local_i <- localmoran(df_c$high_prob,W)
        local_i <- as.data.frame(local_i)
        local_i$gene <- df_c$gene
        # if there are multiple of same gene, just sample 1
        g_sig <- lapply(g_unq,function(gtest){
          ndx_samp <- which(local_i$gene==gtest)
          local_i_sub <- local_i[ndx_samp,c('Ii','Pr(z != E(Ii))'),drop=FALSE]
          if (any(local_i_sub[,'Ii'] > 0 & local_i_sub[,'Pr(z != E(Ii))']<.1)) { # use less stringent p here
            return(TRUE)
          } else {
            return(FALSE)
          }
        })
        names(g_sig) <- g_unq

      } else {
        g_sig <- rep(FALSE,length(g_unq))
        names(g_sig) <- g_unq
      }

      return(g_sig)
    },mc.preschedule=TRUE,n.cores=n.cores,progress=TRUE)
    # })
    names(cell_sig) <- cells_test

    cell_sig <- cell_sig[!is.na(cell_sig)]

    g_all <- unique(unlist(sapply(cell_sig,function(x){
      return(names(x))
    })))
    g_sig_frac_all <- c()
    for (g in g_all) {
      g_sig <- sapply(cell_sig,function(gvec) {
        if (g %in% names(gvec)) {
          return(gvec[g])
        } else {
          return(NA)
        }
      })
      g_sig <- unlist(g_sig[!is.na(g_sig)])
      g_sig_frac <- sum(g_sig)/length(g_sig)
      g_sig_frac_all <- c(g_sig_frac_all,g_sig_frac)
    }
    names(g_sig_frac_all) <- g_all

    sig_frac_full[[ct]] <- g_sig_frac_all
  }
  return(sig_frac_full)
}