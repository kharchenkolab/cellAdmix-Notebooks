



# function to select markers from the full cell type DE results
#' @export
select_markers_nsclc <- function(marker_genes,ctype) {
  markers <- marker_genes[marker_genes$cluster==ctype,]
  markers <- markers[markers$p_val_adj<.01,]
  markers <- markers[markers$avg_log2FC>log2(1.5),]
  markers <- markers[markers$AUC>.6,]
  markers <- markers$Gene

  ct_all <- unique(marker_genes$cluster)
  num_expressing_all_ct <- data.frame(matrix(nrow=length(markers),ncol=length(ct_all)))
  rownames(num_expressing_all_ct) <- markers
  colnames(num_expressing_all_ct) <- ct_all
  for (ct in ct_all) {
    if (ctype=='malignant' & ct %in% c('Club cells','epithelial','Ciliated cells')) {
      next
    } else if (ctype=='macrophage' & ct %in% c('monocyte','DC','bMyeloid precursor-like','mast')) {
      next
    } else if (ctype=='fibroblast' & ct %in% c('Smooth muscle cells')) {
      next
    }
    if (ct!=ctype) {
      marker_genes_sub <- marker_genes[marker_genes$cluster==ct,]
      rownames(marker_genes_sub) <- marker_genes_sub$Gene
      marker_genes_sub <- marker_genes_sub[markers,]
      num_expressing_all_ct[,ct] <- marker_genes_sub$ExpressionFraction
    }
  }
  col_keep <- which(!is.na(num_expressing_all_ct[1,]))
  num_expressing_all_ct <- num_expressing_all_ct[,col_keep]
  marker_expression <- rowSums(num_expressing_all_ct>.2)
  cleaned.marker.genes <- names(marker_expression)[marker_expression==0]
  cleaned.marker.genes <- cleaned.marker.genes[cleaned.marker.genes %in% spatial_genes]
  return(cleaned.marker.genes)
}





select_markers_gut <- function(de_genes_full1,sc_obj) {
  # subset to significant results only
  de_genes_full1 <- de_genes_full1[de_genes_full1$padj<.05,]
  de_genes_full1 <- de_genes_full1[abs(de_genes_full1$Precision)>.3,]

  # order by AUC
  de_genes_full1 <- de_genes_full1[order(de_genes_full1$AUC,decreasing = TRUE),]

  ct_list <- lapply(unique(de_genes_full1$cell_type),function(ct){
    return(de_genes_full1[de_genes_full1$cell_type==ct,])
  })
  names(ct_list) <- unique(de_genes_full1$cell_type)

  # compute fraction of cells expressing each DE gene
  markers_plot <- unique(de_genes_full1$Gene)
  all_expr_counts <- list()
  for (ct in unique(Idents(sc_obj))) {
    cells_keep <- rownames(sc_obj@meta.data)[as.character(Idents(sc_obj))==ct]
    sc_obj1_sub <- subset(sc_obj,cells=cells_keep)
    expr <- sc_obj1_sub[['RNA']]$data[markers_plot,]
    cell_expr_counts <- Matrix::rowSums(expr>0) / ncol(expr)
    all_expr_counts[[ct]] <- cell_expr_counts
  }
  all_expr_counts <- do.call(cbind,all_expr_counts)
  for (ct in unique(de_genes_full1$cell_type)) {
    marker_df_sub <- de_genes_full1[de_genes_full1$cell_type==ct,]
    potential_markers <- marker_df_sub$Gene

    p <- DotPlot(object = sc_obj, features = potential_markers,
                 col.min = 0,col.max = 5)
    p$data$avg.exp.scaled <- NA
    p_dat_lst <- lapply(potential_markers,function(g) {
      p_dat_sub <- p$data[p$data$features.plot==g,]
      X <- p_dat_sub$avg.exp
      min_val <- min(X)
      max_val <- max(X)
      X_std = (X - min_val) / (max_val - min_val)
      p_dat_sub$avg.exp.scaled <- X_std
      return(p_dat_sub)
    })

    p_dat_lst2 <- do.call(rbind.data.frame,p_dat_lst)

    # now remove these if they are more highly expressed in other cell types
    g_rem_all <- c()
    for (mark in potential_markers) {
      de_sub_efrac <- all_expr_counts[mark,]
      if (ct=='CD4+ T cell') {
        de_sub_efrac <- de_sub_efrac[names(de_sub_efrac)!='CD8+ T cell']
      } else if (ct=='CD8+ T cell') {
        de_sub_efrac <- de_sub_efrac[names(de_sub_efrac)!='CD4+ T cell']
      } else if (ct=='B cell') {
        de_sub_efrac <- de_sub_efrac[names(de_sub_efrac)!='Plasma cell']
      } else if (ct=='Plasma cell') {
        de_sub_efrac <- de_sub_efrac[names(de_sub_efrac)!='B cell']
      }

      de_sub_efrac2_other <- de_sub_efrac[names(de_sub_efrac)!=ct]
      if (de_sub_efrac[ct]<max(de_sub_efrac2_other)) {
        g_rem_all <- c(g_rem_all,mark)
      }


      if (sum(de_sub_efrac2_other>.6)>1) {
        g_rem_all <- c(g_rem_all,mark)
      }

      scaled_vals <- p_dat_lst2[p_dat_lst2$features.plot==mark,]
      scaled_vals_nm <- scaled_vals[,'id']
      scaled_vals <- scaled_vals[,'avg.exp.scaled']
      names(scaled_vals) <- scaled_vals_nm
      scaled_vals_other <- scaled_vals[names(scaled_vals)!=ct]
      if (max(scaled_vals_other)>.8) {
        g_rem_all <- c(g_rem_all,mark)
      }

    }
    ct_dat <- ct_list[[ct]]
    ct_dat <- ct_dat[!(ct_dat$Gene %in% g_rem_all),]
    ct_list[[ct]] <- ct_dat
  }
  return(ct_list)
}


select_markers_brain <- function(de_genes_full1,sc_seurat_clean) {
  # subset to significant results only
  de_genes_full1 <- de_genes_full1[de_genes_full1$padj<.05,]
  de_genes_full1 <- de_genes_full1[abs(de_genes_full1$Precision)>.5,]

  # order by AUC
  de_genes_full1 <- de_genes_full1[order(de_genes_full1$AUC,decreasing = TRUE),]

  ct_list <- lapply(unique(de_genes_full1$cell_type),function(ct){
    return(de_genes_full1[de_genes_full1$cell_type==ct,])
  })
  names(ct_list) <- unique(de_genes_full1$cell_type)

  ### doing a version of the above but using expression not only if called a marker
  markers_plot <- unique(de_genes_full1$Gene)
  all_expr_counts <- list()
  for (ct in unique(Idents(sc_seurat_clean))) {
    cells_keep <- rownames(sc_seurat_clean@meta.data)[as.character(Idents(sc_seurat_clean))==ct]
    sc_seurat_clean_sub <- subset(sc_seurat_clean,cells=cells_keep)
    expr <- sc_seurat_clean_sub[['RNA']]$data[markers_plot,]
    cell_expr_counts <- Matrix::rowSums(expr>0) / ncol(expr)
    all_expr_counts[[ct]] <- cell_expr_counts
  }
  all_expr_counts <- do.call(cbind,all_expr_counts)
  for (ct in unique(de_genes_full1$cell_type)) {
    marker_df_sub <- de_genes_full1[de_genes_full1$cell_type==ct,]
    potential_markers <- marker_df_sub$Gene

    p <- DotPlot(object = sc_seurat_clean, features = potential_markers,
                 col.min = 0,col.max = 5)
    p$data$avg.exp.scaled <- NA
    p_dat_lst <- lapply(potential_markers,function(g) {
      p_dat_sub <- p$data[p$data$features.plot==g,]
      X <- p_dat_sub$avg.exp
      min_val <- min(X)
      max_val <- max(X)
      X_std = (X - min_val) / (max_val - min_val)
      p_dat_sub$avg.exp.scaled <- X_std
      return(p_dat_sub)
    })

    p_dat_lst2 <- do.call(rbind.data.frame,p_dat_lst)


    # now remove these if they are more highly expressed in other cell types
    g_rem_all <- c()
    for (mark in potential_markers) {
      de_sub_efrac <- all_expr_counts[mark,]
      de_sub_efrac2_other <- de_sub_efrac[names(de_sub_efrac)!=ct]
      if (de_sub_efrac[ct]<max(de_sub_efrac2_other)) {
        g_rem_all <- c(g_rem_all,mark)
      }

      scaled_vals <- p_dat_lst2[p_dat_lst2$features.plot==mark,]
      scaled_vals_nm <- scaled_vals[,'id']
      scaled_vals <- scaled_vals[,'avg.exp.scaled']
      names(scaled_vals) <- scaled_vals_nm
      scaled_vals_other <- scaled_vals[names(scaled_vals)!=ct]

    }
    ct_dat <- ct_list[[ct]]
    ct_dat <- ct_dat[!(ct_dat$Gene %in% g_rem_all),]
    ct_list[[ct]] <- ct_dat
  }
  return(ct_list)
}



