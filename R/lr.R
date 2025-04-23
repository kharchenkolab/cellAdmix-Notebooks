#' @export
plot_lr_results <- function(spatial_all_scores,p_thresh,fc_thresh) {
  ##### making dotplot visualization for LR results
  ### first select ligands that have a significant result with at least one receptor
  selected_spat1 <- spatial_all_scores[spatial_all_scores$LR_cell_comb %in% c('macrophage--fibroblast'),]
  selected_spat1 <- selected_spat1[selected_spat1$lig_nr >= 2 & selected_spat1$rec_nr >= 2,]
  unq_lig <- unique(selected_spat1$ligand)
  lig_min_pvals <- sapply(unq_lig,function(lig){
    spat_sub <- selected_spat1[selected_spat1$ligand==lig,]
    spat_sub <- spat_sub[spat_sub$p.adj<p_thresh & spat_sub$log2fc > fc_thresh,]
    if (nrow(spat_sub)>0) {
      return('good')
    } else {
      return('bad')
    }
  })
  lig_sig_any <- names(lig_min_pvals)[lig_min_pvals=='good']
  selected_spat1 <- selected_spat1[selected_spat1$ligand %in% lig_sig_any,]

  spat_cast <- acast(selected_spat1, ligand ~ receptor, value.var = "pvalue", fill = NA)

  ## removing columns with less than N connections
  lr_counts <- table(selected_spat1[,c('ligand','receptor')])
  r_counts <- colSums(lr_counts)
  ndx_remove <- which(r_counts<3)
  r_rem <- names(ndx_remove)

  # dont remove if they're still significant
  selected_spat1_sub <- selected_spat1[selected_spat1$receptor%in%r_rem,]
  selected_spat1_sub <- selected_spat1_sub[selected_spat1_sub$p.adj>p_thresh,]
  selected_spat1_sub <- selected_spat1_sub[selected_spat1_sub$log2fc<fc_thresh,]
  r_rem <- unique(selected_spat1_sub$receptor)
  selected_spat1 <- selected_spat1[!(selected_spat1$receptor %in% r_rem),]

  # create an expression heatmap to show below the LR dotplot

  # first order genes by clustering them
  gene_order <- rownames(spat_cast)
  ct_sub <- ct_exp %>%
    filter(gene %in% gene_order)
  wide_df <- dcast(ct_sub, gene ~ celltype, value.var = "z")
  rownames(wide_df) <- wide_df$gene
  wide_df$gene <- NULL
  row_clusters <- hclust(as.dist(cor(t(wide_df))),method = 'complete')
  gene_order <- rownames(wide_df)[row_clusters[["order"]]]
  p_ct = ct_exp %>%
    filter(gene %in% gene_order) %>%
    mutate(gene = factor(gene, gene_order)) %>%
    ggplot(
      aes(y = celltype, x = gene, fill = z)
    ) +
    theme_classic() +
    geom_tile(show.legend = F) +
    scale_fill_gradient2(
      low = 'blue', mid = 'white',
      high = 'red',
      oob = scales::oob_squish,
      limits = c(-2, 2),
      na.value = 'whitesmoke'
    ) +
    ylab('Cell type') +
    xlab('scRNA expression') +
    scale_y_discrete(drop=FALSE, expand = expansion(add = 0), position = 'left') +
    scale_x_discrete(drop=FALSE, expand = expansion(add = 0), position = 'bottom') +
    p_theme +
    theme(
      axis.text.y = element_text(size = 5, angle = 0, vjust = 0.5, hjust = 1),
      panel.border = element_rect(fill = NA, color = 'gray'),
      axis.line = element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()
    )

  # reordering the lr matrix to match the clustered expression matrix
  spat_cast <- spat_cast[gene_order,]
  fib_marker_labs <- sapply(rownames(spat_cast), function(x){
    if (x %in% fibroblast.marker.genes) {
      return('fibroblast marker')
    } else {
      return('not fibroblast marker')
    }
  })

  # Create a named vector mapping categories to colors
  category_colors <- setNames(c('grey20','red4'), c('not fibroblast marker','fibroblast marker'))
  color_vec <- category_colors[fib_marker_labs]
  names(color_vec) <- rownames(spat_cast)

  selected_spat1$ligand <- factor(selected_spat1$ligand,levels=rownames(spat_cast))

  pl <- ggplot2::ggplot(data = selected_spat1) +
    ggplot2::geom_point(aes_string(x = 'ligand',
                                    y = 'receptor', size = 'p.adj', color = 'log2fc')) +
    xlab('ligand (macrophages)') +
    ylab('receptor (fibroblasts)') +
    scale_x_discrete(guide = guide_axis(angle = 45,position = 'top')) +
    scale_size_continuous(range = c(2, 0.5)) +
    scale_fill_gradient(low = 'lightgray', high = 'purple',aesthetics = 'color') +
    theme_light(base_line_size = gg_line_thickness) +
    p_theme +
    theme(axis.text.x = ggtext::element_markdown(color = color_vec),
          axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1,color = 'black'),
          plot.title = element_text(hjust = 0.5,vjust = -5,size = 8),
          axis.title=element_text(size=8),
          legend.text=element_text(size=6),
          legend.key.size = unit(6, 'pt'),
          legend.key.height = unit(6, 'pt'),
          legend.key.width = unit(6, 'pt'),
          legend.title=element_text(size=6))

  pl <- cowplot::plot_grid(pl,p_ct,ncol=1,align = 'v',rel_heights = c(.7,.3))

  ### hypergeometric test
  tester <- spatial_all_scores[spatial_all_scores$LR_cell_comb %in% c('macrophage--fibroblast'),]
  num_fib_tested <- length(unique(tester$ligand[tester$ligand %in% fibroblast.marker.genes]))
  bg_tested <- length(unique(tester$ligand))
  tester_sub <- tester[tester$p.adj<p_thresh,]
  tester_sub <- tester[abs(tester$log2fc)>fc_thresh,]
  num_fib_sig <- length(unique(tester_sub$ligand[tester_sub$ligand %in% fibroblast.marker.genes]))
  bg_sig <- length(unique(tester_sub$ligand))


  pval <- phyper(num_fib_sig-1, num_fib_tested, bg_tested-num_fib_tested, bg_sig, lower.tail=FALSE)

  return(list(pl,pval))
}