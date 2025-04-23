#' @import dplyr
#' @importFrom magrittr %>% %<>% %$% set_colnames
#' @import Seurat
NULL

#' @export
downsample_spat <- function(so_spat,sc_obj) {
  dat_orig <- so_spat[['RNA']]$counts
  tb <- as_tibble(dat_orig)
  dat_orig_tmp <- lapply(tb, function(x) setNames(x, rownames(dat_orig)))
  dat_orig2 <- dat_orig_tmp
  meta_orig <- so_spat@meta.data
  all_ct <- unique(meta_orig$cell_type)
  num_dsampled_all <- c()
  for (ct in all_ct) {
    num_dsampled <- 0
    cells_keep <- rownames(meta_orig)[meta_orig$cell_type==ct]
    dat_sub <- dat_orig[,cells_keep]
    cells_keep <- rownames(sc_obj@meta.data)[sc_obj@meta.data$cell_type==ct]
    sc_sub <- subset(sc_obj,cells=cells_keep)
    lib_sub_spat <- sparseMatrixStats::colSums2(dat_sub)
    lib_sub_sc <- sparseMatrixStats::colSums2(sc_sub[['RNA']]$counts)
    names(lib_sub_spat) <- colnames(dat_sub)
    names(lib_sub_sc) <- colnames(sc_sub[['RNA']]$counts)
    lib_sub_spat <- lib_sub_spat[order(lib_sub_spat,decreasing = TRUE)]
    lib_sub_sc <- lib_sub_sc[order(lib_sub_sc,decreasing = TRUE)]
    sc_quantiles <- sapply(1:length(lib_sub_sc),function(i) {
      return((length(lib_sub_sc)-i+1)/length(lib_sub_sc))
    })
    for (i in 1:length(lib_sub_spat)) {
      spat_cell_nm <- names(lib_sub_spat)[i]
      percentile_cell <- (length(lib_sub_spat)-i+1)/length(lib_sub_spat)
      p_diff <- abs(sc_quantiles-percentile_cell)
      sc_select <- which(p_diff==min(p_diff))[1]
      if (lib_sub_spat[i]<lib_sub_sc[sc_select]) {
        next
      } else {
        num_dsampled <- num_dsampled + 1
        scaling_factor <- lib_sub_sc[sc_select] / lib_sub_spat[i]
        dat_orig2[[spat_cell_nm]] <- rbinom(nrow(dat_orig), dat_orig_tmp[[spat_cell_nm]], scaling_factor)
      }
    }
    num_dsampled_all <- c(num_dsampled_all,num_dsampled)
  }
  names(num_dsampled_all) <- all_ct
  print('Number of cells that had lib sizes downsampled:')
  print(num_dsampled_all)
  dat_orig2 <- do.call(cbind,dat_orig2)
  seurat_spatial <- CreateSeuratObject(dat_orig2,meta.data = meta_orig)
  seurat_spatial <- NormalizeData(seurat_spatial)
  Idents(seurat_spatial) <- as.factor(seurat_spatial@meta.data$cell_type)
  return(seurat_spatial)
}

#' @export
estimateAdmixtureFractions <- function(cell.assignment, is.admixture, breaks=c(0.02, 0.05, 0.1, 1)) {
  mask = (!is.na(cell.assignment) & !is.na(is.admixture))
  cell.assignment = cell.assignment[mask]
  is.admixture = is.admixture[mask]

  table(cell.assignment, is.admixture) %>%
    {.[,2] / rowSums(.)} %>%
    {data.frame(cell=names(.), admix_frac=.)} %>%
    mutate(admix_cat = cut(admix_frac, breaks=breaks, include.lowest=TRUE)) %>%
    tibble::as_tibble()
}

#' @export
estimateAdmixtureComparisonDf <- function(assay.dfs, cell.type.include=NULL) {
  comp.df <- assay.dfs %>%
    bind_rows(.id="assay") %>%
    mutate(assay=factor(assay, names(assay.dfs)))

  if (!is.null(cell.type.include)) {
    comp.df %<>% filter(cell_type %in% cell.type.include)
  }

  comp.df %<>%
    dplyr::count(cell_type, admix_cat, assay) %>%
    group_by(cell_type, assay) %>%
    mutate(
      frac = n/sum(n),
      frac_total = sum(frac[!is.na(admix_cat)])
    ) %>%
    ungroup() %>%
    arrange(-frac_total) %>%
    mutate(cell_type = factor(cell_type, unique(cell_type))) %>%
    na.omit()

  return(comp.df)
}

#' @export
transform_sc <- function(sc_obj) {
  sc_obj@meta.data$cell <- rownames(sc_obj@meta.data)
  df_sc <- as.data.frame(as.table(as.matrix(sc_obj[['RNA']]$counts))) %>%
    filter(Freq > 0) %>%          # Optionally remove rows with zero count
    tidyr::uncount(weights = Freq) %>%   # Expand rows by the count in Freq
    rename(gene = Var1, cell = Var2)
  return(df_sc)
}

#' @export
plotAdmixtureComparison <- function(comp.df) {
  breaks <- unique(as.character(comp.df$admix_cat))
  breaks_assay <- c(paste0(breaks,' ST'), paste0(breaks,' scRNA'))

  comp.df$admix_cat2 <- paste0(comp.df$admix_cat,' ',comp.df$assay)
  comp.df$admix_cat2 <- factor(comp.df$admix_cat2,levels=breaks_assay)
  comp.df$assay <- factor(comp.df$assay,levels=c('ST','scRNA'))

  myColors <- c(RColorBrewer::brewer.pal(3,'Reds'),RColorBrewer::brewer.pal(3,'Blues'))

  p <- comp.df %>%
    ggplot(
      aes(x = assay, y = frac, fill = admix_cat2)
    ) +
    geom_col() +
    ylab('Fraction of cells') +
    theme_classic(base_line_size = gg_line_thickness) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.x = element_blank(),
          panel.spacing = unit(0.1, "lines")) +
    geom_hline(yintercept=0,linewidth=gg_line_thickness,alpha=gg_line_alpha) +
    scale_fill_manual(breaks = breaks_assay, values = myColors) +
    scale_y_continuous(expand = expansion(0)) +
    scale_x_discrete(drop = FALSE) +
    facet_wrap(~cell_type, scales = 'free_x', nrow = 1,strip.position = 'bottom') +
    labs(fill='Admixture\ncategory') +
    p_theme +
    theme(strip.text.x = element_text(angle = 50, size=6),
          strip.background = element_blank(),
          strip.clip = "off",
          legend.text=element_text(size=6),
          legend.key.size = unit(6, 'pt'),
          legend.key.height = unit(6, 'pt'),
          legend.key.width = unit(6, 'pt')
    )
  return(p)
}

#' @export
plotAdmixtureCorrected <- function(so_spat_fracs,sc_fracs,source_ct) {
  # for each admix marker get the average frac in spatial as well as sc
  # then, we multiply the frac in spatial by (frac_sc/frac_spatial) from source cells only
  all_ct <- unique(so_spat@meta.data$cell_type)
  admix_fracs_spat_new <- list()
  admix_fracs_sc <- list()
  for (mark in admix_markers) {
    spat_cell_mal <- rownames(so_spat@meta.data)[so_spat@meta.data$cell_type==source_ct]
    sc_cell_mal <- rownames(sc_obj@meta.data)[sc_obj@meta.data$cell_type==source_ct]

    fracs_spat <- so_spat_fracs[mark,spat_cell_mal]
    fracs_sc <- sc_fracs[mark,sc_cell_mal]
    fracs_ratio <- mean(fracs_sc)/mean(fracs_spat)

    # now loop through each cell type
    for (ct in all_ct) {
      spat_cell_ct <- rownames(so_spat@meta.data)[so_spat@meta.data$cell_type==ct]
      sc_cell_ct <- rownames(sc_obj@meta.data)[sc_obj@meta.data$cell_type==ct]

      fracs_spat <- so_spat_fracs[mark,spat_cell_ct]
      fracs_sc <- sc_fracs[mark,sc_cell_ct]

      # now reduce spatial fracs by the above ratio
      fracs_spat2 <- fracs_spat*fracs_ratio

      if (!(ct %in% names(admix_fracs_spat_new))) {
        admix_fracs_spat_new[[ct]] <- list()
        admix_fracs_sc[[ct]] <- list()
      }
      admix_fracs_spat_new[[ct]][[mark]] <- fracs_spat2
      admix_fracs_sc[[ct]][[mark]] <- fracs_sc
    }
  }


  spat_total_fracs <- sapply(all_ct,function(ct){
    admix_fracs_spat_new_ct <- admix_fracs_spat_new[[ct]]
    spat_all_mark <- do.call(cbind,admix_fracs_spat_new_ct)
    spat_cmeans <- apply(spat_all_mark,MARGIN = 2,mean) # mean across cells, done for each gene
    spat_total <- sum(spat_cmeans)
    return(spat_total)
  })

  sc_total_fracs <- sapply(all_ct,function(ct){
    admix_fracs_sc_ct <- admix_fracs_sc[[ct]]
    sc_all_mark <- do.call(cbind,admix_fracs_sc_ct)
    sc_cmeans <- apply(sc_all_mark,MARGIN = 2,mean) # mean across cells, done for each gene
    sc_total <- sum(sc_cmeans)
    return(sc_total)
  })

  tmp <- cbind.data.frame(c(spat_total_fracs,sc_total_fracs),
                          c(rep('ST',length(spat_total_fracs)),rep('scRNA',length(sc_total_fracs))),
                          c(names(spat_total_fracs),names(sc_total_fracs)))
  colnames(tmp) <- c('total_av_frac','assay','cell_type')

  tmp_sub <- tmp[tmp$assay=='ST',]
  ct_ord <- tmp_sub$cell_type[order(tmp_sub$total_av_frac,decreasing = TRUE)]
  tmp$cell_type <- factor(tmp$cell_type,levels=ct_ord)
  tmp$assay <- factor(tmp$assay,levels=c('ST','scRNA'))
  p <- ggplot(tmp,aes(x = cell_type, y = total_av_frac, fill = assay)) +
    geom_bar(position="dodge", stat="identity") +
    ylab('Summed av.\nexpression fraction') +
    p_theme +
    ggtitle(paste0(toupper(substr(source_ct, 1, 1)), substr(source_ct, 2, nchar(source_ct)),
                   ' marker expression adjusted for\nspatial-sc differences in ',source_ct,' cells')) +
    theme_classic(base_line_size = gg_line_thickness) +
    xlab('Cell type') +
    scale_fill_manual(values = c("ST" = "red", "scRNA" = "blue")) +
    theme(axis.text.x = element_text(size = 7, angle = 30, hjust = 1),
          axis.text.y = element_text(size = 7),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          legend.key.size = unit(0.4, "cm"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          plot.title = element_text(hjust = 0.5,size = 10))
  return(p)
}










