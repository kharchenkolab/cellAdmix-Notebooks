

# to normalize scores by total distance (number of grid points) between all pairs
compute_cell_score <- function(coords_1, coords_2, im, params) {
  coords_1 <- floor(coords_1)
  coords_2 <- floor(coords_2)

  rownames(coords_1) <- paste0('mol',1:nrow(coords_1))
  rownames(coords_2) <- paste0('mol',1:nrow(coords_2))

  ## getting the combinations of molecules to test
  # need to do this within z-slices
  all_z <- as.numeric(unique(c(coords_1$z_index,coords_2$z_index)))
  all_comb_test <- list()
  for (myz in all_z) {
    coords_1_z <- coords_1[coords_1$z_index==myz,]
    coords_2_z <- coords_2[coords_2$z_index==myz,]

    if (params$compare_type=='same_groups') {
      if (nrow(coords_1_z)<2) {
        next
      }
      mol_combs_z <- as.data.frame(t(combn(rownames(coords_1_z),m = 2,simplify = TRUE)))
    } else {
      mol_combs_z <- expand.grid(rownames(coords_1_z), rownames(coords_2_z),
                                 stringsAsFactors = FALSE)
    }

    ### downsampling number of molecule pairs if specified
    if (!is.null(params$downsamp_pairs_to)) {
      if (nrow(mol_combs_z) > params$downsamp_pairs_to) {
        rkeep <- sample(1:nrow(mol_combs_z),params$downsamp_pairs_to)
        mol_combs_z <- mol_combs_z[rkeep,]
      }
    }
    all_comb_test[[length(all_comb_test)+1]] <- mol_combs_z
  }
  all_comb_test <- do.call(rbind,all_comb_test)

  if (is.null(all_comb_test)) {
    return(NA)
  }
  if (nrow(all_comb_test)<params$min_mols) {
    return(NA)
  }

  scores <- c()
  all_dists <- c()
  for (rndx in 1:nrow(all_comb_test)) {
    i <- all_comb_test[rndx,1]
    j <- all_comb_test[rndx,2]
    if (params$same_z && (coords_1[i, 3] != coords_2[j, 3])) next
    if (all(coords_1[i, ] == coords_2[j, ])) next
    score <- cellAdmix:::compute_summary(coords_1[i, ], coords_2[j, ], im)
    if (params$norm_scores) {
      pt_dist <- sqrt(sum((coords_1[i, ] - coords_2[j, ])^2))
      score <- score/pt_dist
      all_dists <- c(all_dists,pt_dist)
    } else if (params$balance_dists) {
      pt_dist <- sqrt(sum((coords_1[i, ] - coords_2[j, ])^2))
      all_dists <- c(all_dists,pt_dist)
    }
    scores <- c(scores, score)
  }

  if (length(scores)==0) { # all pairs were the same molecules
    return(NA)
  }

  if (params$balance_dists) {
    return(list(scores,all_dists))
  } else {
    if (params$sep_type == "pairwise_avg") {
      return(mean(scores))
    } else if (params$sep_type == "pairwise_max") {
      return(max(scores))
    } else {
      stop("Separation type not implemented")
    }
  }
}


# 'compare_type' can be either 'cross_groups' (i.e., admix-pure vs pure-pure)
# or it can be 'same_groups' (i.e.,admix-admix vs pure-pure)
compute_membrane_sep_per_cell <- function(df, im, fov,
                                          params,
                                          ct_test=NULL,n.cores=1) {
  # subset to only the focal cell type
  if (!is.null(ct_test)) {
    df_ct <- df[df$celltype==ct_test,]
  } else {
    df_ct <- df
  }

  # subset to only the fov we loaded images for
  df_ct <- df_ct[df_ct$fov==fov,]

  # checking if there are no cells with admixture labels or if all molecules are labeled admixture
  has_nat_admix <- length(unique(df_ct$is_admixture))==2
  if (!has_nat_admix) {
    print("No valid cells found")
    return(c())
  }

  # prefilter cells by whether they have sufficient molecules, so I can downsample cells tested per fov if desires
  f_counts_all <- table(df_ct[,c('cell','is_admixture')])
  cells <- names(which(rowSums(f_counts_all<params$min_mols)==0))

  # if there are no cells, return empty vector
  if (length(cells)==0) {
    print("No valid cells found")
    return(c())
  }

  # downsample cells if specified
  if (!is.null(params$downsamp_cells_to)) {
    if (length(cells)>params$downsamp_cells_to) {
      cells <- sample(cells,params$downsamp_cells_to)
    }
  }

  results <- plapply(1:length(cells),function(cur_cell_ndx) {
    print(cur_cell_ndx)
    cur_cell <- cells[cur_cell_ndx]
    cur_df <- df_ct %>% filter(cell == cur_cell)
    mols_admix_coords <- cur_df %>% filter(is_admixture) %>% select(x, y, z_index)
    mols_pure_coords <- cur_df %>% filter(!is_admixture) %>% select(x, y, z_index)

    admix_frac <- mean(cur_df$is_admixture)

    ## if balancing distance distributions...
    if (params$balance_dists) {
      if (params$compare_type=='cross_groups') {
        score_admix <- compute_cell_score(mols_admix_coords, mols_pure_coords, im, params)
        score_other <- compute_cell_score(mols_pure_coords, mols_pure_coords, im, params)
      } else if (params$compare_type=='same_groups') {
        score_admix <- compute_cell_score(mols_admix_coords, mols_admix_coords, im, params)
        score_other <- compute_cell_score(mols_pure_coords, mols_pure_coords, im, params)
      } else {
        stop("The compare_type parameter must be one of either 'cross_groups' or 'same_groups'")
      }

      if (is.na(score_admix[[1]][1]) || is.na(score_other[[1]][1])) { # not enough molecules weren't on same z levels
        return(NA)
      }

      dists_admix <- score_admix[[2]]
      dists_other <- score_other[[2]]
      score_admix <- score_admix[[1]]
      score_other <- score_other[[1]]
      mset_sizes <- c(nrow(mols_admix_coords),nrow(mols_pure_coords))
      ndx_smaller <- which(mset_sizes==min(mset_sizes))[1]
      dist_range <- range(c(dists_other,dists_admix))
      dist_buffer <- .02*abs(dist_range[2]-dist_range[1])
      ndx_match <- c()
      if (ndx_smaller==1) {
        for (mydist in dists_admix) {
          mymatch <- which(abs(dists_other-mydist)<dist_buffer)[1]
          if (is.na(mymatch)) {
            ndx_match <- c(ndx_match,NA)
          } else {
            ndx_match <- c(ndx_match,mymatch)
            dists_other[mymatch] <- Inf
          }
        }
        # remove any NAs (no match)
        ndx_keep <- which(!is.na(ndx_match))
        score_admix <- score_admix[ndx_keep]
        ndx_match <- ndx_match[ndx_keep]
        score_other <- score_other[ndx_match]
      } else {
        for (mydist in dists_other) {
          mymatch <- which(abs(dists_admix-mydist)<dist_buffer)[1]
          if (length(mymatch)==0) {
            ndx_match <- c(ndx_match,NA)
          } else {
            ndx_match <- c(ndx_match,mymatch)
            dists_admix[mymatch] <- Inf
          }
        }
        # remove any NAs (no match)
        ndx_keep <- which(!is.na(ndx_match))
        score_other <- score_other[ndx_keep]
        ndx_match <- ndx_match[ndx_keep]
        score_admix <- score_admix[ndx_match]
      }

      # if there aren't enough pair scores with matching distances, return NA
      if (length(score_other)<params$min_mols | length(score_admix)<params$min_mols) {
        return(NA)
      }

      if (params$sep_type == "pairwise_avg") {
        score_other <- mean(score_other)
        score_admix <- mean(score_admix)
      } else if (params$sep_type == "pairwise_max") {
        score_other <- max(score_other)
        score_admix <- max(score_admix)
      } else {
        stop("Separation type not implemented")
      }

    } else { # if not balancing dists
      if (params$compare_type=='cross_groups') {
        score_admix <- compute_cell_score(mols_admix_coords, mols_pure_coords, im, params)
        score_other <- compute_cell_score(mols_pure_coords, mols_pure_coords, im, params)
      } else if (params$compare_type=='same_groups') {
        score_admix <- compute_cell_score(mols_admix_coords, mols_admix_coords, im, params)
        score_other <- compute_cell_score(mols_pure_coords, mols_pure_coords, im, params)
      } else {
        stop("The compare_type parameter must be one of either 'cross_groups' or 'same_groups'")
      }

      if (is.na(score_admix) || is.na(score_other)) { # molecules weren't on same z levels
        return(NA)
      }
    }

    return(list(admix_frac = admix_frac,
                score_admix = score_admix,
                score_other = score_other))
  },mc.preschedule=FALSE, n.cores=n.cores, progress=TRUE, fail.on.error=TRUE)
  results <- results[!is.na(results)]

  if (length(results) == 0) {
    print("No valid cells found")
  }

  return(results)
}


get_scores_all_fov <- function(df,score_params,img_load_fn,ct_test=NULL,n.cores=1,par_inner=FALSE) {
  all_fov <- unique(df$fov)
  all_fov <- all_fov[order(all_fov,decreasing=FALSE)]
  cores_saved <- 1
  if (par_inner) {
    cores_saved <- n.cores
    n.cores <- 1
  }
  ratio_res <- plapply(all_fov,function(fov) {
    # load images for the fov
    im = img_load_fn(fov=fov, normalize=TRUE)

    # calculate membrane scores in fibroblast cells
    scores <- compute_membrane_sep_per_cell(
      df, im, fov, score_params, ct_test=ct_test, n.cores=cores_saved
    )

    if (length(scores)==0) {
      return(c())
    } else {
      # calculate ratio of admixture-pure pairs to pure-pure molecule pairs
      all_ratios <- sapply(1:length(scores),function(i){
        x <- scores[[i]]
        ad_score <- x[[2]]
        pure_score <- x[[3]]
        ratio <- ad_score / pure_score
        return(ratio)
      })
      return(all_ratios)
    }
  },mc.preschedule=FALSE, n.cores=n.cores, progress=TRUE, fail.on.error=TRUE)

  ratio_res_un <- unlist(ratio_res)
  return(ratio_res_un)
}

plot_scores <- function(scores,title,score_params,log_transform=TRUE,
                        binwidth=NULL,capval=NULL) {
  if (log_transform) {
    scores <- log(scores)
    myintercept <- 0
  } else {
    myintercept <- 1
  }

  if (is.null(binwidth)) {
    binwidth <- (2 * IQR(scores)) / (length(scores)^(1/3))
  }

  scores_df <- as.data.frame(scores)
  colnames(scores_df) <- 'sep_ratio'

  if (log_transform) {
    if (score_params$compare_type=='cross_groups') {
      xlabel <- 'Log membrane separation ratio\n(pure vs admix / pure vs pure)'
    } else {
      xlabel <- 'Log membrane separation ratio\n(admix vs admix / pure vs pure)'
    }
  } else {
    if (score_params$compare_type=='cross_groups') {
      xlabel <- 'Membrane separation ratio\n(native vs admix / native vs native)'
    } else {
      xlabel <- 'Membrane separation ratio\n(admix vs admix / native vs native)'
    }
  }

  scores_df$capped_value <- scores_df$sep_ratio
  if (!is.null(capval)) {
    scores_df$capped_value[scores_df$capped_value>capval] <- capval
  }

  if (!is.null(capval)) {
    p <- ggplot(scores_df, aes(x = capped_value)) +
      # Create histogram with customized bin width, fill, and border color
      geom_histogram(binwidth = binwidth,    # adjust binwidth to best display your data
                     color = "black",
                     fill = "white",
                     alpha = 0.8,
                     size=.2) +
      geom_vline(xintercept = myintercept,color='red',linetype='dashed',alpha=.9,linewidth=.65) +
      scale_x_continuous(
        breaks = seq(floor(min(scores_df$capped_value)), capval, by = .5),
        labels = function(x) ifelse(x == capval, paste0(">", capval), x),
        expand = c(0, 0),
        guide = guide_axis(check.overlap = TRUE)
      ) +
      labs(title = title,
           x = xlabel, y = 'Number of cells') +
      theme_classic(base_line_size = .23) +
      # Fine-tune theme elements for clarity and aesthetics
      theme(
        text = element_text(family = "Arial"),
        plot.title   = element_text(size = 8, hjust = 0.5),
        axis.title   = element_text(size = 7),
        axis.text    = element_text(size = 6, color = "black"),
        legend.position = "none",  # remove legend if not necessary
        plot.margin  = margin(10, 10, 10, 10)
      )
  } else {
    p <- ggplot(scores_df, aes(x = sep_ratio)) +
      # Create histogram with customized bin width, fill, and border color
      geom_histogram(binwidth = binwidth,    # adjust binwidth to best display your data
                     color = "black",
                     fill = "white",
                     alpha = 0.8,
                     size=.2) +
      geom_vline(xintercept = myintercept,color='red',linetype='dashed',alpha=.9,linewidth=.65) +
      labs(title = title,
           x = xlabel, y = 'Number of cells') +
      theme_classic(base_line_size = .23) +
      # Fine-tune theme elements for clarity and aesthetics
      theme(
        text = element_text(family = "Arial"),
        plot.title   = element_text(size = 8, hjust = 0.5),
        axis.title   = element_text(size = 7),
        axis.text    = element_text(size = 6, color = "black"),
        legend.position = "none",  # remove legend if not necessary
        plot.margin  = margin(10, 10, 10, 10)
      )
  }

  return(p)
}
