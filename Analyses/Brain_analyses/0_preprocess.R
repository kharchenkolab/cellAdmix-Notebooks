
library(dplyr)
library(readr)
library(reshape2)
library(data.table)
library(cowplot)
library(ggplot2)
library(dataorganizer)
library(Seurat)


########################################################
## Subset molecule data
########################################################

print('Subsetting molecule data')

## read in full dataset
dat <- DatasetPath('mouse_hypothalamus', 'merfish_barcodes.csv') %>%
  read.csv()

# save a smaller version of the data with just naive mice and relevant bregma slices
# Behavior=='Naive' & Animal_sex == 'Female' & Bregma %in% c('0.16', '0.21')
dat_sub <- dat[dat$Behavior=='Naive' & dat$Animal_sex == 'Female' & dat$Bregma %in% c('0.16', '0.21'),]
print(dim(dat))
print(dim(dat_sub))
print(unique(dat_sub$Animal_ID))
DatasetPath('mouse_hypothalamus', 'merfish_subset.rds') %>% write_rds(dat_sub, .)


########################################################
## Prepare spatial metadata
########################################################

print('Preparing spatial metadata')

# read in cell type metadata
spatial_dat <- DatasetPath('mouse_hypothalamus', 'Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv') %>%
  read.csv()

# split up expression from metadata
spatial_meta <- spatial_dat[,1:9]
rownames(spatial_meta) <- spatial_meta$Cell_ID
spatial_meta$Cell_ID <- NULL
rm(spatial_dat)
gc()

spatial_meta <- spatial_meta[spatial_meta$Behavior=='Naive' & spatial_meta$Animal_sex == 'Female' & spatial_meta$Bregma %in% c('0.16', '0.21'),]
spatial_meta$cell <- rownames(spatial_meta)
ndx_change <- which(colnames(spatial_meta)%in%c('Centroid_X','Centroid_Y'))
colnames(spatial_meta)[ndx_change] <- c('x','y')
ndx_change <- which(colnames(spatial_meta)%in%c('Cell_class'))
colnames(spatial_meta)[ndx_change] <- c('cell_type')

# remove ambiguous cells
spatial_meta <- spatial_meta[spatial_meta$cell_type!='Ambiguous',]

# combine subtypes of OD cells and endothelial cells
spatial_meta$cell_type <- sapply(spatial_meta$cell_type,function(x){
  if (x=='Endothelial 1' | x=='Endothelial 2' | x=='Endothelial 3' ) {
    return('Endothelial')
  } else if (x=='OD Immature 1' | x=='OD Immature 2') {
    return('OD Immature')
  } else if (x=='OD Mature 1' | x=='OD Mature 2' | x=='OD Mature 3' | x=='OD Mature 4') {
    return('OD Mature')
  } else {
    return(x)
  }
})

all_mice <- unique(spatial_meta$Animal_ID)
all_slices <- unique(spatial_meta$Bregma)
ct_all <- unique(spatial_meta$cell_type)
ct_all <- paste0('n_',ct_all)

# now compute knn matrix per animal per slice
ct_counts_all <- list()
for (mouse in all_mice) {
  for (sl in all_slices) {
    # subset the metadata
    meta_sub <- spatial_meta[spatial_meta$Animal_ID==mouse & spatial_meta$Bregma==sl,]

    xn = meta_sub %>% select(x, y) %>% FNN::get.knn(k = 8) %>% .$nn.index
    ct_dict = meta_sub %>% mutate(id = 1:n()) %>% {setNames(.$cell_type, .$id)}
    cell_dict = meta_sub %>% mutate(id = 1:n()) %>% {setNames(.$cell, .$id)}
    nmat = apply(xn, 2, function(i){ct_dict[i]})

    ct_counts = nmat %>%
      reshape2::melt() %>%
      select(id = Var1, ct = value) %>%
      mutate(
        cell = cell_dict[id]
      ) %>%
      select(-id) %>%
      as.data.table %>%
      mutate(ct = paste0('n_', ct)) %>%
      count(cell, ct, .drop = T) %>%
      dcast(cell ~ ct, value.var = 'n', fill = 0)

    # add column of 0's for cell types not present
    for (ct in ct_all) {
      if (!(ct %in% colnames(ct_counts))) {
        ct_counts <- cbind(ct_counts,rep(0,nrow(ct_counts)))
        colnames(ct_counts)[ncol(ct_counts)] <- ct
      }
    }

    ct_counts <- as.data.frame(ct_counts)
    rownames(ct_counts) <- ct_counts$cell
    ct_counts$cell <- NULL

    # reorder columns
    ct_counts <- ct_counts[,ct_all]

    # store knn counts matrix
    ct_counts_all[[length(ct_counts_all)+1]] <- ct_counts
  }
}
ct_counts_all <- do.call(rbind.data.frame,ct_counts_all)

x = ct_counts_all %>% as.matrix
hc = hclust(dist(x,method = 'euclidean'),method = 'ward.D')
clust_dict = cutree(hc, 3)
clust_dict <- clust_dict[rownames(spatial_meta)]

# append clusters to metadata
spatial_meta$spat_clust <- clust_dict
spatial_meta$spat_clust <- as.factor(spatial_meta$spat_clust)

# save metadata
DatasetPath('mouse_hypothalamus', 'meta_with_spatial_clusts.rds') %>%
  write_rds(spatial_meta, .)


########################################################
## Plot results
########################################################

print('Plotting results')

# plot results to see if it worked as expected
all_plots_ct <- list()
all_plots_clust <- list()
for (mouse in all_mice) {
  for (sl in all_slices) {
    meta_sub <- spatial_meta[spatial_meta$Animal_ID==mouse & spatial_meta$Bregma==sl,]
    p <- ggplot(meta_sub,aes(x=x,y=y,color=cell_type)) +
      geom_point() +
      ggtitle(paste0(as.character(mouse),' ',as.character(sl))) +
      theme(plot.title = element_text(hjust = 0.5))
    all_plots_ct[[length(all_plots_ct)+1]] <- p
    p <- ggplot(meta_sub,aes(x=x,y=y,color=spat_clust)) +
      geom_point() +
      ggtitle(paste0(as.character(mouse),' ',as.character(sl))) +
      theme(plot.title = element_text(hjust = 0.5))
    all_plots_clust[[length(all_plots_clust)+1]] <- p
  }
}

f1 <- plot_grid(plotlist = all_plots_ct)
f2 <- plot_grid(plotlist = all_plots_clust)
f <- plot_grid(f1,f2,nrow=1)
f
f1
f2

# save these plots
pdf(OutputPath('hypo.pdf'), useDingbats=FALSE, width=30, height=8)
f
dev.off()


########################################################
## Process scRNA data
########################################################

print('Processing scRNA data')

DatasetPath('mouse_hypothalamus_rna', 'scRNA_seurat.rds') %>%
  process_brain_scrna()