#' @import dataorganizer
NULL

###################
## NSCLC dataset ##
###################

#' @export
prepare_nsclc_metadata <- function(gem=NULL, ct_annot=NULL, reps=c('one','all')) {
  reps <- match.arg(reps)
  if (is.null(gem)) {
    load(DatasetPath('nsclc', 'SMI_Giotto_Object.RData'))
  }

  if (is.null(ct_annot)) {
    ct_annot <- DatasetPath('nsclc_lung5_1', 'annotation_adj.csv') %>% read.csv()
  }

  # extract relevant metadata from giotto object
  cell_meta <- gem@cell_metadata$rna
  cell_locs <- gem@spatial_locs$raw
  cell_meta <- cbind.data.frame(cell_locs[,c(1:2)],cell_meta)

  # subset to a single donor/replicate
  if (reps=='one') {
    cell_meta <- cell_meta[cell_meta$Run_Tissue_name=='Lung5_Rep1',]
  } else if (reps=='all') {
    cell_meta <- cell_meta[cell_meta$Run_Tissue_name %in% c('Lung5_Rep1','Lung5_Rep2','Lung5_Rep3'),]
  }

  colnames(cell_meta)[3] <- 'cell'
  colnames(cell_meta)[34] <- 'celltype'

  ## modify celltype names
  if (reps=='one' && !(length(ct_annot) == 1 && ct_annot == FALSE)) {
    cell_meta$celltype <- ct_annot$cell_type
  }

  # rename tumor cells
  cell_meta$celltype %<>% {ifelse(startsWith(., 'tumor'), 'malignant', .)}

  # group all immune cells under one annotation for the visualization
  immune_cell_types <- c(
    'B-cell', 'NK', 'T CD4 memory', 'T CD4 naive', 'T CD8 memory', 'T CD8 naive', 'Treg',
    'plasmablast', 'mast', 'mDC', 'monocyte', 'pDC', 'neutrophil', "DC", "B cells", "CD4+ T cells", "CD8+ T cells"
  )

  cell_meta$cell_type_coarse <- cell_meta$celltype %>% {ifelse(. %in% immune_cell_types, 'immune other', .)}
  rownames(cell_meta) <- cell_meta$cell

  rownames(cell_meta) <- cell_meta$cell

  colnames(cell_meta)[which(colnames(cell_meta)=='sdimx')] <- 'x'
  colnames(cell_meta)[which(colnames(cell_meta)=='sdimy')] <- 'y'
  cell_meta$z <- 1

  return(cell_meta)
}

#' @export
prepare_nsclc_transcript_data <- function(cell_meta,tx_dat=NULL,reps=c('one','all')) {
  if (is.null(tx_dat)) {
    if (reps=='one') {
      tx_dat <- DatasetPath('nsclc_lung5_1', 'Lung5_Rep1_tx_file.csv') %>%
        readr::read_csv(show_col_types=FALSE) %>% as.data.frame()
      tx_dat$cell <- paste0('c_1_',tx_dat$fov,'_',tx_dat$cell_ID)

    } else if (reps=='all') {
      tx_dat <- DatasetPath('nsclc', 'tx_dat_all_reps.rds') %>% readRDS()
    }
  }

  # subsetting data to same cells we have annotations for
  df <- tx_dat[tx_dat$cell %in% cell_meta$cell,]

  # append cell type annotations to molecule-level data
  match_ndx <- match(df$cell,cell_meta$cell)
  df$celltype <- cell_meta$celltype[match_ndx]
  df$cell_type <- df$celltype

  # Change x and y coordinate column names
  colnames(df)[3:4] <- c('x','y')

  # Change gene column name
  colnames(df)[8] <- c('gene')

  # remove negative probes
  negprb_ndx <- grep("^NegPrb", df$gene)
  df <- df[-negprb_ndx,]

  # adding a column for molecule ID
  df$mol_id <- as.character(1:nrow(df))
  rownames(df) <- df$mol_id

  # converting coordinate units from pixels to microns based on information from their publication
  df$x <- (df$x * 180) / 1000
  df$y <- (df$y * 180) / 1000

  # assigning z-coordinates to an expected height in microns
  df$z <- (df$z * 800) / 1000

  return(df)
}

#' @export
prepare_nsclc_scrna <- function(so_rna=NULL) {
  if (is.null(so_rna)) {
    so_rna <- DatasetPath('nsclc_rna', 'so.rds') %>% readRDS()
  }
  so_rna$cell_type[endsWith(so_rna$cell_type, 'specific')] <- 'malignant'

  cell_types_keep <- c(
    'malignant','fibroblast','macrophage','monocyte','DC','mast',
    'Treg','CD4+ T cells','CD8+ T cells',
    'NK','B cells','plasmablast','neutrophil',
    'endothelial','epithelial'
  )

  cells_keep <- colnames(so_rna)[so_rna$cell_type %in% cell_types_keep]
  so_rna <- subset(so_rna, cells=cells_keep)

  # reordering cell types to look nicer on the dotplots
  so_rna$cell_type <- factor(so_rna$cell_type, levels=cell_types_keep)
  Idents(so_rna) <- so_rna$cell_type
  return(so_rna)
}


#' @export
load_markers_nsclc <- function(filter_markers=TRUE) {
  markers = DatasetPath('nsclc_rna', 'nsclc_markers_list_full.tsv') %>%
    read.table(sep = '\t',header=TRUE)

  markers <- markers[,c('Gene','cluster')]
  colnames(markers) <- c('gene','marker_of')

  g_rem <- markers$gene[duplicated(markers$gene)]
  markers <- markers[!(markers$gene %in% g_rem),]

  if (filter_markers) {
    # filtering out cell types with very few marker genes
    marker_counts <- table(markers$marker_of)
    ct_keep <- names(marker_counts)[marker_counts>2]
    markers <- markers[markers$marker_of %in% ct_keep,]
  }

  return(markers)
}


#' @export
load_images_nsclc <- function(fov, normalize=TRUE) {
  # List all file names in the directory matching the FOV pattern
  load_base <- DatasetPath('nsclc_lung5_1', 'RawMorphology') # for bayes
  fnames <- list.files(load_base, pattern=paste0("_F", sprintf("%03d", fov), ".*\\.TIF$"))

  # Extract Z-axes
  zvals <- sapply(fnames, function(fname) {
    # Split the filename by "_" and find the Z part (ends with .TIF)
    parts <- strsplit(fname, "_")[[1]]
    z_part <- parts[length(parts)]
    # Extract the numeric part from Z (e.g., "Z10.TIF" -> "10")
    as.numeric(gsub("\\D", "", z_part))  # Replace non-digits with an empty string
  })

  # Sort files by Z-values
  idx_srt <- order(zvals)
  fnames_srt <- fnames[idx_srt]

  # Load images in sorted order
  im_list <- lapply(fnames_srt, function(fname) {
    cur_im <- file.path(load_base, fname) %>% tiff::readTIFF(as.is = TRUE,all=TRUE)
    cur_im <- cur_im[[1]] # Membrane channel is the first

    # reverse the y-coordinates
    cur_im <- cur_im[nrow(cur_im):1,]

    return(cur_im)
  })
  im <- abind::abind(im_list, along = 3)
  im <- aperm(im, perm = c(3, 1, 2))


  # Normalize images if requested
  if (normalize) {
    im <- normalize_images(im)
  }

  return(im)
}

#' @export
load_scrnaseq_matrix_nsclc <- function() {
  gn_nm <- DatasetPath('nsclc_rna', 'GSE127465_gene_names_human_41861.tsv') %>% read.table(sep="\t")

  sc_meta <- DatasetPath('nsclc_rna', 'GSE127465_human_cell_metadata_54773x25.tsv') %>%
    read.table(sep="\t", header=TRUE)

  sc_meta$cell_id <- paste0(sc_meta$Library,'_',sc_meta$Barcode)
  rownames(sc_meta) <- sc_meta$cell_id

  # dims are cells x genes
  sc_dat <- DatasetPath('nsclc_rna', 'GSE127465_human_counts_normalized_54773x41861.mtx') %>% Matrix::readMM()
  sc_dat <- as.matrix(sc_dat)
  sc_dat <- t(sc_dat)
  rownames(sc_dat) <- gn_nm[,1]
  colnames(sc_dat) <- sc_meta$cell_id

  # read in the integrated cell type annotations for the data
  ct_annot <- DatasetPath('nsclc_rna', 'annotation_adj.csv') %>% read.csv()

  sc_meta$cell_type <- ct_annot[,2]

  # rename patient specific cells to tumor
  ndx_change <- which(sc_meta$cell_type %in% c(
    'Patient1-specific','Patient2-specific', 'Patient3-specific','Patient4-specific',
    'Patient5-specific','Patient6-specific', 'Patient7-specific'
  ))

  sc_meta$cell_type[ndx_change] <- 'malignant'

  sc_obj <- CreateSeuratObject(sc_dat, meta.data=sc_meta)

  # only keeping cell types that are also annotated in the spatial data
  ct_annot_spatial <- DatasetPath('nsclc_lung5_1', 'annotation_adj.csv') %>% read.csv()
  ct_keep <- c(unique(ct_annot_spatial$cell_type),'malignant')
  cells_keep <- colnames(sc_dat)[sc_meta$cell_type %in% ct_keep]
  sc_obj <- subset(sc_obj,cells=cells_keep)

  ct_levels <- c(
    'macrophage','monocyte','DC','mast',
    'Treg','CD4+ T cells','CD8+ T cells',
    'NK','B cells','plasmablast','neutrophil',
    'fibroblast','endothelial','epithelial','malignant'
  )

  sc_obj$cell_type <- factor(sc_obj$cell_type, levels=ct_levels)
  Idents(sc_obj) <- sc_obj$cell_type

  return(sc_obj)
}

#################
## Gut dataset ##
#################

#' @export
prepare_gut_tx_and_meta <- function(
    rm_small_cells=TRUE, append_regions=FALSE, min_assignment_confidence=0.525
  ) {
  df <- DatasetPath('mouse_gut', 'segmentation', 'segmentation.csv') %>%
    data.table::fread()
  cell_annot <- DatasetPath('mouse_gut', 'clustering', 'cell_assignment.csv') %>%
    data.table::fread() %>% as.data.frame() %>% set_colnames(c('cell', 'celltype'))

  df <- df %>% filter(cell != 0)

  rownames(cell_annot) <- paste0('cell_', cell_annot$cell)

  # condense cell type annotations by coarse annotations
  cell_annot$cell_type_coarse <- sapply(cell_annot$celltype,function(x){
    if (x=='B (Follicular, Circulating)' | x=='B (Plasma)' | x=='T (CD8+)' | x=='T (CD4+)' | x=='Macrophage + DC') {
      return('Immune')
    } else if (x=='Enterocyte (Top Villus)' | x=='Enterocyte (Bottom Villus)' | x=='Enterocyte (Mid Villus)') {
      return('Enterocyte')
    } else if (x=='ICC' | x=='Myenteric Plexus' | x=='Pericyte' | x=='Telocyte' | x=='Tuft' | x=='Paneth' | x=='Endothelial') {
      return('Other')
    } else {
      return(x)
    }
  })

  cell_annot$celltype <- sapply(cell_annot$celltype,function(x){
    if (x == 'B (Follicular, Circulating)' | x == 'B (Plasma)') {
      return('B')
    } else if (x == 'Enterocyte (Bottom Villus)' | x == 'Enterocyte (Mid Villus)' | x == 'Enterocyte (Top Villus)') {
      return('Enterocyte')
    } else if (x == 'T (CD4+)' | x == 'T (CD8+)') {
      return('T')
    } else {
      return(x)
    }
  })

  cell_annot <- cell_annot[cell_annot$celltype!='Removed',]

  # make spatial data into seurat object
  df <- df[df$assignment_confidence>=min_assignment_confidence,]
  df$cell <- paste0('cell_',df$cell)
  df <- df[df$cell %in% rownames(cell_annot),]
  df$celltype <- cell_annot[df$cell,'celltype']

  cell_annot <- cell_annot[rownames(cell_annot) %in% unique(df$cell),]
  cell_annot$cell <- rownames(cell_annot)

  df <- as.data.frame(df)

  ## append gut region annotations
  if (append_regions) {
    gut_regions <- DatasetPath('mouse_gut', 'gut_regions.rds') %>% readRDS()
    gut_regions <- gut_regions[cell_annot$cell]
    cell_annot$niche <- as.factor(gut_regions)
  }

  # removing cells with very few molecules
  if (rm_small_cells) {
    cell_counts <- table(df$cell)
    cells_keep <- names(cell_counts)[cell_counts>=10]
    df <- df[df$cell %in% cells_keep,]
    cell_annot <- cell_annot[rownames(cell_annot) %in% cells_keep,]
  }

  df$mol_id <- 1:nrow(df)
  rownames(df) <- df$mol_id

  # need to add 1 to x,y coordinates because they're 0 indexed
  df$x <- df$x + 1
  df$y <- df$y + 1

  centroids <- df %>%
    group_by(cell) %>%
    summarise(
      x = mean(x, na.rm = TRUE),
      y = mean(y, na.rm = TRUE)
    )
  centroids <- as.data.frame(centroids)
  rownames(centroids) <- centroids$cell
  cell_annot <- cbind.data.frame(cell_annot,centroids[rownames(cell_annot),c('x','y')])
  cell_annot$z <- 1

  return(list(cell_meta=cell_annot, df=df))
}

#' @export
adjust_gut_cell_types <- function(cell.types) {
  ct_rename = c(
    'B (Follicular, Circulating)' = 'B', 'B (Plasma)' = 'B',
    'Myenteric Plexus' = 'Other', 'ICC' = 'Other',
    'Pericyte' = 'Stromal', 'Telocyte' = 'Stromal',
    'Enterocyte (Bottom Villus)' = 'Enterocyte', 'Enterocyte (Mid Villus)' = 'Enterocyte',
    'Enterocyte (Top Villus)' = 'Enterocyte',
    'T (CD8+)' = 'T', 'T (CD4+)' = 'T'
  )

  cell.types %<>% {ifelse(. %in% names(ct_rename), ct_rename[.], .)} %>%
    setNames(names(cell.types))

  return(cell.types)
}

#' @export
prepare_gut_scrna <- function() {
  # load the scRNA data
  so_rna <- DatasetPath('mouse_gut_rna', 'so.rds') %>% readRDS()
  so_rna <- so_rna %>% subset(cells=colnames(.)[!(.$cell_type %in% c('EP', 'Enteroendocrine'))])
  so_rna$cell_type %<>% adjust_gut_cell_types()
  Idents(so_rna) <- factor(so_rna$cell_type,levels=c('Goblet','Stem + TA','Enterocyte',
                                                     'Paneth','Tuft'))
  return(so_rna)
}

#' @export
load_markers_gut <- function(filter_markers=TRUE) {
  # load markers
  markers <- DatasetPath('mouse_gut_rna', 'aviv_gut_markers_sub2.csv') %>% read.csv()
  markers <- markers[,c('Gene','cell_type')]
  colnames(markers) <- c('gene','marker_of')

  g_rem <- markers$gene[duplicated(markers$gene)]
  markers <- markers[!(markers$gene %in% g_rem),]

  if (filter_markers) {
    # filtering out cell types with very few marker genes
    marker_counts <- table(markers$marker_of)
    ct_keep <- names(marker_counts)[marker_counts>2]
    markers <- markers[markers$marker_of %in% ct_keep,]
  }
  return(markers)
}


load_images_gut <- function(fov,normalize) {
  file_nm <- DatasetPath('mouse_gut', 'raw_data', 'membrane_stack.tif') # for bayes

  im_list <- tiff::readTIFF(file_nm,as.is = TRUE,all=TRUE)

  # combine z stacks
  im <- abind::abind(im_list, along = 3)

  # put z-dimension first
  im <- aperm(im, perm = c(3, 1, 2))

  # Normalize images if requested
  if (normalize) {
    im <- normalize_images(im)
  }

  return(im)
}


#########################
## Mouse hypothalamus ##
#########################

#' @export
prepare_brain_tx_and_meta <- function(rm_small_cells=TRUE) {

  ## read in molecule level data
  dat_sub <- DatasetPath('mouse_hypothalamus', 'merfish_subset.rds') %>% readRDS()
  spatial_genes <- unique(dat_sub$Gene_name)

  # load metadata
  spatial_meta <- DatasetPath('mouse_hypothalamus', 'meta_with_spatial_clusts.rds') %>% readRDS()

  # reduce data to same cells in metadata
  dat_sub <- dat_sub[dat_sub$Cell_name %in% unique(spatial_meta$cell),]
  colnames(dat_sub)[1:2] <- c('gene','cell')
  colnames(dat_sub)[7:9] <- c('x','y','z')

  # appending cell types to dat_sub
  spatial_meta$celltype <- spatial_meta$cell_type
  ndx_match <- match(dat_sub$cell,spatial_meta$cell)
  dat_sub$celltype <- spatial_meta[ndx_match,'celltype']

  # remove blank probes
  negprb_ndx <- grep("Blank", dat_sub$gene)
  dat_sub <- dat_sub[-negprb_ndx,]

  # removing Ependymal ct because there is only 1 cell
  ndx_rm <- which(dat_sub$celltype=='Ependymal')
  dat_sub <- dat_sub[-ndx_rm,]
  ndx_rm <- which(spatial_meta$celltype=='Ependymal')
  spatial_meta <- spatial_meta[-ndx_rm,]
  rownames(spatial_meta) <- spatial_meta$cell

  spatial_meta$z <- 1

  # removing cells with very few molecules
  if (rm_small_cells) {
    cell_counts <- table(dat_sub$cell)
    cells_keep <- names(cell_counts)[cell_counts>=10]
    dat_sub <- dat_sub[dat_sub$cell %in% cells_keep,]
    spatial_meta <- spatial_meta[rownames(spatial_meta) %in% cells_keep,]
  }

  dat_sub$mol_id <- as.character(1:nrow(dat_sub))
  rownames(dat_sub) <- dat_sub$mol_id

  return(list(spatial_meta,dat_sub))
}

process_brain_scrna <- function(dir_save=NULL) {
  ## load the scRNA data to get marker genes
  sc_dat <- DatasetPath('mouse_hypothalamus_rna') %>% Read10X()

  # read in scRNA-seq metadata (downloaded from publication supplements)
  sc_meta <- DatasetPath('mouse_hypothalamus_rna', 'scRNA_metadata.xlsx') %>% readxl::read_excel()
  sc_meta <- as.data.frame(sc_meta)
  colnames(sc_meta) <- sc_meta[1,]
  sc_meta <- sc_meta[2:nrow(sc_meta),]
  rownames(sc_meta) <- sc_meta[,'Cell name']

  # editing column names to comply with seurat
  colnames(sc_meta)[c(1,3,4:6)] <- c('Cell_name','Replicate_number','Cell_class','Non_neuronal_cluster','Neuronal_cluster')

  sc_seurat <- CreateSeuratObject(counts = sc_dat, min.cells = 0, meta.data = sc_meta)

  sc_seurat <- NormalizeData(sc_seurat)

  # removing cell types not in spatial
  cells_keep <- rownames(sc_seurat@meta.data)[!(sc_seurat@meta.data$Cell_class %in% c("Ambiguous","Unstable","Newly formed oligodendrocyte"))]
  sc_seurat_clean <- subset(sc_seurat,cells = cells_keep)
  Idents(sc_seurat_clean) <- sc_seurat_clean@meta.data$Cell_class
  cells_keep <- rownames(sc_seurat_clean@meta.data)[!(sc_seurat_clean@meta.data$Cell_class %in% c("Macrophage","Fibroblast"))] # these ctypes not in spatial
  sc_seurat_clean <- subset(sc_seurat_clean,cells = cells_keep)
  if (!is.null(dir_save)) {
    saveRDS(sc_seurat_clean,file=dir_save)
  }
  return(sc_seurat_clean)
}


#' @export
prepare_brain_scrna <- function(short_idents=TRUE) {
  sc_obj <- DatasetPath('mouse_hypothalamus_rna', 'scRNA_seurat.rds') %>% readRDS()
  # harmonize ct names between sc and spatial
  if (short_idents) {
    levels(Idents(sc_obj)) <- c("Excitatory","Inhibitory","OD Mature","Endothelial",
                                "OD Immature","Microglia","Mural","Astrocyte","Ependymal")
  }
  sc_obj$cell_type <- Idents(sc_obj)

  return(sc_obj)
}

load_markers_brain <- function(filter_markers=TRUE) {
  # read in marker genes
  markers <- CachePath('brain_markers_sub.csv') %>% read.csv()
  markers <- markers[,c('Gene','cell_type')]
  colnames(markers) <- c('gene','marker_of')

  markers$marker_of <- as.factor(markers$marker_of)
  levels(markers$marker_of) <- c("Astrocyte","Endothelial","Excitatory",
                                 "OD Immature","Inhibitory","OD Mature",
                                 "Microglia","Mural")
  markers$marker_of <- as.character(markers$marker_of)

  g_dup <- markers$gene[duplicated(markers$gene)]
  markers <- markers[!(markers$gene %in% g_dup),]

  # filtering out cell types with very few marker genes
  if (filter_markers) {
    marker_counts <- table(markers$marker_of)
    ct_keep <- names(marker_counts)[marker_counts>2]
    markers <- markers[markers$marker_of %in% ct_keep,]
  }
  return(markers)
}



## Breast cancer dataset

#' @export
prepare_BC_sc_spatial <- function(load_molecules=TRUE) {
  # # load scRNA data
  so_rna <- DatasetPath('human_breast_cancer_rna', 'so.rds') %>% readRDS()

  # load spatial in counts form
  so_spatial <- DatasetPath('human_breast_cancer', 'so.rds') %>% readRDS()


  so_spatial$cell_type[so_spatial$cell_type == 'Prolif_Invasive_Tumor'] <- 'Invasive_Tumor'
  so_rna$cell_type_adj[so_rna$cell_type_adj == 'Prolif_Invasive_Tumor'] <- 'Invasive_Tumor'
  so_rna$cell_type <- so_rna$cell_type_adj

  common_genes <- intersect(rownames(so_rna), rownames(so_spatial))

  so_spatial %<>% subset(features=common_genes)
  so_rna %<>% subset(features=common_genes)
  so_spatial %<>% subset(cells=colnames(.)[.$cell_type %in% unique(so_rna$cell_type)])

  so_rna %<>% subset(cells=colnames(.)[.$cell_type_adj != 'doublets'])

  so_spatial$cell_type[so_spatial$seurat_clusters == '7'] <- 'Plasma'
  so_spatial$cell_type[so_spatial$seurat_clusters == '8'] <- 'B_Cells'

  if (!load_molecules) {
    return(list(so_rna=so_rna, so_spatial=so_spatial))
  }

  # load spatial molecule data
  df_spatial <- DatasetPath('human_breast_cancer', 'transcripts.csv.gz') %>%
    data.table::fread() %>% as_tibble() %>%
    rename(x=x_location, y=y_location, z=z_location, gene=feature_name) %>%
    mutate(cell=as.character(cell_id), celltype=so_spatial$cell_type[cell]) %>%
    filter(!is.na(celltype), !grepl('Hybrid', celltype))

  df_spatial <- as.data.frame(df_spatial)
  df_spatial$cell <- paste0('cell_',df_spatial$cell)

  ## removing negative control probes
  negprb_ndx1 <- grep("Control", df_spatial$gene)
  negprb_ndx2 <- grep("antisense", df_spatial$gene)
  negprb_ndx3 <- grep("BLANK", df_spatial$gene)
  negprb_ndx <- unique(c(negprb_ndx1,negprb_ndx2,negprb_ndx3))
  df_spatial <- df_spatial[-negprb_ndx,]

  rownames(df_spatial) <- df_spatial$mol_id <- 1:nrow(df_spatial)

  # making cell_annot
  cell_annot <- unique(df_spatial[,c('cell','celltype')])
  rownames(cell_annot) <- cell_annot$cell
  centroids <- df_spatial %>%
    group_by(cell) %>%
    summarise(
      x = mean(x, na.rm = TRUE),
      y = mean(y, na.rm = TRUE)
    )
  centroids <- as.data.frame(centroids)
  rownames(centroids) <- centroids$cell
  cell_annot <- cbind.data.frame(cell_annot,centroids[rownames(cell_annot),c('x','y')])
  cell_annot$z <- 1

  df_spatial$transcript_id <- NULL

  return(list(df_spatial=df_spatial, cell_annot=cell_annot, so_rna=so_rna, so_spatial=so_spatial))
}

#' @export
load_markers_BC <- function(filter_markers=TRUE) {
  # load markers
  markers <- read.csv(file=CachePath('BC_markers_sub.csv'))
  markers <- markers[,c('Gene','cell_type')]
  colnames(markers) <- c('gene','marker_of')

  # remove duplicate genes if there are any
  g_rem <- markers$gene[duplicated(markers$gene)]
  markers <- markers[!(markers$gene %in% g_rem),]

  # filtering out cell types with very few marker genes
  if (filter_markers) {
    marker_counts <- table(markers$marker_of)
    ct_keep <- names(marker_counts)[marker_counts>2]
    markers <- markers[markers$marker_of %in% ct_keep,]
  }

  return(markers)
}



## Ovarian cancer dataset

prepare_OC_sc_spatial <- function(load_molecules=TRUE) {
  ### loading up scRNA
  so_rna <- DatasetPath('human_ovarian_cancer_rna', 'processed', 'so.rds') %>% readRDS()

  ### load up spatial data in seurat object
  so_spatial <- DatasetPath('human_ovarian_cancer', 'processed', 'so.rds') %>% readRDS()


  so_spatial$cell_type <- so_spatial$cell_type_full
  so_rna$cell_type <- so_rna$cell_type_full

  common_genes <- intersect(rownames(so_rna), rownames(so_spatial))

  so_spatial %<>% subset(features=common_genes)
  so_rna %<>% subset(features=common_genes)
  so_spatial %<>% subset(cells=colnames(.)[.$cell_type %in% unique(so_rna$cell_type)])

  if (!load_molecules) {
    return(list(so_rna=so_rna, so_spatial=so_spatial))
  }

  cell_annot <- so_spatial@meta.data
  ndx_change <- which(colnames(cell_annot) %in% c('x_centroid', 'y_centroid'))
  colnames(cell_annot)[ndx_change] <- c('x','y')
  cell_annot$celltype <- cell_annot$cell_type
  cell_annot$cell <- rownames(cell_annot)
  cell_annot$z <- 1

  # load spatial molecule data
  df_spatial <- DatasetPath('human_ovarian_cancer', 'transcripts.parquet') %>%
    arrow::read_parquet() %>% as_tibble() %>%
    rename(x=x_location, y=y_location, z=z_location, gene=feature_name) %>%
    mutate(cell=as.character(cell_id), cell_type=so_spatial$cell_type[cell]) %>%
    filter(!is.na(cell_type))

  df_spatial <- as.data.frame(df_spatial)

  # removing non-gene probes
  ndx_keep <- which(df_spatial$codeword_category %in% c('predesigned_gene','custom_gene'))
  df_spatial <- df_spatial[ndx_keep,]

  # removing cells with very few molecules
  cell_counts <- table(df_spatial$cell)
  cells_keep <- names(cell_counts)[cell_counts>=10]
  df_spatial <- df_spatial[df_spatial$cell %in% cells_keep,]
  cell_annot <- cell_annot[rownames(cell_annot) %in% cells_keep,]

  rownames(df_spatial) <- df_spatial$mol_id <- 1:nrow(df_spatial)
  df_spatial$celltype <- df_spatial$cell_type

  return(list(df_spatial,cell_annot,so_rna))
}



load_markers_OC <- function(filter_markers=TRUE) {
  # load markers
  markers <- DatasetPath('human_ovarian_cancer', 'processed', 'OC_markers_sub.csv') %>% read.csv()
  markers <- markers[,c('Gene','cell_type')]
  colnames(markers) <- c('gene','marker_of')

  g_rem <- markers$gene[duplicated(markers$gene)]
  markers <- markers[!(markers$gene %in% g_rem),]

  # filtering out cell types with very few marker genes
  if (filter_markers) {
    marker_counts <- table(markers$marker_of)
    ct_keep <- names(marker_counts)[marker_counts>2]
    markers <- markers[markers$marker_of %in% ct_keep,]
  }

  return(markers)
}


load_images_oc <- function(fov,normalize) {
  file_nm = DatasetPath('human_ovarian_cancer', 'membrane_lowres.tif')

  im <- EBImage::readImage(file_nm) %>% imageData() %>%
    t() %>% .[rev(seq_len(nrow(.))),]

  # reverse y coordinates
  im <- im[nrow(im):1,]

  # combine z stacks
  im <- abind::abind(list(im), along = 3)

  # put z-dimension first
  im <- aperm(im, perm = c(3, 1, 2))

  # Normalize images if requested
  if (normalize) {
    im <- normalize_images(im)
  }

  return(im)
}



## Pancreas data

load_pancreas_sc_rna <- function() {
  sc_obj <- DatasetPath("human_pancreas_rna", "tosti_2020", "p2.rds") %>% readRDS()
  sc_annot <- DatasetPath("human_pancreas_rna", "tosti_2020", "annotation_adj.csv") %>% read.csv()

  rownames(sc_annot) <- sc_annot$cell
  sc_annot <- sc_annot[rownames(sc_obj[["counts"]]),]

  ## making a seurat object to make dotplot
  counts <- Matrix::t(sc_obj[["misc"]][["rawCounts"]][rownames(sc_obj[["counts"]]),])
  sc_seurat <- Seurat::CreateSeuratObject(counts,meta.data = sc_annot)
  sc_seurat <- Seurat::NormalizeData(sc_seurat)
  sc_seurat@meta.data$cell_type <- sc_annot$annotation
  Idents(sc_seurat) <- sc_seurat@meta.data$cell_type

  return(list(sc_seurat=sc_seurat, sc_annot=sc_annot, sc_obj=sc_obj))
}

prepare_pancreas_sc_spatial <- function() {
  ### loading up scRNA
  sc_seurat <- load_pancreas_sc_rna()$sc_seurat

  ### load up spatial data in seurat object
  spatial_seurat <- DatasetPath("human_pancreas", "so.rds") %>% readRDS()
  spatial_annot <- DatasetPath("human_pancreas", "annotation.csv") %>% read.csv()
  colnames(spatial_annot)[2] <- 'celltype'
  spatial_seurat@meta.data$celltype <- spatial_annot$celltype
  spatial_annot <- cbind.data.frame(spatial_annot,spatial_seurat@meta.data[,c('x_centroid', 'y_centroid')])
  colnames(spatial_annot)[3:4] <- c('x','y')
  spatial_annot$z <- 1

  # load the molecule level data
  tx_dat <- DatasetPath("human_pancreas", "transcripts.parquet") %>% arrow::read_parquet() %>%
    rename(x=x_location, y=y_location, z=z_location, gene=feature_name, cell=cell_id)

  # subsetting data to same cells we have annotations for
  df <- tx_dat[tx_dat$cell %in% colnames(spatial_seurat),]
  match_ndx <- match(df$cell,colnames(spatial_seurat))
  df$celltype <- spatial_seurat@meta.data$celltype[match_ndx]

  df <- as.data.frame(df)

  # removing negative control probes
  negprb_ndx1 <- grep("Control", df$gene)
  negprb_ndx2 <- grep("Codeword", df$gene)
  negprb_ndx <- unique(c(negprb_ndx1,negprb_ndx2))
  df <- df[-negprb_ndx,]

  df$mol_id <- 1:nrow(df)
  rownames(df) <- df$mol_id

  ## need to remove all '/' in cell type names
  df$celltype[df$celltype=="Alpha/Beta/Delta/Gamma"] <- "Alpha-Beta-Delta-Gamma"
  spatial_annot$celltype[spatial_annot$celltype=="Alpha/Beta/Delta/Gamma"] <- "Alpha-Beta-Delta-Gamma"
  sc_levs <- levels(Idents(sc_seurat))
  ndx_change <- which(sc_levs=="Alpha/Beta/Delta/Gamma")
  sc_levs[ndx_change] <- "Alpha-Beta-Delta-Gamma"
  levels(Idents(sc_seurat)) <- sc_levs
  sc_seurat@meta.data$cell_type <- Idents(sc_seurat)

  return(list(df,spatial_annot,sc_seurat))
}

load_markers_pancreas <- function(filter_markers=TRUE) {
  # load markers
  markers <- DatasetPath("human_pancreas", "pancreas_markers_sub.csv") %>% read.csv()
  markers <- markers[,c('Gene','cell_type')]
  colnames(markers) <- c('gene','marker_of')

  g_rem <- markers$gene[duplicated(markers$gene)]
  markers <- markers[!(markers$gene %in% g_rem),]

  # filtering out cell types with very few marker genes
  if (filter_markers) {
    marker_counts <- table(markers$marker_of)
    ct_keep <- names(marker_counts)[marker_counts>2]
    markers <- markers[markers$marker_of %in% ct_keep,]
  }

  markers$marker_of[markers$marker_of=="Alpha/Beta/Delta/Gamma"] <- "Alpha-Beta-Delta-Gamma"
  return(markers)
}














