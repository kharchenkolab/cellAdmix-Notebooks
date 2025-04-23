process_seurat <- function(so, n.od.genes=2000, build.knn=TRUE) {
    so %<>%
      Seurat::NormalizeData()

    if (!is.null(n.od.genes)) {
      so %<>% Seurat::FindVariableFeatures(selection.method="vst", nfeatures=n.od.genes)
    }

    so %<>%
      Seurat::ScaleData() %>%
      Seurat::RunPCA(., features=VariableFeatures(.))

    if (!build.knn) {
      return(so)
    }

    so %<>%
      Seurat::FindNeighbors(dims=1:30) %>%
      Seurat::RunUMAP(dims=1:30, n.epochs=500) %>%
      Seurat::FindClusters(resolution=2, algorithm=2)

    return(so)
}

transfer_annotation_knn <- function(pca, cell.types, ref.prefix, k=7, n.cores=50) {
  pca.ref <- pca[startsWith(rownames(pca), ref.prefix),]
  pca.target <- pca[!(rownames(pca) %in% rownames(pca.ref)),]

  nn.mat <- N2R::crossKnn(pca.target, pca.ref, k=k, nThreads=n.cores, verbose=FALSE, indexType="angular") %>%
    as('TsparseMatrix')
  ref.ids.per.cell <- nn.mat %>%
      {split(.@i + 1, rownames(pca.target)[.@j + 1])}

  cell.types.ref <- cell.types[rownames(pca.ref)]

  cell.types <- sapply(ref.ids.per.cell, \(cids) {
      cell.types.ref[cids] %>% table() %>% which.max() %>% names()
  })

  return(cell.types)
}