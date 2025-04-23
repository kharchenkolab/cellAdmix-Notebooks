#' Helper function for creating color palettes
#'
#' @param name character A palette name; please refer to RColorBrewer::brewer.pal()
#' @param n integer (default=NULL) Number of different colors in the palette, minimum 3, maximum depending on palette. Please refer to RColorBrewer::brewer.pal() for more details
#' @param rev boolean (default=TRUE) Whether to reverse the palette order
#' @return palette function
#'
#' @export
brewerPalette <- function(name, n=NULL, rev=TRUE) {
  sccore::checkPackageInstalled("RColorBrewer", cran=TRUE)
  if (is.null(n)) n <- RColorBrewer::brewer.pal.info[name,]$maxcolors
  pal <- RColorBrewer::brewer.pal(n, name)
  if (rev) pal <- rev(pal)
  return(grDevices::colorRampPalette(pal))
}

adjustColorsByBrightness <- function(rgb.vector, brightness.vector) {
  brightness.vector <- as.vector(brightness.vector)
  color.vec <- t(rgb.vector) * brightness.vector
  bg.vec <- matrix(rep(255, 3 * length(brightness.vector)), ncol=3) * (1 - brightness.vector)
  color.vec <- t(color.vec + bg.vec) %>%
    pmin(255) %>% pmax(0) # to avoid floating problems

  return(color.vec)
}

# stain values below clip.min are set to 0, and values above clip.max are set to 1.
#' @export
combinedStainImage <- function(stains, clip.min=c(0.25, 0.25), clip.max=c(1.0, 1.0), alpha=0.75, palettes=c("Greens", "Purples"), verbose=TRUE) {
  if (length(stains) != 2) {
    stop("Only two stains are supported")
  }

  if (verbose) {
    message("Normalizing stains")
  }


  stains.norm <- mapply(\(st, cl.min, cl.max) {
    st %>% {. - quantile(., cl.min)} %>%
      {. / quantile(., cl.max)} %>%
      pmax(0) %>% pmin(1) %>% as.matrix()
  }, stains, clip.min, clip.max, SIMPLIFY=FALSE)

  if (is.character(palettes)) {
    palettes <- lapply(palettes, brewerPalette, rev=FALSE)
  }

  # if (is.character(palettes)) {
  #   palettes <- lapply(palettes,function(x){
  #     if (x=='Purples') {
  #       return(colorRampPalette(colors = c("black", "#020af2")))
  #     } else if (x=='Greens') {
  #       return(colorRampPalette(colors = c("black", "#05f400")))
  #     }
  #   })
  # }

  if (verbose) {
    message("Creating combined stain image")
  }

  n.colors <- 100
  stain.cols <- mapply(\(stain, pal) {
    stain %>% {round(. * (n.colors - 1)) + 1} %>% {pal(n.colors)[.]} %>% col2rgb()
  }, stains.norm, palettes, SIMPLIFY=FALSE)

  if (verbose) {
    message("Adjusting colors by brightness")
  }

  stain.cols <- mapply(adjustColorsByBrightness, stain.cols, stains.norm, SIMPLIFY=FALSE)

  weights <- lapply(stains.norm, \(st) pmax(as.vector(st), 0.01))

  if (verbose) {
    message("Combining stains")
  }

  stain.cols %<>% {t((t(.[[1]]) * weights[[1]] + t(.[[2]]) * weights[[2]]) / (weights[[1]] + weights[[2]])) / 255} %>%
    pmin(1) %>% pmax(0) %>%
    {rgb(.[1,], .[2,], .[3,], alpha)} %>%
    matrix(ncol=ncol(stains[[1]])) %>% t()

  if (verbose) {
    message("Done!")
  }

  return(stain.cols)
}


#' @export
extractStainRegionAnnot <- function(stain, region.df, expand=0, clip=1.0, inv.y.global=FALSE, inv.y.local=FALSE) {
  xr <- range(region.df$x) %>% pmax(1)
  yr <- range(region.df$y) %>% pmax(1)

  if (expand > 0) {
    dx <- diff(xr) * expand
    xr <- xr + c(-dx, dx)
    dy <- diff(yr) * expand
    yr <- yr + c(-dy, dy)
  }

  if (inv.y.global) {
    stain %<>% .[(nrow(.)-yr[2]):(nrow(.)-yr[1]), xr[1]:xr[2]]
  } else if (inv.y.local) {
    stain %<>% .[seq(yr[2], yr[1]), xr[1]:xr[2]]
  } else {
    stain %<>% .[yr[1]:yr[2], xr[1]:xr[2]]
  }


  if (clip < 1.0) {
    stain %<>%
      {. - quantile(., 1 - clip)} %>% pmax(0) %>%
      {. / quantile(., clip)} %>% pmin(1)
  }

  return(annotation_raster(stain, xmin=xr[1], xmax=xr[2], ymin=yr[1], ymax=yr[2]))
}


#' @export
extractAdjacentMolecules <- function(df.spatial, cell.id, offset=0.25, add.cell.flag=TRUE) {
  bbox = df.spatial %>% filter(cell == cell.id) %$% rbind(range(x), range(y))

  offset = t(bbox) %>% diff() %>% {. * offset}
  bbox = t(offset) %>% cbind(-., .) %>% {. + bbox}

  res.mols = df.spatial %>%
    filter(x >= bbox[1, 1], x <= bbox[1, 2], y >= bbox[2, 1], y <= bbox[2, 2])

  if (add.cell.flag) {
    res.mols %<>% mutate(cell_flag=ifelse(cell == cell.id, cell.id, "other"))
  }

  return(res.mols)
}

#' @export
plotMolecules <- function(
    df.spatial, color=NULL, fill=NULL, annotation=NULL, shape=NULL, stain=NULL,
    alpha=1.0, scale.z=TRUE, size.range=c(1, 3), size=NULL, show.ticks=FALSE,
    tooltip="gene", use.raw.color=FALSE, unclean_cell_for_bounds=NULL, expand=c(0, 0), color.title=NULL
) {
  if (!is.null(annotation)) {
    message("The 'annotation' argument is deprecated. Use 'color' instead.")
    color <- annotation
  }

  # if (!is.null(color) && !is.null(annotation))
  #   stop("Only one of 'color' and 'annotation' can be provided")

  if (!("z" %in% colnames(df.spatial))) {
    scale.z <- FALSE
  }

  p.aes <- aes(group=.data[[tooltip]])
  if (scale.z) {
    p.aes %<>% modifyList(aes(size=z))
  }

  if (length(shape) == 1) {
    shape <- df.spatial[[shape]]
  }

  if (!is.null(shape)) {
    p.aes %<>% modifyList(aes(shape=shape))
  }

  if (length(color) == 1) {
    color <- df.spatial[[color]]
  }

  if (!is.null(fill)) {
    p.aes %<>% modifyList(aes(fill=.data[[fill]]))
  }

  if (!is.null(size) & !is.numeric(size)) {
    size <- df.spatial[[size]]
    p.aes %<>% modifyList(aes(size=size))
  }

  if (!use.raw.color) {
    p.aes %<>% modifyList(aes(color=color))
  }

  g.point <- geom_point(p.aes, alpha=alpha)
  if (is.null(size)) {
    g.point$aes_params$size <- mean(size.range)
  } else if (is.numeric(size)) {
    g.point$aes_params$size <- size
  }

  if (use.raw.color) {
    g.point$aes_params$colour <- color
  }

  gg <- ggplot(df.spatial, aes(x=x, y=y))

  if (!is.null(stain)) {
    # gg <- gg + annotation_raster(stain, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
    gg <- gg + stain
  }

  gg <- gg + g.point
  if (scale.z) {
    gg$layers[[length(gg$layers)]]$aes_params$size <- NULL # Remove fixed size set above
    gg <- gg + scale_size_continuous(range=size.range)
  }

  if (!show.ticks) {
    gg <- gg + theme(axis.ticks=element_blank(), axis.text=element_blank())
  }

  if (!is.null(unclean_cell_for_bounds)) {
    gg <- gg +
      scale_x_continuous(limits=range(unclean_cell_for_bounds$x), expand=expand) +
      scale_y_continuous(limits=range(unclean_cell_for_bounds$y), expand=expand)
  } else {
    gg <- gg +
      scale_x_continuous(limits=range(df.spatial$x), expand=expand) +
      scale_y_continuous(limits=range(df.spatial$y), expand=expand)
  }

  if (!is.null(color.title)) {
    gg <- gg + guides(color=guide_legend(title=color.title))
  }

  return(gg)
}


# Function to compute alpha shape polygon vertices and extend them outward by a buffer
#' @export
alpha_polygon <- function(x, y, alpha_val = 1.5, buffer = 0.2) {
  # Compute the alpha shape
  ashape_obj <- ashape(x, y, alpha = alpha_val)

  # Convert edges to a data frame
  edges <- as.data.frame(ashape_obj$edges)

  # Extract the unique endpoints of the alpha shape edges
  vertices <- unique(rbind(
    data.frame(x = edges$x1, y = edges$y1),
    data.frame(x = edges$x2, y = edges$y2)
  ))

  # Compute the centroid using the original points
  center <- c(mean(x), mean(y))

  # Order the vertices by angle around the centroid to form a continuous polygon
  vertices$order <- atan2(vertices$y - center[2], vertices$x - center[1])
  vertices <- vertices[order(vertices$order), ]

  # For each vertex, compute the distance to the centroid, then push the vertex outward by the buffer distance.
  vertices <- vertices %>%
    mutate(
      d = sqrt((x - center[1])^2 + (y - center[2])^2),
      x_new = x + buffer * (x - center[1]) / d,
      y_new = y + buffer * (y - center[2]) / d
    )

  # Return the new, buffered vertices
  return(data.frame(x = vertices$x_new, y = vertices$y_new))
}

#' @export
plotCellAdmixture <- function(
    df.spatial, stains, polygons=NULL, cell.id, z.id,
    filter.z=TRUE, size=2, offset=0.25, expand=0,
    size.scales=NULL, color.scales=NULL, shape.scales=NULL,
    bg.name="other", scale.z=FALSE, fill="source", shape="source",
    shape.name="Marker type", unclean_cell_for_bounds=NULL, inv.y.global=FALSE, inv.y.local=FALSE, ...
) {

  poly <- polygons[[cell.id]] %>% as.data.frame()
  p.df <- extractAdjacentMolecules(df.spatial, cell.id=cell.id, offset=offset)

  if (is.null(stains)) {
    stains <- NULL
  } else {
    # Subsetting is required for export, as otherwise it stores the whole image
    if (is.list(stains)) {
      stains <- stains[[paste(z.id)]]
    }

    if (!is.null(unclean_cell_for_bounds)) {
      p.tmp <- extractAdjacentMolecules(unclean_cell_for_bounds, cell.id=cell.id, offset=offset)
      stains %<>% extractStainRegionAnnot(p.tmp, expand=expand, inv.y.global=inv.y.global, inv.y.local=inv.y.local)
    } else {
      stains %<>% extractStainRegionAnnot(p.df, expand=expand, inv.y.global=inv.y.global, inv.y.local=inv.y.local)
    }
  }

  if (filter.z) {
    p.df %<>% filter(z==z.id)
  }

  gg <- plotMolecules(
    p.df, fill=fill, shape=shape, size=(if(scale.z) NULL else "source"),
    alpha=ifelse(p.df$source == bg.name, 0.5, 1.0), scale.z=scale.z,
    unclean_cell_for_bounds=unclean_cell_for_bounds, stain=stains,
    expand=c(expand, 0), ...
  )

  if (!is.null(poly) && nrow(poly) > 0) {
    gg <- gg +
      # geom_polygon(aes(x=x, y=y, group=cell), data=poly, fill = NA, color="black")
    geom_polygon(data = poly, aes(x = x, y = y),color='black',
                 alpha = 0.001, show.legend = FALSE)
  }

  if (is.null(size.scales)) {
    size.scales <- rep(size, length(unique(p.df$source))) %>%
      setNames(unique(p.df$source))

    size.scales[[bg.name]] <- 0.25 * size
  }

  if (!scale.z) {
    gg <- gg + scale_size_manual(values=size.scales, name="Marker type")
  }

  if (!is.null(color.scales)) {
    gg <- gg + scale_fill_manual(values=color.scales, name="Marker type")
  }

  if (!is.null(shape.scales)) {
    gg <- gg + scale_shape_manual(values=shape.scales, name=shape.name)
  }

  gg
}