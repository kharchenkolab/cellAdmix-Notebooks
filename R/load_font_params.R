#' @export
load_arial_font <- function(font.dir=system.file('data', package='cellAdmixNotebooks')) {
  if (!dir.exists(font.dir)) {
    warning(
      "Arial fonts don't exist.",
      "Please, place the fonts in the data/fonts directory or provide a valid path.",
      "You can use https://befonts.com/arial-font.html to download them (`arial.ttf`)."
    )
    return(NULL)
  }
  extrafont::font_import(paths=font.dir, prompt=FALSE)
  extrafont::loadfonts(quiet=TRUE)
}

font_choice <- "Arial"

#' @export
p_theme <- ggplot2::theme(
  text = ggplot2::element_text(family = font_choice),
  axis.text=ggplot2::element_text(size=6), axis.title=ggplot2::element_text(size=8),
  legend.title=ggplot2::element_text(size=8), legend.text=ggplot2::element_text(size=6), legend.key.width=grid::unit(8, "pt"),
  plot.title = ggplot2::element_text(hjust = 0.5, size = 8)
)

#' @export
geom_text_size <- round((5/15) * 6)

#' @export
gg_line_thickness <- .5/2.15

#' @export
gg_line_alpha <- .8

#' @export
gg_point_alpha <- .8

#' @export
gg_point_size <- .75

#' @export
legend_mod <- ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=1.5)))


