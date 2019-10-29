#' Plot of meta data with annotation of single cell data in
#' bivariate hexagon cells.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   or \code{\link[Seurat]{Seurat}} object.
#' @param col A string referring to the name of one column in the meta data of
#'   sce by which to make the outlines. Note that this should be a factor or
#'   a character.
#' @param col2 A string referring to the name of one column in the meta data of
#'   sce specifying what to color hexagons by.
#' @param action A string specifying how meta data as specified in col2 of
#'   observations in binned  hexagon cells are to be summarized. Possible
#'   actions are \code{prop}, \code{mode}, \code{mean} and \code{median}
#'   (see details).
#' @param no An integer specifying which level to plot of the column. Only in
#'   effect when \code{action=prop}.
#' @param colors A vector of strings specifying which colors to use for plotting
#'    the different levels in the selected column of the meta data.
#' @param title A string containing the title of the plot.
#' @param xlab A string containing the title of the x axis.
#' @param ylab A string containing the title of the y axis.
#' @param exapnd_hull A numeric value determining the expansion of the line
#'   marking different clusters.
#' @param ... Additional arguments passed on to
#'    \code{\link[geom_mark_hull]{ggforce}}.
#'
#' @details This function plots any gene expresssion in the hexagon cell
#'   representation calculated with \code{\link{make_hexbin}} as well as at the
#'   same time representing outlines of clusters. The chosen gene
#'   expression is summarized by one of four actions \code{prop_0},
#'   \code{mode}, \code{mean} and \code{median}:
#'
#'   \describe{
#'     \item{\code{prop}}{Returns the proportion of each level or unique
#'      character in the bin. The associated meta data column needs to be a
#'      factor or character.}
#'     \item{\code{mode}}{Returns the mode of the observations in the bin. The
#'      associated meta data column needs to be numeric.}
#'     \item{\code{mean}}{Returns the mean of the observations in the bin. The
#'      associated meta data column needs to be numeric.}
#'      \item{\code{median}}{Returns the median of the observations in the bin.
#'      The associated meta data column needs to be numeric.}
#'   }
#'
#' @return A \code{\link{ggplot2}{ggplot}} object.
#' @import Seurat
#' @import SingleCellExperiment
#' @import ggplot2
#' @importFrom dplyr as_tibble
#' @importFrom ggforce geom_mark_hull
#' @export
#'
#' @examples
#' #' # For Seurat object
#' library(Seurat)
#' data("pbmc_small")
#' pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
#' plot_hexbin_meta_plus(pbmc_small, col="RNA_snn_res.1",
#' col2="nFeature_RNA", action="mean")
setGeneric("plot_hexbin_meta_plus", function(sce,
                                             col,
                                             col2,
                                             action,
                                             no = 1,
                                             colors=NULL,
                                             title=NULL,
                                             xlab=NULL,
                                             ylab=NULL,
                                             expand_hull=3,
                                             ...) standardGeneric("plot_hexbin_meta_plus"))

#' @export
#' @describeIn plot_hexbin_meta_plus  Plot of gene expression and meta data
#'   of single cell data into hexagon cell for  SingleCellExperiment object.
setMethod("plot_hexbin_meta_plus", "SingleCellExperiment", function(sce,
                                                                    col,
                                                                    col2,
                                                                    action,
                                                                    no = 1,
                                                                    colors=NULL,
                                                                    title=NULL,
                                                                    xlab=NULL,
                                                                    ylab=NULL,
                                                                    expand_hull=3,
                                                                    ...) {

  out <- sce@metadata$hexbin[[2]]
  cID <- sce@metadata$hexbin[[1]]

  if(is.null(out)){
    stop("Compute hexbin representation before plotting.")
  }

  if (any(!col %in% colnames(colData(sce)))) {
    stop("Column cannot be found in colData(sce).")
  }

  if (any(!col2 %in% colnames(colData(sce)))) {
    stop("Column cannot be found in colData(sce).")
  }

  name_s <- paste0("sce$", col2)
  func_col2 <- paste0("x_col2 <- ", name_s)

  eval(parse(text = func_col2))

  name_s <- paste0("sce$", col)
  func <- paste0("x <- ", name_s)

  eval(parse(text = func))

  hh <- schex:::.make_hexbin_function(x, 'majority', cID)
  out <- as_tibble(out)

  col_hh <-paste0(col, "_", "majority")

  if(is.factor(x)){
    func1 <- paste0("out$", col_hh, " <- factor(hh, levels=",
                    "levels(x))")
  } else {
    func1 <- paste0("out$", col_hh, " <- hh")
  }

  eval(parse(text=func1))

  if(action == "prop"){
    col_hh_2 <- schex:::.make_hexbin_colnames(x_col2,col2)
    func1_col2 <- paste0("out$", col_hh_2,
                  " <- hh[,", seq(1,length(col_hh_2),1),"]")
    for(i in seq_len(length(func1_col2))){
      eval(parse(text=func1_col2[i]))
    }
  }
  } else {
    col_hh_2 <-paste0(col, "_", action)
    func1_col2 <- paste0("out$", col_hh_2, " <- hh")
    eval(parse(text=func1_col2))
  }

  .plot_hexbin_plus(out, colour_by = col_hh, fill_by_gene = col_hh_gene2,
                    colors=colors, expand_hull=expand_hull, title=title,
                    xlab=xlab, ylab=ylab, ...)

})

#' @export
#' @describeIn plot_hexbin_meta_plus Plot of gene expression and meta data of
#'   single cell data into hexagon cell for Seurat object.
setMethod("plot_hexbin_meta_plus", "Seurat", function(sce,
                                                      col,
                                                      col2,
                                                      action,
                                                      no = 1,
                                                      colors=NULL,
                                                      title=NULL,
                                                      xlab=NULL,
                                                      ylab=NULL,
                                                      expand_hull=3,
                                                      ...) {

  out <- sce@misc$hexbin[[2]]
  cID <- sce@misc$hexbin[[1]]

  if(is.null(out)){
    stop("Compute hexbin representation before plotting.")
  }

  if (any(!col %in% colnames(sce@meta.data))) {
    stop("Column cannot be found in slot(sce, 'meta.data').")
  }

  if (any(!col2 %in% colnames(sce@meta.data))) {
    stop("Column cannot be found in slot(sce, 'meta.data').")
  }

  name_s <- paste0("sce$", col2)
  func_col2 <- paste0("x_col2 <- ", name_s)

  eval(parse(text = func_col2))

  name_s <- paste0("sce$", col)
  func <- paste0("x <- ", name_s)

  eval(parse(text = func))

  hh <- schex:::.make_hexbin_function(x, 'majority', cID)
  out <- as_tibble(out)

  col_hh <-paste0(col, "_", "majority")

  if(is.factor(x)){
    func1 <- paste0("out$", col_hh, " <- factor(hh, levels=",
                    "levels(x))")
  } else {
    func1 <- paste0("out$", col_hh, " <- hh")
  }

  eval(parse(text=func1))

  if(action == "prop"){
    col_hh_2 <- .make_hexbin_colnames(x_col2,col2)
    func1_col2 <- paste0("out$", col_hh_2,
                         " <- hh[,", seq(1,length(col_hh_2),1),"]")
    for(i in seq_len(length(func1))){
      eval(parse(text=func1_col2[i]))
    }
  }
  } else {
    col_hh_2 <-paste0(col, "_", action)
    func1_col2 <- paste0("out$", col_hh_2, " <- hh")
    eval(parse(text=func1_col2))
  }

  .plot_hexbin_plus(out, colour_by = col_hh, fill_by_gene = col_hh_2,
                  colors=colors, expand_hull=expand_hull, title=title,
                  xlab=xlab, ylab=ylab, ...)

})
