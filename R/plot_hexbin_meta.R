#' Plot of meta data of single cell data in bivariate hexagon cells.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'    or \code{\link[Seurat]{Seurat-class}} object.
#' @param col A string referring to the name of one column in the meta data of
#'    sce by which to colour the hexagons.
#' @param action A string specifying how meta data of observations in
#'    binned  hexagon cells are to be summarized. Possible actions are
#'    \code{majority}, \code{prop}, \code{prop_0}, \code{mode}, \code{mean} and
#'    \code{median} (see details).
#' @param no An integer specifying which level to plot of the column. Only in
#'    effect when \code{action=prop}.
#' @param colors A vector of strings specifying which colors to use for plotting
#'    the different levels in the selected column of the meta data. Only in
#'    effect when the selected \code{action="majority"}.
#' @param title A string containing the title of the plot.
#' @param xlab A string containing the title of the x axis.
#' @param ylab A string containing the title of the y axis.
#' @param na.rm Logical indicating whether NA values should be removed.
#'
#' @details This function plots any column of the meta data in the hexagon cell
#'    representation calculated with \code{\link{make_hexbin}}. The chosen meta
#'    data column is summarized by one of six actions \code{majority},
#'    \code{prop}, \code{prop_0}, \code{mode}, \code{mean} and \code{median}:
#'
#'    \describe{
#'       \item{\code{majority}}{Returns the value of the majority of
#'       observations in the bin. The associated meta data column needs to be
#'       a factor or character.}
#'       \item{\code{prop}}{Returns the proportion of each level or unique
#'       character in the bin. The associated meta data column needs to be a
#'       factor or character.}
#'       \item{\code{prop_0}}{Returns the proportion of observations in the bin
#'       greater than 0. The associated meta data column needs to be numeric.}
#'       \item{\code{mode}}{Returns the mode of the observations in the bin. The
#'       associated meta data column needs to be numeric.}
#'       \item{\code{mean}}{Returns the mean of the observations in the bin. The
#'       associated meta data column needs to be numeric.}
#'       \item{\code{median}}{Returns the median of the observations in the bin.
#'       The associated meta data column needs to be numeric.}
#'   }
#'
#' @return A \code{\link{ggplot2}{ggplot}} object.
#' @import Seurat
#' @import SingleCellExperiment
#' @import ggplot2
#' @importFrom dplyr as_tibble
#' @export
#'
#' @examples
#' #' # For Seurat object
#' library(Seurat)
#' data("pbmc_small")
#' pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
#' plot_hexbin_meta(pbmc_small, col="RNA_snn_res.1", action="prop", no=1)
#' # For SingleCellExperiment object
#' \dontrun{
#' library(TENxPBMCData)
#' library(scater)
#' tenx_pbmc3k <- TENxPBMCData(dataset = "pbmc3k")
#' rm_ind <- calculateAverage(tenx_pbmc3k) < 0.1
#' tenx_pbmc3k <- tenx_pbmc3k[-rm_ind, ]
#' colData(tenx_pbmc3k) <- cbind(
#'   colData(tenx_pbmc3k),
#'   perCellQCMetrics(tenx_pbmc3k)
#' )
#' tenx_pbmc3k <- logNormCounts(tenx_pbmc3k)
#' tenx_pbmc3k <- runPCA(tenx_pbmc3k)
#' tenx_pbmc3k <- make_hexbin(tenx_pbmc3k, 20, dimension_reduction = "PCA")
#' plot_hexbin_meta(tenx_pbmc3k, col = "total", action = "median")
#' }
plot_hexbin_meta <- function(sce,
    col,
    action,
    no = 1,
    colors=NULL,
    title=NULL,
    xlab=NULL,
    ylab=NULL,
    na.rm=FALSE){
  
  out <- .extract_hexbin(sce)
  cID <- .extract_cID(sce)
  
  x <- .prepare_data_meta(sce, col)
  
  .plot_hexbin_meta_helper(x, out, cID, col, action, no, title, xlab, ylab,
                           colors, na.rm)
}


.plot_hexbin_meta_helper <- function(x, out, cID, col, action, no, title,
                                     xlab, ylab, colors, na.rm) {
  if (is.null(out)) {
    stop("Compute hexbin representation before plotting.")
  }

  hh <- .make_hexbin_function(x, action, cID, na.rm)
  out <- as_tibble(out)

  if (action == "prop" | action == "majority") {
    if (action == "prop") {
      col_hh <- .make_hexbin_colnames(x, col)
      func1 <- paste0(
        "out$", col_hh, " <- hh[,",
        seq(1, length(col_hh), 1), "]"
      )
      for (i in seq_len(length(func1))) {
        eval(parse(text = func1[i]))
      }
    }
    if (action == "majority") {
      col_hh <- paste0(col, "_", action)
      if (is.factor(x)) {
        func1 <- paste0(
          "out$", col_hh, " <- factor(hh, levels=",
          "levels(x))"
        )
      } else {
        func1 <- paste0("out$", col_hh, " <- hh")
      }
      eval(parse(text = func1))
    }
  } else {
    col_hh <- paste0(col, "_", action)
    func1 <- paste0("out$", col_hh, " <- hh")
    eval(parse(text = func1))
  }

  if (action != "prop") {
    if (action == "majority") {
      .plot_hexbin(out,
        colour_by = col_hh, colors = colors,
        title = title, xlab = xlab, ylab = ylab
      )
    } else {
      .plot_hexbin(out,
        colour_by = col_hh, colors = NULL,
        title = title, xlab = xlab, ylab = ylab
      )
    }
  } else {
    .plot_hexbin(out,
      colour_by = col_hh[no], colors = NULL,
      title = title, xlab = xlab, ylab = ylab
    )
  }
}
