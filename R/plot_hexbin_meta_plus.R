#' Plot of meta data with annotation of single cell data in
#'    bivariate hexagon cells.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @param col1 A string referring to the name of one column in the meta data of
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
#' @param expand_hull A numeric value determining the expansion of the line
#'   marking different clusters.
#' @param na.rm Logical indicating whether NA values should be removed.
#' @param ... Additional arguments passed on to
#'    \code{\link{ggforce}{geom_mark_hull}}.
#'
#' @details This function plots any meta data in the hexagon cell
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
#' @import SingleCellExperiment
#' @import ggplot2
#' @importFrom dplyr as_tibble
#' @importFrom ggforce geom_mark_hull
#' @import concaveman
#' @export
#'
#' @examples
#' # For SingleCellExperiment object
#' library(TENxPBMCData)
#' library(scater)
#' tenx_pbmc3k <- TENxPBMCData(dataset = "pbmc3k")
#' rm_ind <- calculateAverage(tenx_pbmc3k) < 0.1
#' tenx_pbmc3k <- tenx_pbmc3k[-rm_ind, ]
#' colData(tenx_pbmc3k) <- cbind(colData(tenx_pbmc3k), perCellQCMetrics(tenx_pbmc3k))
#' tenx_pbmc3k <- logNormCounts(tenx_pbmc3k)
#' tenx_pbmc3k <- runPCA(tenx_pbmc3k)
#' tenx_pbmc3k <- make_hexbin(tenx_pbmc3k, 20, dimension_reduction = "PCA")
#' tenx_pbmc3k$random <- factor(sample(1:3, ncol(tenx_pbmc3k), replace = TRUE))
#' tenx_pbmc3k$random <- as.factor(tenx_pbmc3k$random)
#' plot_hexbin_meta_plus(tenx_pbmc3k, col1 = "random", col2 = "total", action = "median")
plot_hexbin_meta_plus <- function(
    sce,
    col1,
    col2,
    action,
    no = 1,
    colors = NULL,
    title = NULL,
    xlab = NULL,
    ylab = NULL,
    expand_hull = 3,
    na.rm = FALSE,
    ...) {
    out <- .extract_hexbin(sce)
    cID <- .extract_cID(sce)

    if (is.null(out)) {
        stop("Compute hexbin representation before plotting.")
    }

    if (is.null(title)) {
        title <- paste0(col1, "_majority", "_", col2, "_", action)
    }

    x_col2 <- .prepare_data_meta(sce, col2)
    x <- .prepare_data_meta(sce, col1)

    hh <- .make_hexbin_function(x, "majority", cID, na.rm)
    hh2 <- .make_hexbin_function(x_col2, action, cID, na.rm)
    out <- as_tibble(out)

    if (is.factor(x)) {
        out$meta <- factor(hh, levels = levels(x))
    } else {
        out$meta <- hh
    }

    out$meta2 <- hh2

    .plot_hexbin_plus(out,
        colour_by = "meta", fill_by_gene = "meta2",
        colors = colors, expand_hull = expand_hull, title = title,
        xlab = xlab, ylab = ylab, ...
    )
}
