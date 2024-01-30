#' Plot of gene expression and meta data of single cell data in
#' bivariate hexagon cells.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param col A string referring to the name of one column in the meta data of
#'    sce by which to colour the hexagons.
#' @param mod A string referring to the name of the modality used for plotting.
#'     For RNA modality use "RNA". For other modalities use name of alternative
#'     object for the \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'     object.
#' @param type A string referring to the type of assay in the
#'   \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @param feature A string referring to the name of one feature.
#' @param action A string specifying how gene expression of observations in
#'   binned  hexagon cells are to be summarized. Possible actions are
#'  \code{prop_0}, \code{mode}, \code{mean} and
#'   \code{median} (see details).
#' @param colors A vector of strings specifying which colors to use for plotting
#'    the different levels in the selected column of the meta data.
#' @param title A string containing the title of the plot.
#' @param xlab A string containing the title of the x axis.
#' @param ylab A string containing the title of the y axis.
#' @param expand_hull A numeric value determining the expansion of the line
#'   marking different clusters.
#' @param ... Additional arguments passed on to
#'    \code{\link{ggforce}{geom_mark_hull}}.
#'
#' @details This function plots any gene expresssion in the hexagon cell
#'   representation calculated with \code{\link{make_hexbin}} as well as at the
#'   same time representing outlines of clusters. The chosen gene
#'   expression is summarized by one of four actions \code{prop_0},
#'   \code{mode}, \code{mean} and \code{median}:
#'
#'   \describe{
#'     \item{\code{prop_0}}{Returns the proportion of observations in the bin
#'      greater than 0. The associated meta data column needs to be numeric.}
#'     \item{\code{mode}}{Returns the mode of the observations in the bin. The
#'      associated meta data column needs to be numeric.}
#'     \item{\code{mean}}{Returns the mean of the observations in the bin. The
#'      associated meta data column needs to be numeric.}
#'      \item{\code{median}}{Returns the median of the observations in the bin.
#'      The associated meta data column needs to be numeric.}
#'   }
#'
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
#' tenx_pbmc3k <- tenx_pbmc3k[!rm_ind, ]
#' tenx_pbmc3k <- logNormCounts(tenx_pbmc3k)
#' tenx_pbmc3k <- runPCA(tenx_pbmc3k)
#' tenx_pbmc3k <- make_hexbin(tenx_pbmc3k, 10, dimension_reduction = "PCA")
#' tenx_pbmc3k$random <- factor(sample(1:3, ncol(tenx_pbmc3k), replace = TRUE))
#' plot_hexbin_feature_plus(tenx_pbmc3k,
#'     col = "random", type = "counts",
#'     feature = "ENSG00000135250", action = "mean"
#' )
plot_hexbin_feature_plus <- function(
    sce,
    col,
    mod = "RNA",
    type,
    feature,
    action,
    colors = NULL,
    title = NULL,
    xlab = NULL,
    ylab = NULL,
    expand_hull = 3,
    ...) {
    out <- .extract_hexbin(sce)
    cID <- .extract_cID(sce)

    if (is.null(out)) {
        stop("Compute hexbin representation before plotting.")
    }

    x_gene <- .prepare_data_feature(sce, mod, type, feature)

    hh_gene <- .make_hexbin_function(x_gene, action, cID)

    x <- .prepare_data_meta(sce, col)

    hh <- .make_hexbin_function(x, "majority", cID)
    out <- as_tibble(out)

    if (is.factor(x)) {
        out$meta <- factor(hh, levels = levels(x))
    } else {
        out$meta <- hh
    }

    out$feature <- hh_gene

    if (is.null(title)) {
        title <- paste0(col, "_majority_", feature, "_", action)
    }

    .plot_hexbin_plus(out,
        colour_by = "meta", fill_by_gene = "feature",
        colors = colors, expand_hull = expand_hull, title = title,
        xlab = xlab, ylab = ylab, ...
    )
}


.plot_hexbin_plus <- function(drhex, colour_by = "meta", fill_by_gene,
                              colors = NULL, expand_hull = 3, legend = legend,
                              title = NULL, xlab = NULL, ylab = NULL, ...) {
    if (any(!c("x", "y", colour_by) %in% colnames(drhex))) {
        stop("The dataframe must contain columns named 'x', 'y' and col.")
    }

    if (is.null(xlab)) {
        xlab <- "x"
    }

    if (is.null(ylab)) {
        ylab <- "y"
    }

    if (is.null(colors)) {
        ggplot(drhex, aes_string(x = "x", y = "y", fill = fill_by_gene)) +
            geom_hex(stat = "identity") +
            geom_mark_hull(aes_string(label = colour_by, col = colour_by),
                show.legend = FALSE, expand = unit(expand_hull, "mm"),
                fill = NA, size = 2, ...
            ) +
            theme_classic() +
            scale_fill_viridis_c() +
            ggtitle(title) +
            labs(x = xlab, y = ylab) +
            theme(legend.title = element_blank())
    } else {
        ggplot(drhex, aes_string(x = "x", y = "y", fill = fill_by_gene)) +
            geom_hex(stat = "identity") +
            geom_mark_hull(aes_string(label = colour_by, col = colour_by),
                show.legend = FALSE, expand = unit(expand_hull, "mm"),
                fill = NA, size = 2, ...
            ) +
            theme_classic() +
            scale_fill_viridis_c() +
            ggtitle(title) +
            labs(x = xlab, y = ylab) +
            theme(legend.title = element_blank()) +
            scale_color_manual(values = colors)
    }
}
