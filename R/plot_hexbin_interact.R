#' Plot of interaction of expression of single cells in bivariate hexagon cells.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @param mod A vector of strings referring to the names of the modularities.
#'    For \code{\link[SingleCellExperiment]{SingleCellExperiment}} use "RNA" to
#'    access the RNA expression data stored as the main experiment type.
#' @param type A vector of strings referring to the types of assays in the
#'    \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @param feature A vector of strings referring to the names of one features in
#'    the same order as the vector of modularities.
#' @param interact A string specifying how interaction between features is
#'    calculated. Possible interaction measures are
#'    \code{corr_spearman} and \code{mi} (see details).
#' @param title A string containing the title of the plot.
#' @param xlab A string containing the title of the x axis.
#' @param ylab A string containing the title of the y axis.
#'
#' @details This function plots the interaction between any features in the
#'    hexagon cell representation calculated with \code{\link{make_hexbin}}. The
#'    interaction between the chosen features is calculated by one of two
#'    measurers \code{corr_spearman}, \code{ratio} and \code{mi}:
#'
#'    \describe{
#'       \item{\code{mi}}{Returns the mutual information coefficient.}
#'       \item{\code{corr_spearman}}{Returns the Spearman correlation.}
#'       \item{\code{fc}}{Return the log fold change between the features.}
#'    }
#'
#'    Note that \code{fc} should be applied to log normalized values.
#'
#' @return A \code{\link{ggplot2}{ggplot}} object.
#' @import SingleCellExperiment
#' @importFrom methods slotNames
#' @import ggplot2
#' @importFrom dplyr as_tibble
#' @export
#'
#' @examples
#' # For SingleCellExperiment
#' library(TENxPBMCData)
#' library(scater)
#' tenx_pbmc3k <- TENxPBMCData(dataset = "pbmc3k")
#' rm_ind <- calculateAverage(tenx_pbmc3k) < 0.1
#' tenx_pbmc3k <- tenx_pbmc3k[!rm_ind, ]
#' colData(tenx_pbmc3k) <- cbind(colData(tenx_pbmc3k), perCellQCMetrics(tenx_pbmc3k))
#' tenx_pbmc3k <- logNormCounts(tenx_pbmc3k)
#' tenx_pbmc3k <- runPCA(tenx_pbmc3k)
#' tenx_pbmc3k <- make_hexbin(tenx_pbmc3k, 10, dimension_reduction = "PCA")
#' plot_hexbin_interact(tenx_pbmc3k,
#'     type = c("counts", "counts"), mod = c("RNA", "RNA"),
#'     feature = c("ENSG00000146109", "ENSG00000102265"), interact = "fc"
#' )
plot_hexbin_interact <- function(
    sce,
    mod,
    type,
    feature,
    interact,
    title = NULL,
    xlab = NULL,
    ylab = NULL) {
    if (length(mod) != length(feature) | length(feature) != length(type)) {
        stop("Specify the same number of modularities, types and features.")
    }

    out <- .extract_hexbin(sce)
    cID <- .extract_cID(sce)

    if (is.null(out)) {
        stop("Compute hexbin representation before plotting.")
    }

    first_x <- .prepare_data_feature(sce, mod[1], type[1], feature[1])

    second_x <- .prepare_data_feature(sce, mod[2], type[2], feature[2])

    .plot_hexbin_interact_helper(
        first_x, second_x, out, cID, interact,
        feature, title, xlab, ylab
    )
}

.plot_hexbin_interact_helper <- function(
    first_x, second_x, out, cID, interact,
    feature, title, xlab, ylab) {
    hh <- .interact_hexbin_function(first_x, second_x, interact, cID)
    out <- as_tibble(out)

    if (is.null(title)) {
        title <- paste0(interact, "_", feature[1], "_", feature[2])
    }

    out$interact <- hh

    .plot_hexbin(out,
        colour_by = "interact", action = "interact",
        title = title, xlab = xlab, ylab = ylab
    )
}
