#' Plot of fold change of selected gene in single cell data using
#'    bivariate hexagon cells.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'    object.
#' @param col A string referring to the name of one column in the meta data of
#'    sce by which to compare. Note this factor can only contain two levels.
#' @param mod A string referring to the name of one column in the meta data of
#'    sce by which to compare. Note this factor can only contain two levels.
#' @param type A string referring to the name of one column in the meta data of
#'    sce by which to compare. Note this factor can only contain two levels.
#' @param feature A string referring to the name of one feature.
#' @param title A string containing the title of the plot.
#' @param xlab A string containing the title of the x axis.
#' @param ylab A string containing the title of the y axis.
#' @param colors A vector of strings specifying which colors to use for plotting
#'    the different levels in the selected column of the meta data.
#'
#' @details This function plots fold change within each
#'   hexagon, which are calculated with \code{\link{make_hexbin}}.
#'   Note that the fold change is only accurate if the condition
#'   investigated is within the same individual. For conditions across
#'   different individuals different methods that account for
#'   individual-specific effects are required.
#'
#' @return A \code{\link{ggplot2}{ggplot}} object.
#' @import SingleCellExperiment
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
#' tenx_pbmc3k <- make_hexbin(tenx_pbmc3k, 20, dimension_reduction = "PCA")
#' tenx_pbmc3k$random <- factor(sample(1:2, ncol(tenx_pbmc3k), replace = TRUE))
#' plot_hexbin_fc(tenx_pbmc3k, col = "random", feature = "ENSG00000187608", type = "counts")
plot_hexbin_fc <- function(
    sce,
    col,
    mod = "RNA",
    type,
    feature,
    title = NULL,
    xlab = NULL,
    ylab = NULL,
    colors) {
    out <- .extract_hexbin(sce)
    cID <- .extract_cID(sce)

    if (is.null(out)) {
        stop("Compute hexbin representation before plotting.")
    }

    x <- .prepare_data_meta(sce, col)
    x_gene <- .prepare_data_feature(sce, mod, type, feature)

    .plot_hexbin_fc_helper(
        x, x_gene, feature, out, cID, col, title,
        xlab, ylab, colors
    )
}

.plot_hexbin_fc_helper <- function(
    x, x_gene, feature, out, cID, col,
    title, xlab, ylab, colors) {
    hh <- .make_hexbin_fc_function(x, x_gene, cID)
    ind_null <- unlist(lapply(hh, is.null))
    hh[ind_null] <- NA
    hh <- unlist(hh)
    out <- as_tibble(out)

    if (is.null(title)) {
        title <- paste0(feature, "_", "fc")
    }

    out$feature <- hh

    .plot_hexbin(out,
        colour_by = "feature", action = "fc",
        title = title, xlab = xlab, ylab = ylab
    )
}
