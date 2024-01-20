#' Plot of feature expression and uncertainty of single cells in bivariate
#'    hexagon cells.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @param mod A string referring to the name of the modality used for plotting.
#'     For RNA modality use "RNA". For other modalities use name of alternative
#'     object for the \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'     object.
#' @param type A string referring to the type of assay in the
#'    \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @param feature A string referring to the name of one feature.
#' @param fan Logical indicating whether to plot uncertainty surpressing palette.
#' @param title A string containing the title of the plot.
#' @param xlab A string containing the title of the x axis.
#' @param ylab A string containing the title of the y axis.
#'
#' @details This function plots the mean expression and a measure of uncertainty
#'    of any feature in the hexagon cell representation calculated with
#'    \code{\link{make_hexbin}} using a bivariate scale. When \code{fan=FALSE},
#'    the standard deviation and the mean expression are plotted. When
#'    \code{fan=TRUE}, the mean expression and coefficient of variation are
#'    calculated. The coefficient of variation is converted to a percentage of
#'    uncertainty. When using \code{fan=TRUE} the raw count data should be used.
#'    In order to enable the calculation of the coefficient of variation a
#'    pseduo-count of 1 is added to every count.
#'
#' @return A \code{\link{ggplot2}{ggplot}} object.
#' @import SingleCellExperiment
#' @import ggplot2
#' @importFrom dplyr as_tibble
#' @importFrom methods slotNames
#' @importFrom stats quantile
#' @import grid
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
#' tenx_pbmc3k <- make_hexbin(tenx_pbmc3k, 80, dimension_reduction = "PCA")
#' plot_hexbin_bivariate(tenx_pbmc3k, type = "counts", feature = "ENSG00000135250")
#' plot_hexbin_bivariate(tenx_pbmc3k, type = "counts", feature = "ENSG00000135250", fan = TRUE)
plot_hexbin_bivariate <- function(sce,
                                  mod = "RNA",
                                  type,
                                  feature,
                                  fan = FALSE,
                                  title = NULL,
                                  xlab = NULL,
                                  ylab = NULL) {
    out <- .extract_hexbin(sce)
    cID <- .extract_cID(sce)

    if (is.null(out)) {
        stop("Compute hexbin representation before plotting.")
    }

    if (fan) {
        warning("When fan=TRUE the raw count data should be used.")
    }

    x <- .prepare_data_feature(sce, mod, type, feature)

    out <- .plot_hexbin_bivariate_helper_1(x, feature, out, cID, fan)

    .plot_bivariate(out,
        title = title, xlab = xlab, ylab = ylab, fan,
        feature
    )
}
