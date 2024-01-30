#' Plot of density of observations from single cell data
#'    in bivariate hexagon cells.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @param title A string containing the title of the plot.
#' @param xlab A string containing the title of the x axis.
#' @param ylab A string containing the title of the y axis.
#'
#' @return A \code{\link{ggplot2}{ggplot}} object.
#' @import SingleCellExperiment
#' @import ggplot2
#' @importFrom dplyr as_tibble
#' @import rlang
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
#' plot_hexbin_density(tenx_pbmc3k)
plot_hexbin_density <- function(
    sce,
    title = NULL,
    xlab = NULL,
    ylab = NULL) {
    out <- .extract_hexbin(sce)
    .plot_hexbin_density_helper(out, title, xlab, ylab)
}

.plot_hexbin_density_helper <- function(out, title, xlab, ylab) {
    if (is.null(out)) {
        stop("Compute hexbin representation before plotting.")
    }

    if (is.null(title)) {
        title <- "Density"
    }

    if (is.null(xlab)) {
        xlab <- "x"
    }

    if (is.null(ylab)) {
        ylab <- "y"
    }

    out <- as_tibble(out)

    ggplot(out, aes(x = !!sym("x"), y = !!sym("y"), fill = !!sym("number_of_cells"))) +
        geom_hex(stat = "identity") +
        scale_fill_viridis_c() +
        theme_classic() +
        theme(legend.position = "bottom") +
        ggtitle(title) +
        labs(x = xlab, y = ylab) +
        theme(legend.title = element_blank())
}
