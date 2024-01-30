#' Bivariate binning of single cell data into hexagon cells.
#'
#' \code{make_hexbin} returns a
#'    \code{\link[SingleCellExperiment]{SingleCellExperiment}} object of binned hexagon cells.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @param nbins The number of bins partitioning the range of the first
#'    component of the chosen dimension reduction.
#' @param dimension_reduction A string indicating the reduced dimension
#'    result to calculate hexagon cell representation of.
#' @param use_dims A vector of two integers specifying the dimensions used.
#'
#' @details This function bins observations with computed reduced dimension
#'    results into hexagon cells. For a
#'    \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'    as a list in the \code{@metadata}. The list contains two items. The first
#'    item stores a vector specifying the hexagon ID for each
#'    observation. The second item stores a matrix with the x and y positions of
#'    the hexagon cells and the number of observations in each of them.
#'
#' @return A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @importFrom hexbin hexbin hcell2xy
#' @import SingleCellExperiment
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
setGeneric("make_hexbin", function(sce, nbins = 80,
                                   dimension_reduction = "UMAP",
                                   use_dims = c(1, 2)) {
    standardGeneric("make_hexbin")
})

#' @export
#' @describeIn make_hexbin Bivariate binning of SingleCellExperiment
#'   into hexagon cells.
setMethod("make_hexbin", "SingleCellExperiment", function(sce,
                                                          nbins = 80,
                                                          dimension_reduction = "UMAP",
                                                          use_dims = c(1, 2)) {
    if (!is.element(dimension_reduction, reducedDimNames(sce))) {
        stop("Specify existing dimension reduction.")
    }

    dr <- reducedDim(sce, dimension_reduction)

    res <- .make_hexbin_helper(dr, nbins, use_dims)
    sce@metadata$hexbin <- res

    return(sce)
})

.make_hexbin_helper <- function(dr, nbins = 80, use_dims) {
    if (dim(dr)[2] < max(use_dims)) {
        stop("Please specify use_dims that are calculated.")
    }

    xbnds <- range(c(dr[, use_dims[1]]))
    ybnds <- range(c(dr[, use_dims[2]]))

    drhex <- hexbin(dr[, use_dims[1]],
        dr[, use_dims[2]],
        nbins,
        xbnds = xbnds,
        ybnds = ybnds,
        IDs = TRUE
    )
    cID <- drhex@cID
    drhex <- cbind(
        as.numeric(hcell2xy(drhex)$x),
        as.numeric(hcell2xy(drhex)$y),
        as.numeric(drhex@count)
    )

    colnames(drhex) <- c("x", "y", "number_of_cells")

    res <- list(cID = cID, hexbin.matrix = drhex)

    return(res)
}
