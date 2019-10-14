#' Plot of density of observations from single cell data
<<<<<<< HEAD
#'    in bivariate hexagon cells.
=======
#'   in bivariate hexagon cells.
>>>>>>> 73fca473510fbfc1906bc0df2320441b4042cb9b
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'    or \code{\link[Seurat]{Seurat}} object.
#' @param title A string containing the title of the plot.
#' @param xlab A string containing the title of the x axis.
#' @param ylab A string containing the title of the y axis.
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
#' plot_hexbin_density(pbmc_small)
#' \dontrun{
#' library(TENxPBMCData)
#' library(scater)
#' tenx_pbmc3k <- TENxPBMCData(dataset = "pbmc3k")
#' rm_ind <- calcAverage(tenx_pbmc3k)<0.1
#' tenx_pbmc3k <- tenx_pbmc3k[!rm_ind,]
#' tenx_pbmc3k <- calculateQCMetrics(tenx_pbmc3k)
#' tenx_pbmc3k <- normalize(tenx_pbmc3k)
#' tenx_pbmc3k <- runPCA(tenx_pbmc3k)
#' tenx_pbmc3k <- make_hexbin( tenx_pbmc3k, 20, dimension_reduction = "PCA")
#' plot_hexbin_density(tenx_pbmc3k)
#' }
setGeneric("plot_hexbin_density", function(sce, title=NULL,
    xlab=NULL,
    ylab=NULL)
    standardGeneric("plot_hexbin_density"))

#' @export
#' @describeIn plot_hexbin_density  Plot of cell density in hexagon cell for
#'    SingleCellExperiment object.
setMethod("plot_hexbin_density", "SingleCellExperiment", function(sce,
    title=NULL,
    xlab=NULL,
    ylab=NULL) {
  
    out <- sce@metadata$hexbin[[2]]
  
    .plot_hexbin_density_helper(out, title, xlab, ylab)

})

#' @export
#' @describeIn plot_hexbin_density  Plot of cell density in hexagon cell for
#'   Seurat object.
setMethod("plot_hexbin_density", "Seurat", function(sce,
    title=NULL,
    xlab=NULL,
    ylab=NULL) {
  
    out <- sce@misc$hexbin[[2]]
  
    .plot_hexbin_density_helper(out, title, xlab, ylab)
  
})

.plot_hexbin_density_helper <- function(out, title, xlab, ylab){
  
    if(is.null(out)){
      stop("Compute hexbin representation before plotting.")
   }
  
    if(is.null(title)) {
        title <- "Density"
    }
  
    if(is.null(xlab)) {
        xlab <- "x"
    }
    
    if(is.null(ylab)) {
        ylab <- "y"
    }
  
    out <- as_tibble(out)
  
    ggplot(out, aes_string("x", "y", fill="number_of_cells")) +
        geom_hex(stat = "identity") + scale_fill_viridis_c() +
        theme_classic() + theme(legend.position="bottom") + ggtitle(title) +
        labs(x=xlab, y=ylab) + theme(legend.title=element_blank())
  
}