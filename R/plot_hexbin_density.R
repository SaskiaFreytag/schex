#' Plot of density of observations from single cell data
#'    in bivariate hexagon cells.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'    or \code{\link[Seurat]{Seurat-class}} object.
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
plot_hexbin_density <- function(sce, 
    title=NULL,
    xlab=NULL,
    ylab=NULL){
    
    out <- .extract_hexbin(sce)
    .plot_hexbin_density_helper(out, title, xlab, ylab)

}

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
  
    ggplot(out, aes("x", "y", fill="number_of_cells")) +
        geom_hex(stat = "identity") + scale_fill_viridis_c() +
        theme_classic() + theme(legend.position="bottom") + ggtitle(title) +
        labs(x=xlab, y=ylab) + theme(legend.title=element_blank())
  
}
