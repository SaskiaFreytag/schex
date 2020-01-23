#' Plot of feature expression and uncertainty of single cells in bivariate 
#'    hexagon cells.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'    or \code{\link[Seurat]{Seurat-class}} object.
#' @param mod A string referring to the name of the modality used for plotting.
#'     For RNA modality use "RNA". For other modalities use name of alternative 
#'     object for the \code{\link[SingleCellExperiment]{SingleCellExperiment}} 
#'     or the name of the assay for the \code{\link[Seurat]{Seurat-class}} 
#'     object.
#' @param type A string referring to the type of assay in the
#'    \code{\link[SingleCellExperiment]{SingleCellExperiment}} object or the
#'    data transformation in the \code{\link[Seurat]{Seurat-class}} object.
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
#'    \code{fan=TRUE}, the mean expression of  The colour of the hexagon represents a 
#'    combination of both standard deviation and mean expression.
#'
#'    To access the data that has been integrated in the 
#'    \code{\link[Seurat]{Seurat-class}} object specifiy \code{mod="integrated"}.
#'
#' @return A \code{\link{ggplot2}{ggplot}} object.
#' @import Seurat
#' @import SingleCellExperiment
#' @import ggplot2
#' @importFrom dplyr as_tibble
#' @importFrom methods slotNames
#' @importFrom stats quantile
#' @importFrom cowplot ggdraw draw_plot
#' @export
#'
#' @examples
#' # For Seurat object
#' library(Seurat)
#' data("pbmc_small")
#' pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
#' plot_hexbin_bivariate(pbmc_small, type="counts", feature="CD3D")
#' plot_hexbin_bivariate(pbmc_small, type="counts", feature="CD3D", fan=TRUE)
#' \dontrun{
#' library(TENxPBMCData)
#' library(scater)
#' tenx_pbmc3k <- TENxPBMCData(dataset = "pbmc3k")
#' rm_ind <- calcAverage(tenx_pbmc3k)<0.1
#' tenx_pbmc3k <- tenx_pbmc3k[!rm_ind,]
#' colData(tenx_pbmc3k) <- cbind(colData(tenx_pbmc3k),
#'    perCellQCMetrics(tenx_pbmc3k))
#' tenx_pbmc3k <- logNormCounts(tenx_pbmc3k)
#' tenx_pbmc3k <- runPCA(tenx_pbmc3k)
#' tenx_pbmc3k <- make_hexbin( tenx_pbmc3k, 20, dimension_reduction = "PCA")
#' plot_hexbin_bivariate(tenx_pbmc3k, type="counts",
#'    feature="ENSG00000135250", fan=TRUE)
#' }
plot_hexbin_bivariate <- function(sce, 
                                mod="RNA", 
                                type,
                                feature,
                                fan=FALSE,
                                title=NULL,
                                xlab=NULL,
                                ylab=NULL) {
  
  out <- schex:::.extract_hexbin(sce)
  cID <- schex:::.extract_cID(sce)
  
  if(is.null(out)){
    stop("Compute hexbin representation before plotting.")
  }
  
  if(fan){
    warning("When fan=TRUE the raw count data should be used.")
  }
  
  x <- schex:::.prepare_data_feature(sce, mod, type, feature) 
  
  out <- .plot_hexbin_bivariate_helper_1(x, feature, out, cID, fan)
  
  .plot_bivariate(out, title=title, xlab=xlab, ylab=ylab, fan,
                              feature)
}





