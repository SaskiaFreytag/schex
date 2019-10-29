#' Plot of external feature expression of single cells in bivariate hexagon
#'    cells.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'    or \code{\link[Seurat]{Seurat-class}} object.
#' @param mod A string referring to the name of the alternative object in a
#'    \code{\link[SingleCellExperiment]{SingleCellExperiment}} or the assay in a
#'    \code{\link[Seurat]{Seurat-class}} object that stores the protein information.
#' @param type A string referring to the type of assay in the
#'    \code{\link[SingleCellExperiment]{SingleCellExperiment}} object or the
#'    data transformation in the \code{\link[Seurat]{Seurat-class}} object.
#' @param feature A string referring to the name of one external feature.
#' @param action A strings pecifying how meta data of observations in
#'    binned  hexagon cells are to be summarized. Possible actions are
#'    \code{prop_0}, \code{mode}, \code{mean} and \code{median} (see details).
#' @param title A string containing the title of the plot.
#' @param xlab A string containing the title of the x axis.
#' @param ylab A string containing the title of the y axis.
#'
#' @details This function plots the expression of any feature in the hexagon
#'    cell representation calculated with \code{\link{make_hexbin}}. The chosen
#'    gene expression is summarized by one of four actions \code{prop_0},
#'    \code{mode}, \code{mean} and \code{median}:
#'
#'    \describe{
#'      \item{\code{prop_0}}{Returns the proportion of observations in the bin
#'       greater than 0. The associated meta data column needs to be numeric.}
#'      \item{\code{mode}}{Returns the mode of the observations in the bin. The
#'       associated meta data column needs to be numeric.}
#'      \item{\code{mean}}{Returns the mean of the observations in the bin. The
#'       associated meta data column needs to be numeric.}
#'       \item{\code{median}}{Returns the median of the observations in the bin.
#'       The associated meta data column needs to be numeric.}
#'    }
#'
#' @return A \code{\link{ggplot2}{ggplot}} object.
#' @import Seurat
#' @import SingleCellExperiment
#' @import ggplot2
#' @importFrom dplyr as_tibble
#' @importFrom methods slotNames
#' @export
#'
#' @examples
#' # For Seurat object
#' library(Seurat)
#' data("pbmc_small")
#' pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
#' protein <- matrix(rnorm(10* ncol(pbmc_small)), ncol=ncol(pbmc_small))
#' rownames(protein) <- paste0("A", seq(1,10,1))
#' colnames(protein) <- colnames(pbmc_small)
#' pbmc_small[["ADT"]] <- CreateAssayObject(counts = protein)
#' pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
#' plot_hexbin_feature(pbmc_small, type="counts", mod="ADT",
#'     feature="A1", action="prop_0")
setGeneric("plot_hexbin_feature", function(sce, 
    mod, 
    type,
    feature,
    action,
    title=NULL,
    xlab=NULL,
    ylab=NULL) standardGeneric("plot_hexbin_feature"))

#' @export
#' @describeIn plot_hexbin_feature  Plot of gene expression into hexagon
#'   cell for SingleCellExperiment object.
setMethod("plot_hexbin_feature", "SingleCellExperiment", function(sce,
    mod,
    type,
    feature,
    action,
    title=NULL,
    xlab=NULL,
    ylab=NULL) {
  
   out <- sce@metadata$hexbin[[2]]
   cID <- sce@metadata$hexbin[[1]]
  
   if(is.null(out)){
       stop("Compute hexbin representation before plotting.")
    }
  
    x <-.prepare_data_feature(sce, mod, type, feature)
  
    .plot_hexbin_feature_helper(x, feature, out, cID, action, title,
        xlab, ylab)
  
})

#' @export
#' @describeIn plot_hexbin_feature  Plot of gene expression into hexagon cell
#'    for Seurat object.
setMethod("plot_hexbin_feature", "Seurat", function(sce,
    mod,
    type,
    feature,
    action,
    title=NULL,
    xlab=NULL,
    ylab=NULL) {
  
    out <- sce@misc$hexbin[[2]]
    cID <- sce@misc$hexbin[[1]]
    
    if(is.null(out)){
       stop("Compute hexbin representation before plotting.")
    }
  
    x <-.prepare_data_feature(sce, mod, type, feature)
  
    .plot_hexbin_feature_helper(x, feature, out, cID, action, title,
        xlab, ylab)
  
})

.plot_hexbin_feature_helper <- function(x, feature, out, cID, action, title,
    xlab, ylab){
  
    hh <- .make_hexbin_function(x, action, cID)
    out <- as_tibble(out)
  
    if(grepl("^[[:digit:]]", feature )){
        feature <- paste0("F_", feature)
    }
  
    feature <- gsub("-", "_", feature)
  
    col_hh <- paste0(feature, "_", action)
  
    func1 <- paste0("out$", col_hh, " <- hh")
    eval(parse(text=func1))
  
    .plot_hexbin(out, colour_by=col_hh,
        title=title, xlab=xlab, ylab=ylab)
}
