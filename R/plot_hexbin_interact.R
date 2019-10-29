#' Plot of interaction of expression of single cells in bivariate hexagon cells.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'    or \code{\link[Seurat]{Seurat-class}} object.
#' @param mod A vector of strings referring to the names of the modularities.
#'    For \code{\link[SingleCellExperiment]{SingleCellExperiment}} use "RNA" to
#'    access the RNA expression data stored as the main experiment type.
#' @param type A vector of strings referring to the types of assays in the
#'    \code{\link[SingleCellExperiment]{SingleCellExperiment}} or the types of
#'    transformation in \code{\link[Seurat]{Seurat-class}} object.
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
#'    measurers \code{corr_spearman}, and \code{mi}:
#'
#'    \describe{
#'       \item{\code{mi}}{Returns the mutual information coefficient.}
#'       \item{\code{corr_spearman}}{Returns the Spearman correlation.}
#'    }
#'
#' @return A \code{\link{ggplot2}{ggplot}} object.
#' @import Seurat
#' @import SingleCellExperiment
#' @importFrom entropy mi.plugin
#' @importFrom stats cor
#' @importFrom methods slotNames
#' @import ggplot2
#' @importFrom dplyr as_tibble
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
#' plot_hexbin_interact(pbmc_small, type=c("counts", "counts"),
#'     mod=c("RNA", "ADT" ), feature=c("CD7", "A1"), interact="mi")
setGeneric("plot_hexbin_interact", function(sce,
    mod,
    type,
    feature,
    interact,
    title=NULL,
    xlab=NULL,
    ylab=NULL) standardGeneric("plot_hexbin_interact"))

#' @export
#' @describeIn plot_hexbin_interact  Plot of gene expression into hexagon cell
#'   for SingleCellExperiment object.
setMethod("plot_hexbin_interact", "SingleCellExperiment", function(sce,
    mod,
    type,
    feature,
    interact,
    title=NULL,
    xlab=NULL,
    ylab=NULL) {
  
    if(length(mod)!=length(feature)|length(feature)!=length(type)){
        stop("Specify the same number of modularities, types and features.")
    }
  
    out <- sce@metadata$hexbin[[2]]
    cID <- sce@metadata$hexbin[[1]]
  
    if(is.null(out)){
        stop("Compute hexbin representation before plotting.")
    }  
  
    if(mod[1]!="RNA"){
        first_x <- .prepare_data(sce, "RNA", type, gene)
    } else {
        first_x <- .prepare_data(sce, mod, type, gene)
    }
  
    if(mod[2]!="RNA"){
        second_x <- .prepare_data(sce, "RNA", type, gene)
    } else {
        second_x <- .prepare_data(sce, mod, type, gene)
    }
  
    .plot_hexbin_interact_helper(first_x, second_x, out, cID, interact,
        feature, title, xlab, ylab)
})

#' @export
#' @describeIn plot_hexbin_interact  Plot of gene expression into hexagon cell
#'   for Seurat object.
setMethod("plot_hexbin_interact", "Seurat", function(sce,
    mod,
    type,
    feature,
    interact,
    title=NULL,
    xlab=NULL,
    ylab=NULL) {
  
    if(length(mod)!=length(feature)|length(feature)!=length(type)){
        stop("Specify the same number of modularities, types and features.")
    }
  
    if(!mod[1] %in% names(sce)|!mod[2] %in% names(sce)){
        stop("Specify a valid modularity.")
    }
  
    if(!type[1] %in% slotNames(GetAssay(sce, mod[1]))|
        !type[2] %in% slotNames(GetAssay(sce, mod[2]))){
        stop("Specify a valid assay type.")
    }
  
    out <- sce@misc$hexbin[[2]]
    cID <- sce@misc$hexbin[[1]]
  
    if(is.null(out)){
        stop("Compute hexbin representation before plotting.")
    }
  
    if(mod[1]!="RNA"){
      first_x <- .prepare_data(sce, "RNA", type, gene)
    } else {
      first_x <- .prepare_data(sce, mod, type, gene)
    }
    
    if(mod[2]!="RNA"){
      second_x <- .prepare_data(sce, "RNA", type, gene)
    } else {
      second_x <- .prepare_data(sce, mod, type, gene)
    }
  
    .plot_hexbin_interact_helper(first_x, second_x, out, cID, interact,
        feature, title, xlab, ylab)
})


.plot_hexbin_interact_helper <- function(first_x, second_x,  out, cID, interact,
    feature, title, xlab, ylab) {
  
    hh <- .interact_hexbin_function(first_x, second_x, interact, cID)
    out <- as_tibble(out)
  
    if(any(grepl("^[[:digit:]]", feature))){
        feature <- paste0("F_", feature)
    }
  
    feature <- gsub("-", "_", feature)
  
    col_hh <- paste0(interact, "_", feature[1], "_", feature[2])
  
    func1 <- paste0("out$", col_hh, " <- hh")
    eval(parse(text=func1))
  
    .plot_hexbin(out, colour_by=col_hh,
        title=title, xlab=xlab, ylab=ylab)
}
