#' Plot of interaction of expression of single cells in bivariate hexagon cells.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   or \code{\link[Seurat]{Seurat}} object.
#' @param mod A vector of strings referring to the names of the modularities. For
#'     \code{\link[SingleCellExperiment]{SingleCellExperiment}} use "RNA" to access
#'     the RNA expression data stored as the main experiment type.
#' @param type A vector of strings referring to the types of assays in the
#'     \code{\link[SingleCellExperiment]{SingleCellExperiment}} or the types of
#'     transformation.
#' @param feature A vector of strings referring to the names of one features in
#'     the same order as the vector of modularities.
#' @param interact A string specifying how interaction between features is
#'   calculated. Possible interaction measures are
#'   \code{corr_spearman} and \code{mi} (see details).
#' @param title A string containing the title of the plot.
#' @param xlab A string containing the title of the x axis.
#' @param ylab A string containing the title of the y axis.
#'
#' @details This function plots the interaction between any features in the hexagon cell
#'   representation calculated with \code{\link{make_hexbin}}. The interaction
#'   between the chosen features is calculated by one of two measurers \code{corr_spearman},
#'   and \code{mi}:
#'
#'   \describe{
#'     \item{\code{mi}}{Returns the mutual information coefficient.}
#'     \item{\code{corr_spearman}}{Returns the Spearman correlation.}
#'   }
#'
#' @return A \code{\link{ggplot2}{ggplot}} object.
#' @import Seurat
#' @import SingleCellExperiment
#' @importFrom entropy mi.plugin
#' @import ggplot2
#' @importFrom dplyr as_tibble
#' @export
#'
#' @examples
#' # For Seurat object
#' library(Seurat)
#' data("pbmc_small")
#' pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
#' plot_hexbin_gene(pbmc_small, type="counts", gene="TALDO1", action="prop_0")
#' # For SingleCellExperiment object
#' \dontrun{
#' library(TENxPBMCData)
#' library(scater)
#' tenx_pbmc3k <- TENxPBMCData(dataset = "pbmc3k")
#' rm_ind <- calculateAverage(tenx_pbmc3k)<0.1
#' tenx_pbmc3k <- tenx_pbmc3k[!rm_ind,]
#' colData(tenx_pbmc3k) <- cbind(colData(tenx_pbmc3k),
#'      perCellQCMetrics(tenx_pbmc3k))
#' tenx_pbmc3k <- normalize(tenx_pbmc3k)
#' tenx_pbmc3k <- runPCA(tenx_pbmc3k)
#' tenx_pbmc3k <- make_hexbin( tenx_pbmc3k, 20, dimension_reduction = "PCA")
#' plot_hexbin_gene(tenx_pbmc3k, type="logcounts", gene="ENSG00000135250", action="mean")
#' plot_hexbin_gene(tenx_pbmc3k, type="logcounts", gene="ENSG00000135250", action="mode")
#' }
setGeneric("plot_hexbin_interact", function(sce, mod,
                                            type,
                                            feature,
                                           interact,
                                           title=NULL,
                                           xlab=NULL,
                                           ylab=NULL) standardGeneric("plot_hexbin_interact"))

#' @export
#' @describeIn plot_hexbin_interact  Plot of gene expression into hexagon cell for
#'   SingleCellExperiment object.
setMethod("plot_hexbin_interact", "SingleCellExperiment", function(sce,
                                                                   mod,
                                                                   type,
                                                                   feature,
                                                                   interact,
                                                                   title=NULL,
                                                                   xlab=NULL,
                                                                   ylab=NULL) {

  if(length(mod)!=length(feature)|length(feature)!=length(types)){
    stop("Specify the same number of modularities, types and features.")
  }


  if(!mod[1] %in% c(altExpNames(sce), "RNA")|
     !mod[2] %in% c(altExpNames(sce), "RNA")){
    stop("Specify a valid modularity.")
  }

  if(mod[1]!="RNA"){
    first_x <- assays(altExp(sce, mod[1]))
  } else {
    first_x <- assays(sce)
  }

  if(mod[2]!="RNA"){
    second_x <- assays(altExp(sce, mod[2]))
  } else {
    second_x <- assays(sce)
  }

  if(!type[1] %in% names(first_x)|
     !type[2] %in% names(second_x)){
       stop("Specify a valid assay type.")
     }

     out <- sce@metadata$hexbin[[2]]
     cID <- sce@metadata$hexbin[[1]]

     if(is.null(out)){
       stop("Compute hexbin representation before plotting.")
     }

     first_x <- first_x[[which(names(first_x)==type[1])]]
     second_x <- second_x[[which(names(second_x)==type[2])]]

     first_ind <- match(feature[1], rownames(first_x))
     second_ind <- match(feature[2], rownames(second_x))

     if (is.na(second_ind)|is.na(first_ind)) {
       stop("One/two features cannot be found.")
     }

     first_x <- as.numeric(first_x[first_ind,])
     second_x <- as.numeric(second_x[second_ind,])

     hh <- .interact_hexbin_function(first_x, second_x, interact, cID)
     out <- as_tibble(out)

     feature <- gsub("-", "_", feature)

     col_hh <- paste0(interact, "_", feature[1], "_", feature[2])

     func1 <- paste0("out$", col_hh, " <- hh")
     eval(parse(text=func1))

     schex:::.plot_hexbin(out, colour_by=col_hh,
                          title=title, xlab=xlab, ylab=ylab)

})

#' @export
#' @describeIn plot_hexbin_interact  Plot of gene expression into hexagon cell for
#'   Seurat object.
setMethod("plot_hexbin_interact", "Seurat", function(sce,
                                                    mod,
                                                    type,
                                                    feature,
                                                    interact,
                                                    title=NULL,
                                                    xlab=NULL,
                                                    ylab=NULL) {

  if(length(mod)!=length(feature)|length(feature)!=length(types)){
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

  first_x <- GetAssayData(sce, assay=mod[1], type[1])
  second_x <- GetAssayData(sce, assay=mod[2], type[2])

  first_ind <- match(feature[1], rownames(first_x))
  second_ind <- match(feature[2], rownames(second_x))

  if (is.na(second_ind)|is.na(first_ind)) {
    stop("One/two features cannot be found.")
  }

  first_x <- as.numeric(first_x[first_ind,])
  second_x <- as.numeric(second_x[second_ind,])

  hh <- .interact_hexbin_function(first_x, second_x, interact, cID)
  out <- as_tibble(out)

  feature <- gsub("-", "_", feature)

  col_hh <- paste0(interact, "_", feature[1], "_", feature[2])

  func1 <- paste0("out$", col_hh, " <- hh")
  eval(parse(text=func1))

  schex:::.plot_hexbin(out, colour_by=col_hh,
                       title=title, xlab=xlab, ylab=ylab)

})


.interact_hexbin_function<- function(first_x, second_x, interact, cID) {

  if(interact=="corr_spearman"){

    func_if <- !(is.numeric(first_x)|is.numeric(second_x))

    if (func_if) {
      stop(paste0("Features need to be numeric."))

    } else {

      res_first <- tapply(first_x, cID, FUN = function(z) z)
      res_second <- tapply(second_x, cID, FUN = function(z) z)

      res <- unlist(lapply(1:length(res_first), function(x)
        cor(res_first[[x]], res_second[[x]], method="spearman")))

     return(res)
    }
  }

  if(interact=="mi"){

    func_if <- !(is.numeric(first_x)|is.numeric(second_x))

    if (func_if) {
      stop(paste0("Features need to be numeric."))

    } else {

      res_first <- tapply(first_x, cID, FUN = function(z) z)
      res_second <- tapply(second_x, cID, FUN = function(z) z)

      res <- lapply(1:length(res_first), function(x)
        rbind(res_first[[x]], res_second[[x]]))

      res <- unlist(lapply(res, function(x)
        mi.plugin(x)))

     return(res)
    }
  }
}
