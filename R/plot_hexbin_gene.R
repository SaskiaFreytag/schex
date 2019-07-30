#' Plot of gene expression of single cells in bivariate hexagon cells.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   or \code{\link[Seurat]{Seurat}} object.
#' @param type A string referring to the type of expression data plotted.
#'     Possible options are \code{counts} for raw counts and \code{logcounts}
#'     for normalized coun  in the
#'     \code{\link[SingleCellExperiment]{SingleCellExperiment}} object,
#'     additionally for the \code{\link[Seurat]{Seurat}} object
#'     the scaled data can be accessed using \code{scale.data}.
#' @param gene A string referring to the name of one gene.
#' @param action A strings pecifying how meta data of observations in
#'   binned  hexagon cells are to be summarized. Possible actions are
#'   \code{prop_0}, \code{mode}, \code{mean} and \code{median} (see details).
#' @param title A string containing the title of the plot.
#' @param xlab A string containing the title of the x axis.
#' @param ylab A string containing the title of the y axis.
#'
#' @details This function plots the expression of any gene in the hexagon cell
#'   representation calculated with \code{\link{make_hexbin}}. The chosen gene
#'   expression is summarized by one of four actions \code{prop_0}, \code{mode},
#'   \code{mean} and \code{median}:
#'
#'   \describe{
#'     \item{\code{prop_0}}{Returns the proportion of observations in the bin
#'      greater than 0. The associated meta data column needs to be numeric.}
#'     \item{\code{mode}}{Returns the mode of the observations in the bin. The
#'      associated meta data column needs to be numeric.}
#'     \item{\code{mean}}{Returns the mean of the observations in the bin. The
#'      associated meta data column needs to be numeric.}
#'      \item{\code{median}}{Returns the median of the observations in the bin.
#'      The associated meta data column needs to be numeric.}
#'   }
#'
#' @return A \code{\link{ggplot2}{ggplot}} object.
#' @import Seurat
#' @import SingleCellExperiment
#' @import ggplot2
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
#' rm_ind <- calcAverage(tenx_pbmc3k)<0.1
#' tenx_pbmc3k <- tenx_pbmc3k[!rm_ind,]
#' tenx_pbmc3k <- calculateQCMetrics(tenx_pbmc3k)
#' tenx_pbmc3k <- normalize(tenx_pbmc3k)
#' tenx_pbmc3k <- runPCA(tenx_pbmc3k)
#' tenx_pbmc3k <- make_hexbin( tenx_pbmc3k, 20, dimension_reduction = "PCA")
#' plot_hexbin_gene(tenx_pbmc3k, type="logcounts", gene="ENSG00000135250", action="mean")
#' plot_hexbin_gene(tenx_pbmc3k, type="logcounts", gene="ENSG00000135250", action="mode")
#' }
setGeneric("plot_hexbin_gene", function(sce, type,
                                        gene,
                                        action,
                                        title=NULL,
                                        xlab=NULL,
                                        ylab=NULL) standardGeneric("plot_hexbin_gene"))

#' @export
#' @describeIn plot_hexbin_gene  Plot of gene expression into hexagon cell for
#'   SingleCellExperiment object.
setMethod("plot_hexbin_gene", "SingleCellExperiment", function(sce,
                                                               type,
                                                               gene,
                                                               action,
                                                               title=NULL,
                                                               xlab=NULL,
                                                               ylab=NULL) {

  if(!type %in% c("counts", "logcounts", "scale.data")){
    stop("Specify a valid assay type.")
  }

  out <- sce@metadata$hexbin[[2]]
  cID <- sce@metadata$hexbin[[1]]

  if(is.null(out)){
    stop("Compute hexbin representation before plotting.")
  }

  ind <- match(gene, rownames(sce))

  if (is.na(ind)) {
    stop("Gene cannot be found.")
  }

  if(type=="counts"){
    x <- as.numeric(counts(sce[ind,]))
  }
  if(type=="logcounts"){
    x <- as.numeric(logcounts(sce[ind,]))
  }

  hh <- .make_hexbin_function(x, action, cID)
  out <- as_tibble(out)

  gene <- gsub("-", "_", gene)

  col_hh <- paste0(gene, "_", action)

  func1 <- paste0("out$", col_hh, " <- hh")
  eval(parse(text=func1))

  .plot_hexbin(out, colour_by=col_hh,
               title=title, xlab=xlab, ylab=ylab)

})

#' @export
#' @describeIn plot_hexbin_gene  Plot of gene expression into hexagon cell for
#'   Seurat object.
setMethod("plot_hexbin_gene", "Seurat", function(sce,
                                                 type,
                                                 gene,
                                                 action,
                                                 title=NULL,
                                                 xlab=NULL,
                                                 ylab=NULL) {

  if(!type %in% c("counts", "logcounts", "scale.data")){
    stop("Specify a valid assay type.")
  }

  out <- sce@misc$hexbin[[2]]
  cID <- sce@misc$hexbin[[1]]

  if(is.null(out)){
    stop("Compute hexbin representation before plotting.")
  }

  if(type == "counts"){

    x <- GetAssayData(sce, "counts")

  } else if(type == "logcounts"){

    x <- GetAssayData(sce, "data")

  } else{

    x <- GetAssayData(sce, "scale.data")

  }

  ind <- match(gene, rownames(x))

  if (is.na(ind)) {
    stop("Gene cannot be found.")
  }

  x <- as.numeric(x[ind,])

  hh <- .make_hexbin_function(x, action, cID)
  out <- as_tibble(out)

  gene <- gsub("-", "_", gene)

  col_hh <- paste0(gene, "_", action)

  func1 <- paste0("out$", col_hh, " <- hh")
  eval(parse(text=func1))

  .plot_hexbin(out, colour_by=col_hh,
               title=title, xlab=xlab, ylab=ylab)

})
