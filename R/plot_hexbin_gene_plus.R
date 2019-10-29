#' Plot of gene expression and meta data of single cell data in
#' bivariate hexagon cells.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   or \code{\link[Seurat]{Seurat}} object.
#' @param col A string referring to the name of one column in the meta data of
#'   sce by which to make the outlines. Note that this should be a factor or
#'   a character.
#' @param type A string referring to the type of assay in the
#'   \code{\link[SingleCellExperiment]{SingleCellExperiment}} object or the
#'   data transformation in the \code{\link[Seurat]{Seurat}} object.
#' @param gene A string referring to the name of one gene.
#' @param action A string specifying how gene expression of observations in
#'   binned  hexagon cells are to be summarized. Possible actions are
#'  \code{prop_0}, \code{mode}, \code{mean} and
#'   \code{median} (see details).
#' @param colors A vector of strings specifying which colors to use for plotting
#'    the different levels in the selected column of the meta data.
#' @param title A string containing the title of the plot.
#' @param xlab A string containing the title of the x axis.
#' @param ylab A string containing the title of the y axis.
#' @param exapnd_hull A numeric value determining the expansion of the line
#'   marking different clusters.
#' @param ... Additional arguments passed on to
#'    \code{\link[geom_mark_hull]{ggforce}}.
#'
#' @details This function plots any gene expresssion in the hexagon cell
#'   representation calculated with \code{\link{make_hexbin}} as well as at the
#'   same time representing outlines of clusters. The chosen gene
#'   expression is summarized by one of four actions \code{prop_0},
#'   \code{mode}, \code{mean} and \code{median}:
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
#' @importFrom dplyr as_tibble
#' @importFrom ggforce geom_mark_hull
#' @export
#'
#' @examples
#' #' # For Seurat object
#' library(Seurat)
#' data("pbmc_small")
#' pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
#' plot_hexbin_gene_plus(pbmc_small, col="RNA_snn_res.1", type="counts",
#' gene="NRBP1",action="mean")
setGeneric("plot_hexbin_gene_plus", function(sce,
                                        col,
                                        type,
                                        gene,
                                        action,
                                        colors=NULL,
                                        title=NULL,
                                        xlab=NULL,
                                        ylab=NULL,
                                        expand_hull=3,
                                        ...) standardGeneric("plot_hexbin_gene_plus"))

#' @export
#' @describeIn plot_hexbin_gene_plus  Plot of gene expression and meta data
#'   of single cell data into hexagon cell for  SingleCellExperiment object.
setMethod("plot_hexbin_gene_plus", "SingleCellExperiment", function(sce,
                                                               col,
                                                               type,
                                                               gene,
                                                               action,
                                                               colors=NULL,
                                                               title=NULL,
                                                               xlab=NULL,
                                                               ylab=NULL,
                                                               expand_hull=3,
                                                               ...) {

  out <- sce@metadata$hexbin[[2]]
  cID <- sce@metadata$hexbin[[1]]

  if(is.null(out)){
    stop("Compute hexbin representation before plotting.")
  }

  if (any(!col %in% colnames(colData(sce)))) {
    stop("Column cannot be found in colData(sce).")
  }

  ind <- match(gene, rownames(sce))
  x <- assays(sce)

  if(!type %in% names(x)){
    stop("Specify a valid assay type.")
  }

  if (is.na(ind)) {
    stop("Gene cannot be found.")
  }

  x_gene <- as.numeric(x[[which(names(x)==type)]][ind,])

  hh_gene <- schex:::.make_hexbin_function(x_gene, action, cID)

  name_s <- paste0("sce$", col)
  func <- paste0("x <- ", name_s)

  eval(parse(text = func))

  hh <- schex:::.make_hexbin_function(x, 'majority', cID)
  out <- as_tibble(out)

  col_hh <-paste0(col, "_", "majority")

  if(is.factor(x)){
    func1 <- paste0("out$", col_hh, " <- factor(hh, levels=",
                        "levels(x))")
  } else {
        func1 <- paste0("out$", col_hh, " <- hh")
  }

  eval(parse(text=func1))

  if(grepl("^[[:digit:]]", gene )){
    gene <- paste0("G_", gene)
  }

  gene <- gsub("-", "_", gene)

  col_hh_gene <- paste0(gene, "_", action)

  func2 <- paste0("out$", col_hh_gene, " <- hh_gene")
  eval(parse(text=func2))

  .plot_hexbin_plus(out, colour_by = col_hh, fill_by_gene = col_hh_gene,
                    colors=colors, expand_hull=expand_hull, title=title,
                    xlab=xlab, ylab=ylab, ...)

})

#' @export
#' @describeIn plot_hexbin_gene_plus Plot of gene expression and meta data of
#'   single cell data into hexagon cell for Seurat object.
setMethod("plot_hexbin_gene_plus", "Seurat", function(sce,
                                                 col,
                                                 type,
                                                 gene,
                                                 action,
                                                 colors=NULL,
                                                 title=NULL,
                                                 xlab=NULL,
                                                 ylab=NULL,
                                                 expand_hull=3,
                                                 ...) {

  out <- sce@misc$hexbin[[2]]
  cID <- sce@misc$hexbin[[1]]

  if(is.null(out)){
    stop("Compute hexbin representation before plotting.")
  }

  if (any(!col %in% colnames(sce@meta.data))) {
    stop("Column cannot be found in slot(sce, 'meta.data').")
  }

  x_gene <- GetAssayData(sce, type)

  ind <- match(gene, rownames(x_gene))

  if (is.na(ind)) {
    stop("Gene cannot be found.")
  }

  x_gene <- as.numeric(x_gene[ind,])

  hh_gene <- schex:::.make_hexbin_function(x_gene, action, cID)

  name_s <- paste0("sce$", col)
  func <- paste0("x <- ", name_s)

  eval(parse(text = func))

  hh <- schex:::.make_hexbin_function(x, "majority", cID)
  out <- as_tibble(out)

  col_hh <-paste0(col, "_", "majority")

  if(is.factor(x)){
    func1 <- paste0("out$", col_hh, " <- factor(hh, levels=",
                    "levels(x))")
  } else {
    func1 <- paste0("out$", col_hh, " <- hh")
  }

  eval(parse(text=func1))

  if(grepl("^[[:digit:]]", gene )){
    gene <- paste0("G_", gene)
  }

  gene <- gsub("-", "_", gene)

  col_hh_gene <- paste0(gene, "_", action)

  func2 <- paste0("out$", col_hh_gene, " <- hh_gene")
  eval(parse(text=func2))

  .plot_hexbin_plus(out, colour_by = col_hh, fill_by_gene = col_hh_gene,
                    colors=colors, expand_hull=expand_hull, title=title,
                    xlab=xlab, ylab=ylab, ...)

})

.plot_hexbin_plus <- function(drhex, colour_by="Cluster_majority", fill_by_gene,
                            colors=NULL, expand_hull=3, legend=legend,
                         title=NULL, xlab=NULL, ylab=NULL, ...) {

  if (any(!c("x", "y", colour_by) %in% colnames(drhex))) {
    stop("The dataframe must contain columns named 'x', 'y' and label.")
  }

  if(is.null(title)) {
    title <- colour_by
  }

  if(is.null(xlab)) {
    xlab <- "x"
  }

  if(is.null(ylab)) {
    ylab <- "y"
  }

  if(is.null(colors)){

    ggplot(drhex, aes_string(x="x", y="y", fill=fill_by_gene)) +
      geom_hex(stat = "identity") +
      geom_mark_hull(aes_string(label = colour_by, col = colour_by),
      show.legend = FALSE, expand = unit(expand_hull, "mm"),
      fill=NA, size=2, ...) +
      theme_classic() + scale_fill_viridis_c() +
      ggtitle(title) + labs(x=xlab, y=ylab) +
      theme(legend.title=element_blank())

    } else {

      ggplot(drhex, aes_string(x="x", y="y", fill=fill_by_gene)) +
        geom_hex(stat = "identity") +
        geom_mark_hull(aes_string(label = colour_by, col = colour_by),
        show.legend = FALSE, expand = unit(expand_hull, "mm"),
        fill=NA, size=2, ...) + theme_classic() + scale_fill_viridis_c() +
        ggtitle(title) + labs(x=xlab, y=ylab) +
        theme(legend.title=element_blank()) + scale_color_manual(values=colors)

    }

}

