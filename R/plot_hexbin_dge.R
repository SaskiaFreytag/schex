#' Plot of fold change of selected gene in single cell data using 
#'    bivariate hexagon cells.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'    or \code{\link[Seurat]{Seurat-class}} object.
#' @param col A string referring to the name of one column in the meta data of
#'    sce by which to compare. Note this factor can only contain two levels.
#' @param gene A string referring to the name of one gene.
#' @param title A string containing the title of the plot.
#' @param xlab A string containing the title of the x axis.
#' @param ylab A string containing the title of the y axis.
#'
#' @details This function plots differential gene expression within each 
#'   hexagon, which are calculated with \code{\link{make_hexbin}}. 
#'   Note that the fold change is only accurate if the condition 
#'   investigated is within the same individual. For conditions across 
#'   different individuals different methods that account for 
#'   individual-specific effects are required.
#'   
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
#' plot_hexbin_meta(pbmc_small, col="RNA_snn_res.1", action="prop", no=1)
#' plot_hexbin_meta(pbmc_small, col="RNA_snn_res.1", action="prop", no=2)
#' # For SingleCellExperiment object
#' \dontrun{
#' library(TENxPBMCData)
#' library(scater)
#' tenx_pbmc3k <- TENxPBMCData(dataset = "pbmc3k")
#' rm_ind <- calculateAverage(tenx_pbmc3k)<0.1
#' tenx_pbmc3k <- tenx_pbmc3k[-rm_ind,]
#' colData(tenx_pbmc3k) <- cbind(colData(tenx_pbmc3k),
#'      perCellQCMetrics(tenx_pbmc3k))
#' tenx_pbmc3k <- normalize(tenx_pbmc3k)
#' tenx_pbmc3k <- runPCA(tenx_pbmc3k)
#' tenx_pbmc3k <- make_hexbin( tenx_pbmc3k, 20, dimension_reduction = "PCA")
#' plot_hexbin_meta(tenx_pbmc3k, col="total", action="median")
#' }
setGeneric("plot_hexbin_fc", function(sce,
    col,
    gene,
    dge = "fold-change",
    method = "t-test",
    title=NULL,
    xlab=NULL,
    ylab=NULL) standardGeneric("plot_hexbin_fc"))

#' @export
#' @describeIn plot_hexbin_fc  Plot of fold change for specific gene in 
#'   SingleCellExperiment object.
setMethod("plot_hexbin_fc", "SingleCellExperiment", function(sce,
    col,
    gene,
    title=NULL,
    xlab=NULL,
    ylab=NULL) {
  
    out <- sce@metadata$hexbin[[2]]
    cID <- sce@metadata$hexbin[[1]]
  
    if(is.null(out)){
        stop("Compute hexbin representation before plotting.")
    }
  
    if (any(!col %in% colnames(colData(sce)))) {
        stop("Column cannot be found in colData(sce).")
    }
  
    name_s <- paste0("sce$", col)
    func <- paste0("x <- ", name_s)
  
    eval(parse(text = func))
  
    ind <- match(gene, rownames(sce))
    x_gene <- assays(sce)
  
    if(!"logcounts" %in% names(x_gene)){
        stop("Logcounts not found.")
    }
  
    if (is.na(ind)) {
        stop("Gene cannot be found.")
    }
  
    x_gene <- as.numeric(x[[which(names(x_gene)=="logcounts")]][ind,])
  
    .plot_hexbin_fc_helper(x, x_gene, out, cID, col, title, 
         xlab, ylab, colors)
  
})

#' @export
#' @describeIn plot_hexbin_fc  Plot of fold change for specific gene in
#'   Seurat object.
setMethod("plot_hexbin_fc", "Seurat", function(sce,
    col,
    gene,
    title=NULL,
    xlab=NULL,
    ylab=NULL) {
  
    out <- sce@misc$hexbin[[2]]
    cID <- sce@misc$hexbin[[1]]
  
    if(is.null(out)){
       stop("Compute hexbin representation before plotting.")
    }
  
    if (any(!col %in% colnames(sce@meta.data))) {
        stop("Column cannot be found in slot(sce, 'meta.data').")
    }
  
    name_s <- paste0("sce$", col)
    func <- paste0("x <- ", name_s)
  
    eval(parse(text = func))
  
    x_gene <- GetAssayData(sce, "scale.data")
    
    ind <- match(gene, rownames(x_gene))
  
    if (is.na(ind)) {
        stop("Gene cannot be found.")
    }
  
    x_gene <- as.numeric(x_gene[ind,])
  
  
    .plot_hexbin_fc_helper(x, x_gene, out, cID, col, title, 
          xlab, ylab, colors)
  
})

.plot_hexbin_fc_helper <- function(x, x_gene, out, cID, col, 
    title, xlab, ylab, colors) {
  
    hh <- .make_hexbin_fc_function(x, x_gene, cID)
    out <- as_tibble(out)
  
    if(grepl("^[[:digit:]]", gene )){
        gene <- paste0("G_", gene)
    }
  
    gene <- gsub("-", "_", gene)
  
    col_hh <- paste0(gene, "_", dge)
  
    func1 <- paste0("out$", col_hh, " <- hh")
    eval(parse(text=func1))
  
    .plot_hexbin(out, colour_by=col_hh,
        title=title, xlab=xlab, ylab=ylab)
}

