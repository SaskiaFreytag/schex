#' Plot of feature expression of single cells in bivariate hexagon cells.
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
#' @param action A strings pecifying how meta data of observations in
#'    binned  hexagon cells are to be summarized. Possible actions are
#'    \code{prop_0}, \code{mode}, \code{mean} and \code{median} (see details).
#' @param title A string containing the title of the plot.
#' @param xlab A string containing the title of the x axis.
#' @param ylab A string containing the title of the y axis.
#' @param lower_cutoff For \code{mode}, \code{mean} and \code{median} actions, 
#'     remove expression values below this quantile. Expressed as decimal. 
#'     Default: 0
#' @param upper_cutoff For \code{mode}, \code{mean} and \code{median} actions, 
#'     remove expression values above this quantile. Expressed as decimal.
#'     Default: 1
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
#' @return An object that represents the app. 
#' @import Seurat
#' @import SingleCellExperiment
#' @import ggplot2
#' @importFrom dplyr as_tibble
#' @importFrom methods slotNames
#' @importFrom stats quantile
#' @export
#'
#' @examples
#' # For Seurat object
#' library(Seurat)
#' data("pbmc_small")
#' pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
#' plot_hexbin_feature(pbmc_small, type="counts", feature="TALDO1", 
#'    action="median")
#' plot_hexbin_feature(pbmc_small, type="counts", feature="TALDO1", 
#'    action="median", lower_cutoff=0.2, upper_cutoff=0.5)
#' # For SingleCellExperiment object
#' \dontrun{
#' library(TENxPBMCData)
#' library(scater)
#' tenx_pbmc3k <- TENxPBMCData(dataset = "pbmc3k")
#' rm_ind <- calcAverage(tenx_pbmc3k)<0.1
#' tenx_pbmc3k <- tenx_pbmc3k[!rm_ind,]
#' colData(tenx_pbmc3k) <- cbind(colData(tenx_pbmc3k),
#'    perCellQCMetrics(tenx_pbmc3k))
#' tenx_pbmc3k <- normalize(tenx_pbmc3k)
#' tenx_pbmc3k <- runPCA(tenx_pbmc3k)
#' tenx_pbmc3k <- make_hexbin( tenx_pbmc3k, 20, dimension_reduction = "PCA")
#' plot_hexbin_feature(tenx_pbmc3k, type="logcounts",
#'    feature="ENSG00000135250", action="mean")
#' plot_hexbin_feature(tenx_pbmc3k, type="logcounts",
#'    feature="ENSG00000135250", action="mode")
#' }
#' # For other modalities in Seurat object
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
plot_hexbin_feature <- function(sce, 
    mod="RNA", 
    type,
    feature,
    action,
    title=NULL,
    xlab=NULL,
    ylab=NULL,
    lower_cutoff = 0,
    upper_cutoff = 1) {
  
    out <- .extract_hexbin(sce)
    cID <- .extract_cID(sce)
  
    if(is.null(out)){
       stop("Compute hexbin representation before plotting.")
    }
  
    x <-.prepare_data_feature(sce, mod, type, feature)
  
    .plot_hexbin_feature_helper(x, feature, out, cID, action, title,
                              xlab, ylab, lower_cutoff, upper_cutoff)
  
}

.plot_hexbin_feature_helper <- function(x, feature, out, cID, action, title,
    xlab, ylab, lower_cutoff, upper_cutoff){
  
   if (action %in% c("mean", "median", "mode")) {
        lowend <- quantile(x[x > 0], lower_cutoff)
        highend <- quantile(x[x > 0], upper_cutoff)
        x <- replace(x = x,
            list = x < lowend,
            values = lowend)
        x <- replace(x = x,
            list = x > highend,
            values = highend)
    }
  
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
