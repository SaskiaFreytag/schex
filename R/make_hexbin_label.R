#' Group label plot position.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'    or \code{\link[Seurat]{Seurat-class}} object.
#' @param col The name referring to one column in meta data for which
#'    the label position on the plot is calculated for every level. The chosen
#'    column needs to be a factor.
#'
#' @return A dataframe.
#' @importFrom cluster pam
#' @import Seurat
#' @import SingleCellExperiment
#' @export
#'
#' @examples
#' #' # For Seurat object
#' library(Seurat)
#' data("pbmc_small")
#' pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
#' make_hexbin_label(pbmc_small, col = "RNA_snn_res.1")
setGeneric("make_hexbin_label", function(sce, col) {
  standardGeneric("make_hexbin_label")
})

#' @export
#' @describeIn make_hexbin_label  Group label position for
#'    Seurat object.
setMethod("make_hexbin_label", "Seurat", function(sce, col) {
  func <- paste0("!(is.factor(sce$", col, ")| is.character(sce$", col, "))")

  if (any(!col %in% colnames(sce@meta.data))) {
    stop("Column cannot be found in sce@meta.data.")
  }

  if (eval(parse(text = func))) {
    stop("The specified column must be a character or a factor.")
  }

  out <- sce@misc$hexbin[[2]]
  cID <- sce@misc$hexbin[[1]]

  name_s <- paste0("sce$", col)
  func <- paste0("x <- ", name_s)

  eval(parse(text = func))

  .make_hexbin_label_helper(x, out, cID, col)
})

#' @export
#' @describeIn make_hexbin_label  Group label position for
#'   SingleCellExperiment object.
setMethod("make_hexbin_label", "SingleCellExperiment", function(sce, col) {
  func <- paste0("!(is.factor(sce$", col, ")|is.character(sce$", col, "))")

  if (any(!col %in% colnames(colData(sce)))) {
    stop("Column cannot be found in colData(sce).")
  }

  if (eval(parse(text = func))) {
    stop("The specified column must be a character or a factor.")
  }

  out <- sce@metadata$hexbin[[2]]
  cID <- sce@metadata$hexbin[[1]]

  name_s <- paste0("sce$", col)
  func <- paste0("x <- ", name_s)

  eval(parse(text = func))

  .make_hexbin_label_helper(x, out, cID, col)
})

.make_hexbin_label_helper <- function(x, out, cID, col) {
  if (is.null(out)) {
    stop("Compute hexbin representation before plotting.")
  }

  action <- "majority"

  hh <- .make_hexbin_function(x, action, cID)
  out <- as.data.frame(out)

  label <- paste0(col, "_majority")

  func1 <- paste0("out$", label, " <- hh")
  eval(parse(text = func1))

  label.df_2 <- list()
  for (i in levels(out[, label])) {
    label.df_2[[i]] <- out[which(out[, label] == i), c(1, 2)]
  }

  label.df_3 <- lapply(label.df_2, function(x) cluster::pam(x, 1)$medoids)

  label.df_3 <- Reduce(rbind, label.df_3)
  label.df_3 <- data.frame(
    x = label.df_3[, 1],
    y = label.df_3[, 2],
    label = levels(out[, label])
  )

  label.df_3
}
