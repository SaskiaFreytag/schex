#' Group label plot position.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   or \code{\link[Seurat]{Seurat}} object.
#' @param col The name referring to one column in meta data for which
#'   the label position on the plot is calculated for every level. The chosen
#'   column needs to be a factor.
#'
#' @return A tibble.
#' @importFrom tidyr nest
#' @importFrom cluster pam
#' @import Seurat
#' @import SingleCellExperiment
#' @import dplyr
#' @export
#'
#' @examples
#' #' # For Seurat object
#' library(Seurat)
#' data("pbmc_small")
#' pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
#' make_hexbin_label(pbmc_small, col="RNA_snn_res.1")
setGeneric("make_hexbin_label", function(sce, col)
  standardGeneric("make_hexbin_label"))

#' @export
#' @describeIn make_hexbin_label  Group label position for
#'   Seurat object.
setMethod("make_hexbin_label", "Seurat", function(sce, col){

  func <- paste0("!(is.factor(sce$", col, ")|is.character(sce$", col, "))")

  if (any(!col %in% colnames(sce@meta.data))) {
    stop("Column cannot be found in sce@meta.data.")
  }

  if(eval(parse(text=func))){
    stop("The specified column must be a character or a factor.")
  }

  out <- sce@misc$hexbin[[2]]
  cID <- sce@misc$hexbin[[1]]

  if(is.null(out)){
    stop("Compute hexbin representation before plotting.")
  }

  action <- "majority"

  name_s <- paste0("sce$", col)
  func <- paste0("x <- ", name_s)

  eval(parse(text = func))

  hh <- .make_hexbin_function(x, action, cID)
  out <- as_tibble(out)

  label <- paste0(col, "_majority")

  func1 <- paste0("out$", label, " <- hh")
  eval(parse(text=func1))

  label.df_2 <- out %>%
    dplyr::select_("x", "y", label) %>%
    dplyr::group_by_(label) %>%
    nest()

  label.df_3 <- lapply(label.df_2$data, function(x) cluster::pam(x, 1)$medoids)
  names(label.df_3) <- unlist(label.df_2 %>% select_(label))

  label.df_3 <- Reduce(rbind, label.df_3)
  label.df_3 <- data.frame(
    x = label.df_3[, 1],
    y = label.df_3[, 2],
    label = unlist(label.df_2 %>% select_(label))
  )
  rownames(label.df_3) <- NULL

  label.df_3
})

#' @export
#' @describeIn make_hexbin_label  Group label position for
#'   SingleCellExperiment object.
setMethod("make_hexbin_label", "SingleCellExperiment", function(sce, col){

  func <- paste0("!(is.factor(sce$", col, ")|is.character(sce$", col, "))")

  if (any(!col %in% colnames(colData(sce)))) {
    stop("Column cannot be found in colData(sce).")
  }

  if(eval(parse(text=func))){
    stop("The specified column must be a character or a factor.")
  }

  out <- sce@metadata$hexbin[[2]]
  cID <- sce@metadata$hexbin[[1]]

  if(is.null(out)){
    stop("Compute hexbin representation before plotting.")
  }

  action <- "majority"

  name_s <- paste0("sce$", col)
  func <- paste0("x <- ", name_s)

  eval(parse(text = func))

  hh <- .make_hexbin_function(x, action, cID)
  out <- as_tibble(out)

  label <- paste0(col, "_majority")

  func1 <- paste0("out$", label, " <- hh")
  eval(parse(text=func1))

  label.df_2 <- out %>%
    dplyr::select_("x", "y", label) %>%
    dplyr::group_by_(label) %>%
    nest()

  label.df_3 <- lapply(label.df_2$data, function(x) cluster::pam(x, 1)$medoids)
  names(label.df_3) <- unlist(label.df_2 %>% select_(label))

  label.df_3 <- Reduce(rbind, label.df_3)
  label.df_3 <- data.frame(
    x = label.df_3[, 1],
    y = label.df_3[, 2],
    label = unlist(label.df_2 %>% select_(label))
  )
  rownames(label.df_3) <- NULL

  label.df_3
})

