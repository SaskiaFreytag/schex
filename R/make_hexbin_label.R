#' Group label plot position.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @param col The name referring to one column in meta data for which
#'    the label position on the plot is calculated for every level. The chosen
#'    column needs to be a factor.
#'
#' @return A dataframe.
#' @importFrom cluster pam
#' @import SingleCellExperiment
#' @export
#'
#' @examples
#' # For SingleCellExperiment object
#' library(TENxPBMCData)
#' library(scater)
#' tenx_pbmc3k <- TENxPBMCData(dataset = "pbmc3k")
#' rm_ind <- calculateAverage(tenx_pbmc3k) < 0.1
#' tenx_pbmc3k <- tenx_pbmc3k[!rm_ind, ]
#' tenx_pbmc3k <- logNormCounts(tenx_pbmc3k)
#' tenx_pbmc3k <- runPCA(tenx_pbmc3k)
#' tenx_pbmc3k <- make_hexbin(tenx_pbmc3k, 80, dimension_reduction = "PCA")
#' tenx_pbmc3k$random <- factor(sample(1:3, ncol(tenx_pbmc3k), replace = TRUE))
#' make_hexbin_label(tenx_pbmc3k, col = "random")
make_hexbin_label <- function(
    sce,
    col) {
    func <- paste0("!(is.factor(sce$", col, ")| is.character(sce$", col, "))")

    if (eval(parse(text = func))) {
        stop("The specified column must be a character or a factor.")
    }

    out <- .extract_hexbin(sce)
    cID <- .extract_cID(sce)

    x <- .prepare_data_meta(sce, col)

    .make_hexbin_label_helper(x, out, cID, col)
}


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
