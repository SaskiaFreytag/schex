#' @importFrom stats median
#' @importFrom stats sd
.make_hexbin_fc_function <- function(x, x_gene, cID) {
    if (length(levels(x)) != 2) {
        stop("Col needs to be a factor with two levels.")
    }

    x_list <- tapply(x, cID, FUN = function(z) z)
    x_gene_list <- tapply(x_gene, cID, FUN = function(z) z)

    tmp <- lapply(seq(x_list), function(z) {
        data.frame(x_list[[z]], x_gene_list[[z]])
    })

    lapply(tmp, function(z) .fold_change_tapply(z))
}

.fold_change_tapply <- function(tmp) {
    if (any(table(tmp[, 1]) < 1)) {
        return(NULL)
    } else {
        levs <- levels(tmp[, 1])
        x_test <- tmp[which(tmp[, 1] == levs[1]), 2]
        y_test <- tmp[which(tmp[, 1] == levs[2]), 2]

        mean(x_test) - mean(y_test)
    }
}
