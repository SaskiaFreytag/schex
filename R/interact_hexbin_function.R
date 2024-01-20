#' @importFrom entropy mi.plugin
#' @importFrom stats cor
.interact_hexbin_function <- function(first_x, second_x, interact, cID) {
    if (interact == "corr_spearman") {
        func_if <- !(is.numeric(first_x) | is.numeric(second_x))

        if (func_if) {
            stop("Features need to be numeric.")
        } else {
            res_first <- tapply(first_x, cID, FUN = function(z) z)
            res_second <- tapply(second_x, cID, FUN = function(z) z)

            res <- unlist(lapply(seq_len(length(res_first)), function(x) {
                cor(res_first[[x]], res_second[[x]], method = "spearman")
            }))

            return(res)
        }
    }

    if (interact == "mi") {
        func_if <- !(is.numeric(first_x) | is.numeric(second_x))

        if (func_if) {
            stop("Features need to be numeric.")
        } else {
            res_first <- tapply(first_x, cID, FUN = function(z) z)
            res_second <- tapply(second_x, cID, FUN = function(z) z)

            res <- lapply(seq_len(length(res_first)), function(x) {
                rbind(res_first[[x]], res_second[[x]])
            })

            res <- unlist(lapply(res, function(x) {
                mi.plugin(x)
            }))

            return(res)
        }
    }

    if (interact == "fc") {
        func_if <- !(is.numeric(first_x) | is.numeric(second_x))

        if (func_if) {
            stop("Features need to be numeric.")
        } else {
            res_first <- tapply(first_x, cID, FUN = function(z) z)
            res_second <- tapply(second_x, cID, FUN = function(z) z)

            res <- lapply(seq_len(length(res_first)), function(x) {
                rbind(res_first[[x]], res_second[[x]])
            })

            res <- unlist(lapply(res, function(x) {
                mean(x[1, ] - x[2, ])
            }))

            return(res)
        }
    } else {
        stop("Please specify a valid interaction.")
    }
}
