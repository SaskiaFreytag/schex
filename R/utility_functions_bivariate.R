.plot_hexbin_bivariate_helper_1 <- function(x, feature, out, cID, fan) {
    if (fan) {
        hh <- .make_hexbin_function(x + 1, "mean", cID)
        hh_sd <- .make_hexbin_function(x + 1, "sd", cID)
        hh_cv <- hh_sd / hh
    } else {
        hh <- .make_hexbin_function(x, "mean", cID)
        hh_sd <- .make_hexbin_function(x, "sd", cID)
    }

    out <- as_tibble(out)

    if (grepl("^[[:digit:]]", feature)) {
        feature <- paste0("F_", feature)
    }

    feature <- gsub("-", "_", feature)

    col_hh <- paste0(feature, "_", "mean")
    col_hh_sd <- paste0(feature, "_", "sd")

    func1 <- paste0("out$", col_hh, " <- hh")
    eval(parse(text = func1))

    func2 <- paste0("out$", col_hh_sd, " <- hh_sd")
    eval(parse(text = func2))

    if (fan) {
        col_hh_cv <- paste0(feature, "_", "cv")
        func3 <- paste0("out$", col_hh_cv, " <- hh_cv")
        eval(parse(text = func3))

        .plot_hexbin_bivariate_helper_2(out, x = col_hh, y = col_hh_cv, fan)
    } else {
        .plot_hexbin_bivariate_helper_2(out, x = col_hh, y = col_hh_sd, fan)
    }
}

.plot_hexbin_bivariate_helper_2 <- function(out, x, y, fan) {
    na_ind <- which(is.na(out[, y]))
    if (length(na_ind) > 0) {
        out[na_ind, y] <- 0
    }

    if (fan) {
        out[[y]] <- replace(out[[y]], out[[y]] > 1, 1)
    }

    out$bi_class <- .bi_class(out, x = x, y = y, fan)


    if (length(na_ind) > 0) {
        out$bi_class[na_ind] <- NA
        out[na_ind, y] <- NA
        out[na_ind, x] <- NA
    }

    out$bi_color <- .bivariate_colour_scheme(fan)[
        match(out$bi_class, .bivariate_colour_scheme(fan)[, 1]), 2
    ]
    out$bi_color[is.na(out$bi_class)] <- "grey"

    out
}


.bi_class <- function(out, x, y, fan) {
    if (fan) {
        breaks_x <- as.numeric(cut(out[[x]], 8))
        breaks_y <- as.numeric(cut(out[[y]], breaks = c(-1, 0.25, 0.5, 0.75, 1)))

        paste0(breaks_x, "-", breaks_y)
    } else {
        breaks_x <- as.numeric(cut(out[[x]], 4))
        breaks_y <- as.numeric(cut(out[[y]], 4))

        paste0(breaks_x, "-", breaks_y)
    }
}

.bivariate_colour_scheme <- function(fan) {
    if (fan) {
        matrix(ncol = 2, c(
            "1-1", "#421964",
            "2-1", "#3d4080",
            "3-1", "#326389",
            "4-1", "#30818c",
            "5-1", "#399e8b",
            "5-1", "#57ba7e",
            "6-1", "#92d267",
            "7-1", "#dde15c",
            "8-1", "#6b588e",
            "1-2", "#6b588e",
            "2-2", "#638c9f",
            "3-2", "#638c9f",
            "4-2", "#72b99b",
            "5-2", "#72b99b",
            "6-2", "#c4db7d",
            "7-2", "#c4db7d",
            "8-2", "#8f95b1",
            "1-3", "#8f95b1",
            "2-3", "#8f95b1",
            "3-3", "#8f95b1",
            "4-3", "#acd3a8",
            "5-3", "#acd3a8",
            "6-3", "#acd3a8",
            "7-3", "#acd3a8",
            "8-3", "#b7cac9",
            "1-4", "#b7cac9",
            "2-4", "#b7cac9",
            "3-4", "#b7cac9",
            "4-4", "#b7cac9",
            "5-4", "#b7cac9",
            "6-4", "#b7cac9",
            "7-4", "#b7cac9",
            "8-4", "#b7cac9"
        ), byrow = TRUE)
    } else {
        matrix(ncol = 2, c(
            "1-1", "#402d76",
            "2-1", "#6b588f",
            "3-1", "#9283aa",
            "4-1", "#b8b1c3",
            "1-2", "#30728b",
            "2-2", "#638c9f",
            "3-2", "#8da7b4",
            "4-2", "#b6c3c9",
            "1-3", "#43ad86",
            "2-3", "#72ba9c",
            "3-3", "#98c6b1",
            "4-3", "#bbd2c8",
            "1-4", "#b6da5f",
            "2-4", "#c4db7d",
            "3-4", "#cfdd9e",
            "4-4", "#d8ddbd"
        ), byrow = TRUE)
    }
}
