.plot_hexbin <- function(drhex,
                         colour_by = "Cluster_majority",
                         action,
                         colors = NULL,
                         title = NULL,
                         xlab = NULL,
                         ylab = NULL) {
    if (any(!c("x", "y", colour_by) %in% colnames(drhex))) {
        stop("The dataframe must contain columns named 'x', 'y' and label.")
    }

    if (is.null(xlab)) {
        xlab <- "x"
    }

    if (is.null(ylab)) {
        ylab <- "y"
    }

    if (action == "majority") {
        if (is.null(colors)) {
            ggplot(drhex, aes_string("x", "y", fill = colour_by)) +
                geom_hex(stat = "identity") +
                theme_classic() +
                theme(legend.position = "bottom") +
                ggtitle(title) +
                labs(x = xlab, y = ylab) +
                theme(legend.title = element_blank())
        } else {
            ggplot(drhex, aes_string("x", "y", fill = colour_by)) +
                geom_hex(stat = "identity") +
                scale_fill_manual(values = colors) +
                theme_classic() +
                theme(legend.position = "bottom") +
                ggtitle(title) +
                labs(x = xlab, y = ylab) +
                theme(legend.title = element_blank())
        }
    } else {
        ggplot(drhex, aes_string("x", "y", fill = colour_by)) +
            geom_hex(stat = "identity") +
            theme_classic() +
            scale_fill_viridis_c() +
            ggtitle(title) +
            labs(x = xlab, y = ylab) +
            theme(legend.title = element_blank())
    }
}
