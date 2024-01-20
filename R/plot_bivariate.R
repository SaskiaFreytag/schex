#' @importFrom dplyr tibble
.plot_bivariate <- function(
    out, title = title, xlab = xlab, ylab = ylab, fan,
    feature) {
    if (is.null(title)) {
        title <- paste("Bivariate plot", feature)
    }

    if (is.null(xlab)) {
        xlab <- "x"
    }

    if (is.null(ylab)) {
        ylab <- "y"
    }


    xx <- paste0(feature, "_mean")

    if (fan) {
        yy <- paste0(feature, "_cv")

        g1 <- ggplot(out, aes_string("x", "y")) +
            geom_hex(stat = "identity", fill = out$bi_color) +
            theme_classic() +
            ggtitle(title) +
            labs(x = xlab, y = ylab) +
            theme(legend.title = element_blank())

        dat <- .colourfan_grob(nrow = 4, ncol = 8)
        id <- rep(seq(1, (4 * 8)), each = 22)

        label.x.pos <- dat[dat$x == min(dat$x), ]
        label.y.pos <- dat[dat$y == min(dat$y), ]

        label.x.pos <- .transform_radial(tibble(x = seq(0, 1, 1 / 4), y = 1),
            yoff = 0.08
        )
        label.y.pos <- .transform_radial(tibble(x = 0, y = seq(0, 1, 1 / 4)),
            xoff = -0.3, yoff = -0.2
        )

        label.x <- round(seq(1, max(out[[xx]], na.rm = TRUE), length.out = 5), 2)

        colours_fun <- .bivariate_colour_scheme(fan)[, 2]

        grid.newpage()
        pushViewport(viewport(
            width = 0.65, height = 1,
            x = 0, y = 0, just = c("left", "bottom")
        ))
        grid.draw(ggplotGrob(g1))
        popViewport()
        pushViewport(viewport(
            width = 0.25, height = 0.275,
            x = 0.725, y = 0.65, just = c("left", "bottom")
        ))
        grid.draw(polygonGrob(dat$x, dat$y, id, gp = gpar(
            fill = colours_fun,
            col = colours_fun, lwd = 1, lty = 1
        )))

        grid.draw(textGrob(label.x, label.x.pos$x, label.x.pos$y,
            gp = gpar(cex = 0.75)
        ))
        grid.draw(textGrob(c(">100%", "75%", "50%", "25%", "0%"),
            label.y.pos$x, label.y.pos$y,
            gp = gpar(cex = 0.75)
        ))
        popViewport()

        pushViewport(viewport(
            x = 0.78, y = 0.75,
            width = 0.2, height = 0.02
        ))
        grid.draw(textGrob("uncertainty", 0.1, 0.1,
            rot = 300,
            gp = gpar(cex = 0.75)
        ))
        popViewport()

        pushViewport(viewport(
            x = 0.855, y = 0.6,
            width = 0.2, height = 0.02
        ))
        grid.draw(textGrob("expression", 0.5, 0.2, gp = gpar(cex = 0.75)))
        popViewport()
    } else {
        yy <- paste0(feature, "_sd")

        dat_leg <- data.frame(
            y = rep(seq(1, 4), 4), x = rep(seq(1, 4), each = 4),
            col = .bivariate_colour_scheme(fan)[, 2]
        )

        breaks_x <- seq(0, max(out[[xx]], na.rm = TRUE), length.out = 5)
        new_x <- vapply(seq_len(length(breaks_x) - 1), function(xxx) {
            (breaks_x[xxx] + breaks_x[xxx + 1]) / 2
        }, double(1))

        breaks_y <- seq(0, max(out[[yy]], na.rm = TRUE), length.out = 5)
        new_y <- vapply(seq_len(length(breaks_y) - 1), function(xxx) {
            (breaks_y[xxx] + breaks_y[xxx + 1]) / 2
        }, double(1))

        dat_leg$y <- new_y[dat_leg$y]
        dat_leg$x <- new_x[dat_leg$x]

        gl <- ggplot(dat_leg, aes(x = x, y = y)) +
            geom_tile(
                stat = "identity",
                fill = dat_leg$col
            ) +
            theme_classic() +
            scale_x_continuous(
                limits = c(breaks_x[1], breaks_x[5]),
                expand = c(0, 0), breaks = round(breaks_x, 3)
            ) +
            scale_y_continuous(
                limits = c(breaks_y[1], breaks_y[5]),
                expand = c(0, 0), breaks = round(breaks_y, 3)
            ) +
            theme(
                plot.margin = margin(0, 0, 0, 0, "cm"),
                axis.text.x = element_text(angle = 45, hjust = 1),
                axis.title = element_blank()
            )

        g1 <- ggplot(out, aes_string("x", "y")) +
            geom_hex(stat = "identity", fill = out$bi_color) +
            theme_classic() +
            ggtitle(title) +
            labs(x = xlab, y = ylab) +
            theme(legend.title = element_blank())

        grid.newpage()

        pushViewport(viewport(
            x = 0, y = 0, width = 0.75, height = 1,
            just = c("left", "bottom")
        ))
        grid.draw(ggplotGrob(g1))
        popViewport()

        pushViewport(viewport(x = 0.875, y = 0.8, width = 0.2, height = 0.225))
        grid.draw(ggplotGrob(gl))
        popViewport()

        pushViewport(viewport(
            x = 0.9, y = 0.685,
            width = 0.2, height = 0.02
        ))
        grid.draw(textGrob("mean expression", 0.5, 0.2, gp = gpar(cex = 0.75)))
        popViewport()

        pushViewport(viewport(
            x = 0.84, y = 0.8,
            width = 0.2, height = 0.02
        ))
        grid.draw(textGrob("standard deviation", 0.1, 0.1,
            rot = 90,
            gp = gpar(cex = 0.75)
        ))
        popViewport()
    }
}


.transform_radial <- function(data, xoff = 0, yoff = 0) {
    phi <- (data$x * 60 - 30) * (pi / 180)
    Y <- (data$y + yoff) * cos(phi) - xoff * sin(60 * pi / 360)
    X <- (data$y + yoff) * sin(phi) + 0.5 + xoff * cos(60 * pi / 360)
    tibble(x = X, y = Y)
}

.colourfan_grob <- function(nrow, ncol, nmunch = 10) {
    dx <- 1 / ncol
    dy <- 1 / nrow

    x <- rep((seq(0, (ncol - 1))) / ncol, nrow)
    y <- rep((seq((nrow - 1), 0)) / nrow, each = ncol)

    x <- unlist(lapply(x, function(x) {
        c(
            x + dx * (seq(0, nmunch)) / nmunch,
            x + dx * (seq(nmunch, 0)) / nmunch
        )
    }))
    y <- unlist(lapply(y, function(y) {
        c(
            rep(y, nmunch + 1),
            rep(y + dy, nmunch + 1)
        )
    }))

    data <- .transform_radial(tibble(x, y))
}
