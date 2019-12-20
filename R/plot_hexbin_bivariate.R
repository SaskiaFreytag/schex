#' Plot of feature expression and uncertainty of single cells in bivariate 
#'    hexagon cells.
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
#' @param fan Logical indicating whether to plot uncertainty surpressing palette.
#' @param title A string containing the title of the plot.
#' @param xlab A string containing the title of the x axis.
#' @param ylab A string containing the title of the y axis.
#'
#' @details This function plots the mean expression and standard deviation of 
#'    any feature in the hexagon cell representation calculated with 
#'    \code{\link{make_hexbin}}. The colour of the hexagon represents a 
#'    combination of both standard deviation and mean expression.
#'
#'    To access the data that has been integrated in the 
#'    \code{\link[Seurat]{Seurat-class}} object specifiy \code{mod="integrated"}.
#'
#' @return A \code{\link{ggplot2}{ggplot}} object.
#' @import Seurat
#' @import SingleCellExperiment
#' @import ggplot2
#' @importFrom dplyr as_tibble
#' @importFrom methods slotNames
#' @importFrom stats quantile
#' @importFrom cowplot ggdraw draw_plot
#' @export
#'
#' @examples
#' # For Seurat object
#' library(Seurat)
#' data("pbmc_small")
#' pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
#' plot_hexbin_bivariate(pbmc_small, type="counts", feature="TALDO1")
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
#' plot_hexbin_bivariate(tenx_pbmc3k, type="logcounts",
#'    feature="ENSG00000135250")
#' }
plot_hexbin_bivariate <- function(sce, 
                                mod="RNA", 
                                type,
                                feature,
                                fan=FALSE,
                                title=NULL,
                                xlab=NULL,
                                ylab=NULL) {
  
  out <- .extract_hexbin(sce)
  cID <- .extract_cID(sce)
  
  if(is.null(out)){
    stop("Compute hexbin representation before plotting.")
  }
  
  x <- .prepare_data_feature(sce, mod, type, feature) 
  
  .plot_hexbin_bivariate_helper(x, feature, out, cID, title,
                              xlab, ylab, fan)
  
}

.plot_hexbin_bivariate_helper <- function(x, feature, out, cID, title,
                                        xlab, ylab, fan){
  
  if(fan){
    hh <- .make_hexbin_function(x+1, "mean", cID)
    hh_sd <- .make_hexbin_function(x+1, "sd", cID)
  } else{
    hh <- .make_hexbin_function(x, "mean", cID)
    hh_sd <- .make_hexbin_function(x, "sd", cID)
  }
  
  out <- as_tibble(out)
  
  if(grepl("^[[:digit:]]", feature )){
    feature <- paste0("F_", feature)
  }
  
  feature <- gsub("-", "_", feature)
  
  col_hh <- paste0(feature, "_", "mean")
  col_hh_sd <- paste0(feature, "_", "sd")
  
  func1 <- paste0("out$", col_hh, " <- hh")
  eval(parse(text=func1))
  
  func2 <- paste0("out$", col_hh_sd, " <- hh_sd")
  eval(parse(text=func2))
  
  .plot_bivariate(out, title=title, xlab=xlab, ylab=ylab, x=col_hh,
              y=col_hh_sd, fan)
}


.bivariate_colour_scheme <- function(fan) {
  
  if(fan){
    matrix(ncol=2, c(
    "1-1" ,"#421964", 
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
    "8-4", "#b7cac9")
  } else {
    matrix(ncol=2, c(
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
    "4-4", "#d8ddbd"), byrow=TRUE)
  }
}
  

bi_class <- function(out, x, y, fan){
  
  if(fan){
    
    breaks_x <- as.numeric(cut(out[[x]],8))
    breaks_y <- as.numeric(cut(out[[y]],4))
    
    paste0(breaks_x, "-", breaks_y)
    
  } else {
    
    breaks_x <- as.numeric(cut(out[[x]],4))
    breaks_y <- as.numeric(cut(out[[y]],4))
    
    paste0(breaks_x, "-", breaks_y)
    
  }
  
}

.plot_bivariate <- function(out, title=title, xlab=xlab, ylab=ylab, x, y, fan){
  
  na_ind <- which(is.na(out[, y]))
  if(length(na_ind)>0){
    out[na_ind, y] <- 0
  }
  out$bi_class <- bi_class(out, x = x, y = y, fan)
  if(length(na_ind)>0){
    out$bi_class[na_ind] <- NA
  }
  
  out$bi_color <- .bivariate_colour_scheme[
    match(out$bi_class, .bivariate_colour_scheme[,1]), 2]
  out$bi_color[is.na(out$bi_color)] <- "grey"
  
  if (is.null(title)) {
    title <- paste("Bivariate plot", strsplit(x, "_")[[1]][1])
  }
  
  if (is.null(xlab)) {
    xlab <- "x"
  }
  
  if (is.null(ylab)) {
    ylab <- "y"
  }
  
  g1 <- ggplot(out, aes_string("x", "y")) +
    geom_hex(stat = "identity",  fill = out$bi_color) + 
    theme_classic()  + ggtitle(title) +
    labs(x = xlab, y = ylab) + theme(legend.title = element_blank())
  
  g2 <- .legend_builder(out, colours=.bivariate_colour_scheme[,2], x, y, fan)
  
  ggdraw() +
    draw_plot(g1, x = 0, y = 0, width = 0.75, height = 1) +
    draw_plot(g2, x = 0.75, y = 0.65, width = 0.225, height = 0.225)
  
}


.legend_builder <- function(out, colours=.bivariate_colour_scheme[,2], x,
                           y, fan) {
  
  if(fan){
    
    
    
  } else {
  
    dat_leg <- data.frame(y=rep(1:4, 4), x=rep(1:4,each=4), col=colours)
    breaks_x <- seq(range(out[[x]])[1],range(out[[x]])[2], length.out=5)
    new_x <- vapply(seq_len(length(breaks_x)-1), function(xx) 
      (breaks_x[xx]+breaks_x[xx+1])/2,double(1))
    breaks_y <- seq(range(out[[y]])[1],range(out[[y]])[2], length.out=5)
    new_y <- vapply(seq_len(length(breaks_y)-1), function(xx) 
      (breaks_y[xx]+breaks_y[xx+1])/2,double(1))
    dat_leg$y <- new_y[dat_leg$y]
    dat_leg$x <- new_x[dat_leg$x]
  
    ggplot(dat_leg, aes(x=x, y=y)) + geom_tile(stat = "identity", 
        fill=dat_leg$col) + theme_classic() + 
        scale_x_continuous(limits = c(breaks_x[1], breaks_x[5]), 
          expand = c(0, 0), breaks=round(breaks_x,3)) +
      scale_y_continuous(limits = c(breaks_y[1], breaks_y[5]), expand = c(0, 0),
          breaks=round(breaks_y,3))  + xlab(x) + ylab(y) +
      theme(plot.margin = margin(0, 0, 0, 0, "cm"), 
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  }
}

