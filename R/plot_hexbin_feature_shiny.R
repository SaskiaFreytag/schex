#' Plot of feature expression of single cell data in bivariate hexagon cells as
#'     shiny instance.
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
#'    \code{prop_0}, \code{mode}, \code{mean} and \code{median} (see 
#'    \code{\link{plot_hexbin_feature}}).
#' @param min_nbins The miniumum number of bins partitioning the range of the 
#'    first component of the chosen dimension reduction.
#' @param max_nbins The miniumum number of bins partitioning the range of the 
#'    first component of the chosen dimension reduction. 
#' @param dimension_reduction A string indicating the reduced dimension
#'    result to calculate hexagon cell representation of.
#'
#' @details This function opens a shiny instance, which allows to investigate 
#'    the effect of the resolution parameter. The user can change the resolution
#'    using the slider. Each hexagon is clickable, which will plot the 
#'    observations in the chosen hexagons in a histograms below. 
#'
#' @seealso \code{\link{plot_hexbin_feature}}
#' @return An object that represents the app. 
#' @import Seurat
#' @import shiny
#' @import SingleCellExperiment
#' @import ggplot2
#' @importFrom dplyr as_tibble
#' @importFrom methods slotNames
#' @export
#'
#' @examples
#' # For Seurat object
#' \dontrun{
#' library(Seurat)
#' data("pbmc_small")
#' plot_hexbin_feature_shiny(pbmc_small, type="counts", feature="TALDO1", 
#'    action="prop_0", min_nbins=2, max_nbins=10, dimension_reduction="PCA",
#'    mod="RNA")
#' }
plot_hexbin_feature_shiny <- function(sce,
    mod="RNA",
    type, 
    feature, 
    action, 
    min_nbins,
    max_nbins,
    dimension_reduction){
  
    sce <- make_hexbin(sce, min_nbins, dimension_reduction)
    gg <- plot_hexbin_feature(sce, mod, type, feature, action)
    cID <- .extract_cID(sce)
    gg$data$index <- sort(unique(cID))
    x <- .prepare_data_feature(sce, mod, type, feature)
  
    ui <- fluidPage(
        fluidRow(
            column(width = 12,
             plotOutput("plot1", height = 400,
                        click = "plot1_click")
        )),
        fluidRow(
            column(width = 4,
                sliderInput("slider", NULL, min_nbins, max=max_nbins, 
                value=min_nbins)),
            column(width = 6,
                h4("Observations in selected hexagon"),
                plotOutput("click_info", height=150)
        ),
        column(width = 2)
        )
    )

  
    server <- function(input, output) {
    
        output$plot1 <- renderPlot({
            sce <- make_hexbin(sce, input$slider, dimension_reduction)
            gg <- plot_hexbin_feature(sce, mod, type, feature, action)
            gg
        })
    
    
      output$click_info <- renderPlot({
          sce <- make_hexbin(sce, input$slider, dimension_reduction)
          gg <- plot_hexbin_feature(sce, mod, type, feature, action)
          cID <- .extract_cID(sce)
          gg$data$index <- sort(unique(cID))
          x <- .prepare_data_feature(sce, mod, type, feature)
          index <- nearPoints(gg$data, input$plot1_click, threshold=10)$index[1]
          index_col <- cID==index
          x_hist <- x[index_col]
          qplot(x_hist, geom="histogram") + theme_classic()
      })
    
    }
  
    shinyApp(ui, server)
}
 
