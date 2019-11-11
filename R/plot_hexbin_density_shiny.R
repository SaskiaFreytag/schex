#' Plot of density of observations from single cell data
#'    in bivariate hexagon cells as shiny instance.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'    or \code{\link[Seurat]{Seurat-class}} object.
#' @param min_nbins The miniumum number of bins partitioning the range of the 
#'    first component of the chosen dimension reduction.
#' @param max_nbins The miniumum number of bins partitioning the range of the 
#'    first component of the chosen dimension reduction. 
#' @param dimension_reduction A string indicating the reduced dimension
#'    result to calculate hexagon cell representation of.
#'    
#' @details This function opens a shiny instance, which allows to investigate 
#'    the effect of the resolution parameter. The user can change the resolution
#'    using the slider. 
#'
#' @seealso \code{\link{plot_hexbin_density}}
#' @return An object that represents the app. 
#' @import Seurat
#' @import shiny
#' @import SingleCellExperiment
#' @import ggplot2
#' @importFrom dplyr as_tibble
#' @export
#'
#' @examples
#' # For Seurat object
#' \dontrun{
#' library(Seurat)
#' data("pbmc_small")
#' plot_hexbin_density_shiny(pbmc_small,3, 10, dimension_reduction = "PCA")
#' }
plot_hexbin_density_shiny <- function(sce,
                                      min_nbins,
                                      max_nbins,
                                      dimension_reduction){
  
  sce <- make_hexbin(sce, min_nbins, dimension_reduction)
  out <- .extract_hexbin(sce)
  gg <- .plot_hexbin_density_helper(out, NULL, NULL, NULL)
  
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
      column(width = 2)
      )
    )
  
  
  server <- function(input, output) {
    
    output$plot1 <- renderPlot({
      sce <- make_hexbin(sce, input$slider, dimension_reduction)
      out <- .extract_hexbin(sce)
      gg <- .plot_hexbin_density_helper(out, NULL, NULL, NULL)
      gg
    })
    
    
  }
  
  shinyApp(ui, server)
}





