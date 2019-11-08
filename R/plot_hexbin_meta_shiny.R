#' Plot of meta data of single cell data in bivariate hexagon cells as
#'     shiny instance.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'     or \code{\link[Seurat]{Seurat-class}} object.
#' @param col A string referring to the name of one column in the meta data of
#'    sce by which to colour the hexagons.
#' @param action A string specifying how meta data of observations in
#'    binned  hexagon cells are to be summarized. Possible actions are
#'    \code{majority}, \code{prop_0}, \code{mode}, \code{mean} and
#'    \code{median} (see \code{\link{plot_hexbin_meta}}).
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
#'    observations in the chosen hexagons in a histograms/bar plot below. 
#'
#' @seealso \code{\link{plot_hexbin_meta}} 
#' @return An object that represents the app. 
#' @import Seurat
#' @import shiny
#' @importFrom scales hue_pal
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
#' plot_hexbin_meta_shiny(pbmc_small, col="RNA_snn_res.1", action="majority", 
#'    min_nbins=2, max_nbins=10, dimension_reduction="PCA")
#' }
plot_hexbin_meta_shiny <- function(sce, 
    col, 
    action, 
    min_nbins,
    max_nbins,
    dimension_reduction){
  
  if(action=="prop"){
      stop("prop is not a valid action for shiny instances.")
  }
  
  sce <- make_hexbin(sce, min_nbins, dimension_reduction)
  gg <- plot_hexbin_meta(sce, col, action)
  cID <- .extract_cID(sce)
  gg$data$index <- sort(unique(cID))
  x <- .prepare_data_meta(sce, col)
  
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
          gg <- plot_hexbin_meta(sce, col, action)
          gg
      })
    
    
      output$click_info <- renderPlot({
          sce <- make_hexbin(sce, input$slider, dimension_reduction)
          gg <- plot_hexbin_meta(sce, col, action)
          cID <- .extract_cID(sce)
          gg$data$index <- sort(unique(cID))
          x <-  .prepare_data_meta(sce, col)
          index <- nearPoints(gg$data, input$plot1_click, threshold=10)$index[1]
          index_col <- cID==index
          x_data <- data.frame(groups=x[index_col])
          
          if(action=="majority"){
              col <- scales::hue_pal()(length(levels(x)))
              names(col) <- levels(x)
              ggplot(x_data, aes(x=groups, fill=groups)) + geom_bar() + 
                  theme_classic() + scale_fill_manual(values=col)
          } else {
            qplot(x_data$groups, geom="histogram") + theme_classic()
          }
          
      })
    
    }
  
    shinyApp(ui, server)
}
 
