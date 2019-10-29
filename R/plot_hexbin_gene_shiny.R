library(shiny)

library(TENxPBMCData)
library(scater)
tenx_pbmc3k <- TENxPBMCData(dataset = "pbmc3k")
rm_ind <- calcAverage(tenx_pbmc3k)<0.1
tenx_pbmc3k <- tenx_pbmc3k[!rm_ind,]
colData(tenx_pbmc3k) <- cbind(colData(tenx_pbmc3k),
    perCellQCMetrics(tenx_pbmc3k))
tenx_pbmc3k <- normalize(tenx_pbmc3k)
tenx_pbmc3k <- runPCA(tenx_pbmc3k)
tenx_pbmc3k <- make_hexbin( tenx_pbmc3k, 10, dimension_reduction = "PCA")
plot_hexbin_gene_shiny(tenx_pbmc3k, type="counts",
   gene="ENSG00000135250", action="mean")
#' plot_hexbin_gene(tenx_pbmc3k, type="logcounts",
#'    gene="ENSG00000135250", action="mode")

setGeneric("plot_hexbin_gene_shiny", function (sce, 
    type,
    action,
    gene) standardGeneric("plot_hexbin_gene_shiny"))

setMethod("plot_hexbin_gene_shiny", "SingleCellExperiment", function (sce,
    type,
    action,
    gene) {
  
    cID <- sce@metadata$hexbin[[1]]
    
    if(is.null(cID)){
      stop("Compute hexbin representation before plotting.")
    }
    
    ind <- match(gene, rownames(sce))
    x <- assays(sce)
    
    if(!type %in% names(x)){
      stop("Specify a valid assay type.")
    }
    
    if (is.na(ind)) {
      stop("Gene cannot be found.")
    }
    
    x <- as.numeric(x[[which(names(x)==type)]][ind,])
    
    .make_shiny_hist_gene(x, sce, type, action, gene, cID)
  
})

setMethod("plot_hexbin_gene_shiny", "Seurat", function (sce,
                                                        type,
                                                        action,
                                                        gene) {
  
  cID <- sce@misc$hexbin$cID
  
  if(is.null(cID)){
    stop("Compute hexbin representation before plotting.")
  }
  
  
  x <- GetAssayData(sce, type)
  
  ind <- match(gene, rownames(x))
  
  if (is.na(ind)) {
    stop("Gene cannot be found.")
  }
  
  x <- as.numeric(x[ind,])
  
  .make_shiny_hist_gene(x, sce, type, action, gene, cID)
  
})

.make_shiny_hist_gene <- function(x, sce, type, action, gene, cID){
  
  gg <- plot_hexbin_gene(sce, type, gene, action)
  gg$data$index <- sort(unique(cID))
  
  ui <- fluidPage(
    fluidRow(
      column(width = 12,
             plotOutput("plot1", height = 400,
                        click = "plot1_click")
      )),
    fluidRow(
      column(width = 2),
      column(width = 8,
             h4("Observations in selected hexagon"),
             plotOutput("click_info", height=150)
      ),
      column(width = 2)
      )
  )
  
  server <- function(input, output) {
    output$plot1 <- renderPlot({
      gg
    })
    
    output$click_info <- renderPlot({
      index <- nearPoints(gg$data, input$plot1_click, threshold=10)$index[1]
      index_col <- cID==index
      x_hist <- x[index_col]
      qplot(x_hist, geom="histogram") + theme_classic()
    })
    
  }
  
  shinyApp(ui, server)
}
  


