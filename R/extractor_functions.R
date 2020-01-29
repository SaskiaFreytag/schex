setGeneric(".extract_cID", function (sce) standardGeneric(".extract_cID"))

setMethod(".extract_cID", "SingleCellExperiment", function (sce) {

    sce@metadata$hexbin[[1]]
  
}) 

setMethod(".extract_cID", "Seurat", function (sce) {

    sce@misc$hexbin[[1]]
  
}) 

setGeneric(".extract_hexbin", function (sce) standardGeneric(".extract_hexbin"))

setMethod(".extract_hexbin", "SingleCellExperiment", function (sce) {

    sce@metadata$hexbin[[2]]
  
}) 

setMethod(".extract_hexbin", "Seurat", function (sce) {
  
    sce@misc$hexbin[[2]]
  
}) 