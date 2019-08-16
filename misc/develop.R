library(schex)
library(Seurat)
library(scater)
library(here)

cbmc.rna <- as.sparse(read.csv(file = "misc/data/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv", sep = ",",
                               header = TRUE, row.names = 1))

cbmc.rna <- CollapseSpeciesExpressionMatrix(cbmc.rna)

cbmc.adt <- as.sparse(read.csv(file = "misc/data/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv", sep = ",",
                               header = TRUE, row.names = 1))

cbmc.adt <- cbmc.adt[setdiff(rownames(x = cbmc.adt), c("CCR5", "CCR7", "CD10")), ]

cbmc <- CreateSeuratObject(counts = cbmc.rna)
cbmc[["ADT"]] <- CreateAssayObject(counts = cbmc.adt)

cbmc <- NormalizeData(cbmc)
cbmc <- FindVariableFeatures(cbmc)
cbmc <- ScaleData(cbmc)
cbmc <- RunPCA(cbmc, verbose = FALSE)
cbmc <- RunTSNE(cbmc, verbose = FALSE)

cbmc <- NormalizeData(cbmc, assay = "ADT", normalization.method = "CLR")
cbmc <- ScaleData(cbmc, assay = "ADT")

cbmc <- make_hexbin(cbmc, nbins=10, dimension_reduction = "TSNE")
summary(cbmc@misc$hexbin$hexbin.matrix[,3])

plot_hexbin_protein(cbmc, assay="ADT", mod="data", protein="CD19",
        action="median", title="Mean of CD19 Protein Expression")

sce <- SingleCellExperiment(assays = list(counts = cbmc.rna))
sce_adt <- SummarizedExperiment(assays=list(counts=cbmc.adt,
                                            data=scale(cbmc.adt)))

altExp(sce, "adt") <- sce_adt

rm_ind <- calculateAverage(sce)<0.1
sce <- sce[!rm_ind,]
colData(sce) <- cbind(colData(sce), perCellQCMetrics(sce))
sce <- normalize(sce)
sce <- runPCA(sce)
sce <- make_hexbin(sce, 20, dimension_reduction = "PCA")


plot_hexbin_protein(sce, mod="adt", type="counts", protein="CD19",
                    action="median", title="Mean of CD19 Protein Expression")
