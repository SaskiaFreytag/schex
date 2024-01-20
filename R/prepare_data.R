setGeneric(".prepare_data_feature", function(
    sce,
    mod = "RNA",
    type,
    feature) {
    standardGeneric(".prepare_data_feature")
})


setMethod(".prepare_data_feature", "SingleCellExperiment", function(
    sce,
    mod = "RNA",
    type,
    feature) {
    if (mod == "RNA") {
        ind <- match(feature, rownames(sce))

        if (is.na(ind)) {
            stop("Gene cannot be found.")
        }

        x <- assays(sce)

        if (!type %in% names(x)) {
            stop("Specify a valid assay type.")
        }

        x <- as.numeric(x[[which(names(x) == type)]][ind, ])
    } else {
        if (!mod %in% altExpNames(sce)) {
            stop("Specify a valid modularity.")
        }

        if (!type %in% assayNames(altExp(sce))) {
            stop("Specify a valid assay type.")
        }

        x <- assays(altExp(sce, mod))
        x <- x[[which(names(x) == type)]]

        ind <- match(feature, rownames(x))

        if (is.na(ind)) {
            stop("Feature cannot be found.")
        }

        x <- as.numeric(x[ind, ])
    }

    return(x)
})

setGeneric(".prepare_data_meta", function(
    sce,
    col) {
    standardGeneric(".prepare_data_meta")
})


setMethod(".prepare_data_meta", "SingleCellExperiment", function(
    sce,
    col) {
    if (any(!col %in% colnames(colData(sce)))) {
        stop("Column cannot be found in colData(sce).")
    }

    name_s <- paste0("sce$", col)
    func <- paste0("x <- ", name_s)

    eval(parse(text = func))

    return(x)
})
