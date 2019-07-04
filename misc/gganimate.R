library(gganimate)
library(TENxPBMCData)
library(scater)
library(scran)
library(ggplot2)
library(schex)
library(animation)

tenx_pbmc6k <- TENxPBMCData(dataset = "pbmc6k")
rowData(tenx_pbmc6k)$Mito <- grepl("^MT-", rowData(tenx_pbmc6k)$Symbol_TENx)
tenx_pbmc6k <- calculateQCMetrics(tenx_pbmc6k,
                                  feature_controls = list(Mt = rowData(tenx_pbmc6k)$Mito)
)
tenx_pbmc6k <- tenx_pbmc6k[, !colData(tenx_pbmc6k)$pct_counts_Mt > 50]

libsize_drop <- isOutlier(tenx_pbmc6k$total_counts,
                          nmads = 3,type = "lower", log = TRUE
)
feature_drop <- isOutlier(tenx_pbmc6k$total_features_by_counts,
                          nmads = 3, type = "lower", log = TRUE
)
tenx_pbmc6k <- tenx_pbmc6k[, !(libsize_drop | feature_drop)]
rm_ind <- calcAverage(tenx_pbmc6k)<0
tenx_pbmc6k <- tenx_pbmc6k[!rm_ind,]
tenx_pbmc6k <- normalize(tenx_pbmc6k)
tenx_pbmc6k <- runPCA(tenx_pbmc6k)
tenx_pbmc6k <- runUMAP(tenx_pbmc6k, use_dimred = "PCA", n_neighbors = 50)

tenx_pbmc6k <- make_hexbin(tenx_pbmc6k, nbins=40, dimension_reduction = "UMAP")

gene <- "CD19"
gene_id <- match(gene, rowData(tenx_pbmc6k)$Symbol_TENx)
dat1 <- data.frame(umap1 = reducedDim(tenx_pbmc6k, "UMAP")[,1],
                   umap2 = reducedDim(tenx_pbmc6k, "UMAP")[,2],
                   gene = counts(tenx_pbmc6k[gene_id,]),
                   cID = tenx_pbmc6k@metadata$hexbin$cID)

dat2 <- data.frame(umap1 = tenx_pbmc6k@metadata$hexbin$hexbin.matrix[,1],
                   umap2 = tenx_pbmc6k@metadata$hexbin$hexbin.matrix[,2],
                   gene = tapply(counts(tenx_pbmc6k[gene_id,]),
                  tenx_pbmc6k@metadata$hexbin$cID, function(z) mean(z)),
                  cID=sort(unique(tenx_pbmc6k@metadata$hexbin$cID)))
)

dat3 <- data.frame(cID = tenx_pbmc6k@metadata$hexbin$cID)

dat3$umap1 <- dat2$umap1[match(dat1$cID, dat2$cID)]
dat3$umap2 <- dat2$umap2[match(dat1$cID, dat2$cID)]
dat3$gene <- dat2$gene[match(dat1$cID, dat2$cID)]
dat3$time <- 2
dat1$time <-1

dat3 <- dat3[,colnames(dat1)]

dat <- rbind(dat1, dat3)
dat <- as_tibble(dat)

gene_id_id <-rownames(tenx_pbmc6k)[gene_id]
end_pic <- plot_hexbin_gene(tenx_pbmc6k, type="logcounts", gene=gene_id_id,
                 action="mean", xlab="UMAP1", ylab="UMAP2",
                 title=paste0("Mean of ", gene))

end_pic <- plot(gene)
  ggplot(dat2, aes(x=umap1, y=umap2, colour=gene)) +
  geom_hex() + scale_fill_viridis_c() + theme_classic()

end_pic

base_pic <- dat %>%
  ggplot(aes(
    x = umap1,
    y = umap2,
    colour = gene)) +
  geom_point(alpha=0.75, size=3) +
  scale_colour_viridis_c(rescaler = function(x, to = c(0, 1), from = NULL) {
    ifelse(x<0.6,
      scales::rescale(x, to = to,
       from = c(min(x, na.rm = TRUE), 0.6)),1)}) +
  theme_classic()

base_pic

base_anim <- base_pic + transition_time(time = time)\
anim <- base_anim %>% animate()
anim_save(here::here("misc/", "myanimation.gif"), anim)

base_anim_list <- base_anim %>% animate(device = "png",
  renderer = file_renderer(here::here("misc/gganimate"),
        prefix = "gganim_plot", overwrite = TRUE))

for(i in 1:9){
  ggsave(end_pic, file=here::here("misc/gganimate",
      paste0("gganim_plot010", i, ".png")), width = 4, height = 4)
}

tt <- 0.01010101

dat <- data.frame(frame=101:109,
                 nframes=100,
                 progress=seq(1.01:1.09),
                 frame_time=seq(2+tt, 2.09090909, tt),
                 frame_source=NA)

new_frames <- paste0("/Users/saskiafreytag/Desktop/schex/misc/gganimate/gganim_plot010",
                     1:9, ".png")

base_anim_list1 <- c(base_anim_list, new_frames)
attr(base_anim_list1, "frame_vars") <- rbind(attr(base_anim_list, "frame_vars"),
                                            dat)
anim_save(here::here("misc/", "myanimation11.gif"), base_anim_list, device="png" )
