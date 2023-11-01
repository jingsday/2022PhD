library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

#locations of files
mtx <- "/Users/lidiayung/project/resource/GSE174554_RAW/GSM5319548_SF2979_matrix.mtx.gz"
cells <- "/Users/lidiayung/project/resource/GSE174554_RAW/GSM5319548_SF2979_barcodes.tsv.gz"
features <- "/Users/lidiayung/project/resource/GSE174554_RAW/GSM5319548_SF2979_features.tsv.gz"
str(features)
sf2979 <- ReadMtx(mtx = mtx, cells = cells, features = features,feature.column = 1)


# Initialize the Seurat object with the raw (non-normalized data).
sf2979_p <- CreateSeuratObject(counts = sf2979, project = "sf2979", min.cells = 3, min.features = 200)
str(sf2979_p)

#read metadata and creating tumor_vector to reference
metadata <- read.csv("/Users/lidiayung/project/resource/GSE174554_RAW/GSE174554_Tumor_normal_metadata.txt", header = TRUE,sep=" ")
metadata_subset_sf2979 <- metadata[metadata$Sample == "SF2979", ]
metadata_subset <-metadata_subset_sf2979[metadata_subset_sf2979$Tumor_Normal_annotation == "Tumor", ]
head(metadata_subset)
str(metadata_subset)

tumor_vector<- paste0(metadata_subset$Barcode, "-1")
str(tumor_vector)

#dataset prepared to use
tumor_sf2979 <- subset(sf2979_p, cells = tumor_vector)
tumor_sf2979

tumor_sf2979[["percent.mt"]] <- PercentageFeatureSet(tumor_sf2979, pattern = "^MT-")

VlnPlot(tumor_sf2979, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(tumor_sf2979, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(tumor_sf2979, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# filter >2.5% mitochondrial read counts and <200 expressed genes
filter_tumor_sf2979 <- subset(tumor_sf2979, subset = nFeature_RNA > 200 & percent.mt < 2.5)
str(filter_tumor_sf2979)
VlnPlot(filter_tumor_sf2979, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#normalisation
norm_filter_tumor_sf2979 <- NormalizeData(filter_tumor_sf2979, normalization.method = "LogNormalize", scale.factor = 10000)
str(norm_filter_tumor_sf2979)
#Highly variable features
norm_filter_tumor_sf2979 <- FindVariableFeatures(norm_filter_tumor_sf2979, mean.cutoff = c(0.1, 8),
  dispersion.cutoff = c(1, Inf), selection.method = "mvp")

#1062 HVGs found
str(VariableFeatures(norm_filter_tumor_sf2979))
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(norm_filter_tumor_sf2979), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(norm_filter_tumor_sf2979)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Scaling
all.genes <- rownames(norm_filter_tumor_sf2979)
norm_filter_tumor_sf2979 <- ScaleData(norm_filter_tumor_sf2979, features = all.genes)


#Linear dimensional reduction
norm_filter_tumor_sf2979 <- RunPCA(norm_filter_tumor_sf2979, features = VariableFeatures(object = norm_filter_tumor_sf2979))

# Examine and visualize PCA results a few different ways
print(norm_filter_tumor_sf2979[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(norm_filter_tumor_sf2979, dims = 1:15, reduction = "pca")

DimPlot(norm_filter_tumor_sf2979, reduction="pca")

DimHeatmap(norm_filter_tumor_sf2979, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(norm_filter_tumor_sf2979, dims = 1:15, cells = 500, balanced = TRUE)


norm_filter_tumor_sf2979 <- JackStraw(norm_filter_tumor_sf2979, num.replicate = 100)
norm_filter_tumor_sf2979 <- ScoreJackStraw(norm_filter_tumor_sf2979, dims = 1:20)

JackStrawPlot(norm_filter_tumor_sf2979, dims = 1:15)

ElbowPlot(norm_filter_tumor_sf2979)

#Run non-linear dimensional reduction using the first 15 components 
norm_filter_tumor_sf2979 <- FindNeighbors(norm_filter_tumor_sf2979, dims = 1:15)
norm_filter_tumor_sf2979 <- FindClusters(norm_filter_tumor_sf2979, resolution = 0.3, algorithm = 4)

norm_filter_tumor_sf2979 <- RunUMAP(norm_filter_tumor_sf2979, dims = 1:15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(norm_filter_tumor_sf2979, reduction = "umap")

jpeg("/Users/lidiayung/Downloads/umap_sf2979.png")
DimPlot(norm_filter_tumor_sf2979,reduction = "umap")
dev.off()

#saveRDS(norm_filter_tumor_sf2979, file = "../output/sf2979_tutorial.rds")

#Finding differentially expressed features
cluster2.markers <- FindMarkers(norm_filter_tumor_sf2979, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)


#Find markers per cluster using MAST method used in the paper
norm_filter_tumor_sf2979.markers <- FindAllMarkers(norm_filter_tumor_sf2979, only.pos = FALSE, min.pct = 0.25, test.use="MAST")
#filter using avg_log2FC>1 & p_val_adj<0.05
norm_filter_tumor_sf2979.markers <- norm_filter_tumor_sf2979.markers %>%
    filter(abs(avg_log2FC)>1 & p_val_adj<0.05) 

str(norm_filter_tumor_sf2979.markers)

##check markers in each cluster
library(knitr)
library(kableExtra)
norm_filter_tumor_sf2979.markers  %>% 
  count(cluster)%>%
  kable("html") %>%
  kable_styling()



jpeg("/Users/lidiayung/Downloads/sf2979_all.png")
DoHeatmap(norm_filter_tumor_sf2979, features=norm_filter_tumor_sf2979.markers$gene) +NoLegend()
dev.off()
#ordered
ordered <- norm_filter_tumor_sf2979.markers %>% 
  arrange(desc(avg_log2FC))

jpeg("/Users/lidiayung/Downloads/ordered_sf2979_all.png")
DoHeatmap(norm_filter_tumor_sf2979, features=ordered$gene) +NoLegend()
dev.off()

#top30
norm_filter_tumor_sf2979.markers %>% 
  group_by(cluster)  %>% 
  top_n(n=30, wt=avg_log2FC) -> top30
DoHeatmap(norm_filter_tumor_sf2979, features=top30$gene) +NoLegend()
str(norm_filter_tumor_sf2979.markers)


#top 1/3 
jpeg("/Users/lidiayung/Downloads/sf2979_1three.png")
norm_filter_tumor_sf2979.markers %>%
  group_by(cluster) %>%
  top_frac(n = 1/3, wt = avg_log2FC) -> top_frac30
DoHeatmap(norm_filter_tumor_sf2979, features = top_frac30$gene) + NoLegend()
dev.off()

#Top 50 DEGs per cluster
jpeg("/Users/lidiayung/Downloads/sf2979_top50.png")
norm_filter_tumor_sf2979.markers %>% 
  group_by(cluster)  %>% 
  top_n(n=50, wt=avg_log2FC) -> top50
DoHeatmap(norm_filter_tumor_sf2979, features=top50$gene) +NoLegend()
dev.off()
#Top 100 DEGs per cluster
jpeg("/Users/lidiayung/Downloads/sf2979_top100.png",res=1200)
norm_filter_tumor_sf2979.markers %>% 
  group_by(cluster)  %>% 
  top_n(n=100, wt=avg_log2FC) -> top100
DoHeatmap(norm_filter_tumor_sf2979, features=top100$gene) +NoLegend()
dev.off()
#Top 150 DEGs per cluster
jpeg("/Users/lidiayung/Downloads/sf2979_top150.png",res=1200)
norm_filter_tumor_sf2979.markers %>% 
  group_by(cluster)  %>% 
  top_n(n=150, wt=avg_log2FC) -> top150
DoHeatmap(norm_filter_tumor_sf2979, features=top150$gene) +NoLegend()
dev.off()

#overlapped genes 
jpeg("/Users/lidiayung/Downloads/test2.png")
DoHeatmap(norm_filter_tumor_sf2979, features=selected_genes1third_$gene) +NoLegend()
dev.off()

selected_genes1third_ <- norm_filter_tumor_sf2979.markers %>%
  filter(gene %in% overlapped) %>%
  group_by(cluster) %>%
  top_frac(n = 3/4, wt = avg_log2FC)

selected_genes1third_  %>% 
  count(cluster)


ggplot(norm_filter_tumor_sf2979.markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = ifelse(abs(avg_log2FC) > 1, "Significant", "Not Significant")), alpha = 0.6) +
  #geom_point(aes(color = ifelse(abs(p_val_adj) > 1, "Significant", "Not Significant")), alpha = 0.6) +
  geom_vline(xintercept = 1,)+
  geom_vline(xintercept = -1)+
  #scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  labs(x = "Log2(Fold Change)", y = "-log10(P-value)") +
  theme_classic()


library(svglite)
ggsave("umap_sf2979.svg")
DimPlot(norm_filter_tumor_sf2979, reduction = "umap", width=4, height=4)
dev.off()


#Retrieve union genes
clustered_genes <- norm_filter_tumor_sf2979.markers %>%
   group_by(cluster) %>%
   summarize(genes = list(unique(gene)))

all_genes <- unique(unlist(clustered_genes$genes))
cluster_combinations <- combn(unique(norm_filter_tumor_sf2979.markers$cluster), 2, simplify = FALSE)

intersection_clusters <- purrr::map(cluster_combinations, function(pair) {
  Reduce(intersect, clustered_genes$genes[clustered_genes$cluster %in% pair])
})
intersection_clusters
cat("Union of all genes:", length(all_genes), "genes\n")
write.table(all_genes, file = "/Users/lidiayung/Downloads/SF2979_904all_genes.txt", col.names = FALSE, row.names = FALSE)

primary <- read.table("/Users/lidiayung/Downloads/SF2777_565all_genes.txt", header = FALSE)
str(primary)

overlapped <- intersection <- intersect(primary$V1, all_genes)
write.table(overlapped, file = "/Users/lidiayung/Downloads/SF27772979_overlapped286_genes.txt", col.names = FALSE, row.names = FALSE)


for (i in 1:length(intersection_clusters)) {
  cluster_pair <- cluster_combinations[[i]]
  cat("Intersection of genes between clusters", cluster_pair[1], " and ", cluster_pair[2], ":",
      length(intersection_clusters[[i]]), "genes\n")
}

for (i in 1:length(intersection_clusters)) {
  cluster_pair <- cluster_combinations[[i]]
  intersection_genes <- intersection_clusters[[i]]
  
  cat("Intersection of genes between clusters", cluster_pair[1], " and ", cluster_pair[2], ":",
      length(intersection_genes), "genes: ", paste(intersection_genes, collapse = ", "), "\n")
}
#clusterprofiler
#enrichR
#singleR



