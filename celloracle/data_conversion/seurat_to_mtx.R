
library(Seurat)
library(Matrix)
library(colorspace)

args <- commandArgs(trailingOnly = T)
file_path_seurat_object = args[1]
file_path_seurat_object = "PS_seurat_object_190725.Rds"
# load seurat_object
message("loading seurat object ...")

SO <- readRDS(file = file_path_seurat_object)

# check version
if (unlist(SO@version)[1] == 3){
message("  seurat object version is 3x")

# process data
message("processing data ...")
meta <- SO@meta.data
meta$seurat_clusters <- factor(SO@active.ident)

if ("FItSNE" %in% names(SO@reductions)) {
  tsne <- as.data.frame(SO@reductions$FItSNE@cell.embeddings)
  colnames(tsne) <- c("tsne_1", "tsne_2")
  meta <- cbind(meta, tsne)
}


if ("tsne" %in% names(SO@reductions)) {
  tsne <- as.data.frame(SO@reductions$tsne@cell.embeddings)
  colnames(tsne) <- c("tsne_1", "tsne_2")
  meta <- cbind(meta, tsne)
}

if ("umap" %in% names(SO@reductions)) {
  umap <- as.data.frame(SO@reductions$umap@cell.embeddings)
  colnames(umap) <- c("umap_1", "umap_2")
  meta <- cbind(meta, umap)
}



if (length(intersect(c("umap", "tsne", "FItSNE"), names(SO@reductions))) == 0) {
  message("the Seurat object has neither umap data nor tsne data.")
}

# write data
message("making matrix files ...")

active_assay <- SO@active.assay

if (active_assay=="RNA"){
print("active_assay is RNA")
  n <- Matrix::writeMM(obj = SO@assays$RNA@counts, file = "tmp/data.mtx")
  write.csv(colnames(SO@assays$RNA@counts), file="tmp/cells_index.csv")
  write.csv(rownames(SO@assays$RNA@counts), file="tmp/variables_index.csv")

} else if (active_assay=="SCT"){
  print("active_assay is SCT")
  n <- Matrix::writeMM(obj = SO@assays$SCT@data, file = "tmp/data.mtx")
  n <- Matrix::writeMM(obj = SO@assays$SCT@counts, file = "tmp/raw_data.mtx")
  
  write.csv(colnames(SO@assays$SCT@data), file="tmp/cells_index.csv")
  write.csv(rownames(SO@assays$SCT@data), file="tmp/variables_index.csv")

}
write.csv(meta, file = "tmp/meta_data.csv")

# write data type of meta.data
df <- data.frame(row.names = colnames(meta))
for (i in colnames(meta)){
  df[i, "dtype"] = class(meta[,i])
}
write.csv(df, "tmp/meta_data_dtype.csv")

# get and write cluster color information
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n_cluster = length(levels(meta$seurat_clusters))
colors_hex = gg_color_hue(n_cluster)
colors_hex_df = data.frame(colors_hex, row.names = levels(meta$seurat_clusters))
colors_hex_df["cluster"] <- levels(meta$seurat_clusters)
write.csv(colors_hex_df, "tmp/cluster_color_hex.csv")

# variable genes
var_genes <- VariableFeatures(SO)
write.csv(var_genes, file="tmp/var_genes.csv")


}else if(unlist(SO@version)[1] == 2) {


message("  seurat object version is 2x")


# process data
message("processing data ...")
meta <- SO@meta.data
meta$seurat_clusters <- factor(SO@ident)


if ("FItSNE" %in% names(SO@dr)) {
  tsne <- as.data.frame(SO@dr$FItSNE@cell.embeddings)
  colnames(tsne) <- c("tsne_1", "tsne_2")
  meta <- cbind(meta, tsne)
}


if ("tsne" %in% names(SO@dr)) {
  tsne <- as.data.frame(SO@dr$tsne@cell.embeddings)
  colnames(tsne) <- c("tsne_1", "tsne_2")
  meta <- cbind(meta, tsne)
}

if ("umap" %in% names(SO@dr)) {
  umap <- as.data.frame(SO@dr$umap@cell.embeddings)
  colnames(umap) <- c("umap_1", "umap_2")
  meta <- cbind(meta, umap)
}

if (length(intersect(c("umap", "tsne"), names(SO@dr))) == 0) {
  message("the Seurat object has neither umap data nor tsne data.")
}

# write data
message("making matrix files ...")


n <- Matrix::writeMM(obj = SO@data, file = "tmp/data.mtx")
n <- Matrix::writeMM(obj = SO@raw.data, file = "tmp/raw_data.mtx")

write.csv(colnames(SO@data), file="tmp/cells_index.csv")
write.csv(rownames(SO@data), file="tmp/variables_index.csv")


write.csv(meta, file = "tmp/meta_data.csv")

# write data type of meta.data
df <- data.frame(row.names = colnames(meta))
for (i in colnames(meta)){
  df[i, "dtype"] = class(meta[,i])
}
write.csv(df, "tmp/meta_data_dtype.csv")

# get and write cluster color information
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n_cluster = length(levels(meta$seurat_clusters))
colors_hex = gg_color_hue(n_cluster)
colors_hex_df = data.frame(colors_hex, row.names = levels(meta$seurat_clusters))
colors_hex_df["cluster"] <- levels(meta$seurat_clusters)
write.csv(colors_hex_df, "tmp/cluster_color_hex.csv")

# write vargenes
var_genes <- SO@var.genes
write.csv(var_genes, file="tmp/var_genes.csv")
}
