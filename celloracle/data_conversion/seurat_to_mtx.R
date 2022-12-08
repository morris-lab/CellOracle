
library(Seurat)
library(Matrix)
library(colorspace)

args <- commandArgs(trailingOnly = T)
file_path_seurat_object = args[1]

# Extract cluster color information
get_cluster_colors <- function(clusters){

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  n_cluster = length(levels(clusters))
  colors_hex = gg_color_hue(n_cluster)
  colors_hex_df = data.frame(colors_hex, row.names = levels(clusters))
  colors_hex_df["cluster"] <- levels(clusters)
  return(colors_hex_df)
}

export_SeuratObjectV3 <- function(SO){

  # Make directory to save files
  dir.create("tmp")

  # Save gene expression matrix
  for (assay_name in names(SO@assays)){
    tryCatch({
    message(paste0("Processing an assay: ", assay_name))

    dim_ <- dim(SO@assays[assay_name][[1]]@data)
    if ((dim_[1] > 0) & (dim_[2] > 0)){
      n <- Matrix::writeMM(obj = t(SO@assays[assay_name][[1]]@data),
                           file = paste0("tmp/assay_", assay_name,  "_data.mtx"))
      write.csv(colnames(SO@assays[assay_name][[1]]@data), file = paste0("tmp/assay_", assay_name,  "_cells.csv"))
      write.csv(rownames(SO@assays[assay_name][[1]]@data), file = paste0("tmp/assay_", assay_name,  "_genes.csv"))

    }
    dim_ <- dim(SO@assays[assay_name][[1]]@counts)
    #print(dim_)
    if ((dim_[1] > 0) & (dim_[2] > 0)){
      n <- Matrix::writeMM(obj = t(SO@assays[assay_name][[1]]@counts),
                       file = paste0("tmp/assay_", assay_name,  "_rawdata.mtx"))
      #write.csv(colnames(SO@assays[assay_name][[1]]@counts), file = paste0("tmp/assay_", assay_name,  "_cells.csv"))
      #write.csv(rownames(SO@assays[assay_name][[1]]@counts), file = paste0("tmp/assay_", assay_name,  "_genes.csv"))
    }
    }, error = function(e) {
    message(paste0("Error. The following assay was skipped: ", assay_name))

    }
  )
  }

  for (reduction_name in names(SO@reductions)){
    write.csv(as.data.frame(SO@reductions[reduction_name][[1]]@cell.embeddings),
              file = paste0("tmp/reduction_", reduction_name,  ".csv"))
  }

  # Meta data
  meta <- SO@meta.data
  meta$active_ident <- factor(SO@active.ident)
  write.csv(meta, file = "tmp/meta_data.csv")

   # write data type of meta.data
  meta_data_dtype <- data.frame(row.names = colnames(meta))
  for (i in colnames(meta)){
    if (length(class(meta[,i])) == 1){
      class_ = class(meta[,i])
    } else {
      n_class <- length(class(meta[,i]))
      class_ <- class(meta[,i])[n_class]
    }
          meta_data_dtype[i, "dtype"] <- class_
  }
  write.csv(meta_data_dtype, "tmp/meta_data_dtype.csv")

  # Color info
  clusters <- factor(SO@active.ident)
  colors_hex_df <- get_cluster_colors(clusters = clusters)
  write.csv(colors_hex_df, "tmp/cluster_color_hex.csv")

  # Variable gene list
  var_genes <- VariableFeatures(SO)
  write.csv(var_genes, file="tmp/var_genes.csv")

}


export_SeuratObjectV4 <- export_SeuratObjectV3

export_SeuratObjectV2 <- function(SO){

  # Make directory to save files
  dir.create("tmp")

  # Save gene expression matrix
  n <- Matrix::writeMM(obj = t(SO@data),
                       file = paste0("tmp/assay_RNA_data.mtx"))
  n <- Matrix::writeMM(obj = t(SO@raw.data[rownames(SO@data), colnames(SO@data)]),
                       file = paste0("tmp/assay_RNA_rawdata.mtx"))
  write.csv(colnames(SO@data), file = paste0("tmp/assay_RNA_cells.csv"))
  write.csv(rownames(SO@data), file = paste0("tmp/assay_RNA_genes.csv"))


  for (reduction_name in names(SO@dr)){
    write.csv(as.data.frame(SO@dr[reduction_name][[1]]@cell.embeddings),
              file = paste0("tmp/reduction_", reduction_name,  ".csv"))
  }

  # Meta data
  meta <- SO@meta.data
  meta$active_ident <- factor(SO@ident)
  write.csv(meta, file = "tmp/meta_data.csv")

  # write data type of meta.data
  meta_data_dtype <- data.frame(row.names = colnames(meta))
  for (i in colnames(meta)){
    meta_data_dtype[i, "dtype"] = class(meta[,i])
  }
  write.csv(meta_data_dtype, "tmp/meta_data_dtype.csv")

  # Color info
  clusters <- factor(SO@ident)
  colors_hex_df <- get_cluster_colors(clusters = clusters)
  write.csv(colors_hex_df, "tmp/cluster_color_hex.csv")

  # Variable gene list
  var_genes <- SO@var.genes
  write.csv(var_genes, file="tmp/var_genes.csv")

}


main <- function(){
  # load seurat_object
  message("loading seurat object ...")
  suppressMessages(SO <- readRDS(file = file_path_seurat_object))

  # check version
  if (unlist(SO@version)[1] == 4){
    message("  seurat object version is 4x")
    export_SeuratObjectV4(SO)

  }else if (unlist(SO@version)[1] == 3){
    message("  seurat object version is 3x")
    export_SeuratObjectV3(SO)

  }else if(unlist(SO@version)[1] == 2) {

    message("  seurat object version is 2x")
    export_SeuratObjectV2(SO)

  }
}


main()
