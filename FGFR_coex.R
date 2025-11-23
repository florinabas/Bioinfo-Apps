##### Step 0: Load libraries #####
library(Seurat)
library(hdf5r)
library(dplyr)
library(Matrix)
library(ggplot2)
library(pheatmap)

##### Step 1: Read H5AD file #####
path <- "C:/Users/Floarea/Bioinfo app/Bioinfo-Apps/FGFR co-expression/input data"
h5ad_file <- file.path(path, "MCA1.1_adata.h5ad")
adata <- H5File$new(h5ad_file, mode = "r")
h5ls(h5ad_file) # X: 34947 genes x 333778 cells

# Cell names
cell_names <- adata[["obs"]]$read()[["index"]]
# Cell types: substring before first "_"
cell_types <- sub("_.*$", "", cell_names)

# Read gene names and make a dataframe
gene_names <- adata[["var"]]$read()[["index"]]
gene_names_df <- data.frame(gene_name = gene_names, stringsAsFactors = FALSE)

##### Step 2: Identify FGFR genes 
# Detect any occurrence of "fgfr" (case-insensitive) anywhere in the string
fgfr_genes_all <- gene_names[grepl("fgfr", gene_names, ignore.case = TRUE)]
print(fgfr_genes_all)
# output: "Fgfr1"      "Fgfr1op"    "Fgfr1op2"   "Fgfr2"      "Fgfr3"      "Fgfrl1"     "Fgfr4"      "Fgfr3-ps"   "Fgfr3-ps.1"

## focus on the classic gene names also used in the study trying to replicate
fgfr_genes <- c("Fgfr1", "Fgfr2", "Fgfr3", "Fgfr4")  # mouse names
fgfr_rows <- which(gene_names_df$gene_name %in% fgfr_genes)
# Fgfr genes (rows) 4802  4805  4806 23156
fgfr_genes_df <- gene_names_df[fgfr_rows, , drop = FALSE]
message("FGFR genes present in data: ", paste(fgfr_genes_df$gene_name, collapse = ", "))

##### Step 3: Subset X to FGFR genes and FGFR+ cells #####
# Full expression matrix
X <- adata[["X"]]

# Only FGFR rows
fgfr_expr <- X[fgfr_rows, ]


#### check how many tissues left 
# Identify cells where expression > 1 for any FGFR gene
cells_high_expr <- which(colSums(fgfr_expr > 1) >= 1)

# Extract corresponding cell types
cell_types_high_expr <- sub("_.*$", "", cell_names[cells_high_expr])

# Get unique cell types
unique_cell_types <- unique(cell_types_high_expr)

# Count
length(unique_cell_types)
length(unique(cell_types))
# 40 tissues left out of 48


# Identify cells expressing at least one FGFR gene
fgfr_expr_per_cell <- colSums(fgfr_expr > 0)
fgfr_cells <- which(fgfr_expr_per_cell >= 1)
length(fgfr_cells)  # number of FGFR+ cells

# Keep all metadata for these cells
obs_df_fgfr <- adata[["obs"]]$read()[fgfr_cells, , drop = FALSE]
obs_df_fgfr$cell_type <- sub("_.*$", "", obs_df_fgfr$index)  # keep cell type column

# Subset expression matrix to FGFR genes x FGFR+ cells
fgfr_counts_small <- fgfr_expr[, fgfr_cells]
rownames(fgfr_counts_small) <- fgfr_genes_df$gene_name

##### Step 4: Create Seurat object #####
seurat_fgfr <- CreateSeuratObject(
  counts = fgfr_counts_small,
  meta.data = obs_df_fgfr
)

# Normalize and scale
seurat_fgfr <- NormalizeData(seurat_fgfr)
seurat_fgfr <- ScaleData(seurat_fgfr)

##### Step 5: PCA -----
seurat_fgfr <- RunPCA(seurat_fgfr, features = fgfr_genes, npcs = 4)

# Add tiny noise to PCA embeddings to avoid identical points
pca_data <- Embeddings(seurat_fgfr, "pca")
pca_data <- pca_data + matrix(rnorm(length(pca_data), sd = 1e-6), nrow = nrow(pca_data))
seurat_fgfr[["pca"]] <- CreateDimReducObject(embeddings = pca_data, key = "PC_", assay = "RNA")

# Quick check: Fgfr1 vs Fgfr2
FeatureScatter(seurat_fgfr, feature1 = "Fgfr1", feature2 = "Fgfr2", slot="count", raster = FALSE)

##### Step 6: UMAP -----
num_pcs <- ncol(Embeddings(seurat_fgfr, "pca"))
seurat_fgfr <- RunUMAP(seurat_fgfr, dims = 1:num_pcs)

# DimPlot colored by cell type (preserved from metadata)
DimPlot(seurat_fgfr, reduction = "umap", group.by = "cell_type", pt.size = 0.5)

# Extract embeddings for ggplot
umap_df <- as.data.frame(seurat_fgfr@reductions$umap@cell.embeddings)
umap_df$cell_type <- seurat_fgfr$cell_type

# Plot with ggplot
#if(!is.null(dev.list())) dev.off()
#dev.new()
ggplot(umap_df, aes(x = umap_1, y = umap_2, color = cell_type)) +
  geom_point(size = 0.5, alpha = 0.6) +
  theme_minimal() +
  labs(title = "UMAP of FGFR+ cells colored by cell type")

##### Highlight Fgfr1 expression 
Fgfr1_expr <- GetAssayData(seurat_fgfr, layer = "data")["Fgfr1", ]
seurat_fgfr$Fgfr1_positive <- Fgfr1_expr > 0

DimPlot(
  seurat_fgfr,
  reduction = "umap",
  group.by = "Fgfr1_positive",
  cols = c("grey80", "red"),
  pt.size = 0.5
) + ggtitle("Cells expressing Fgfr1")


##### Highlight Fgfr2 expression 
Fgfr2_expr <- GetAssayData(seurat_fgfr, layer = "data")["Fgfr2", ]
seurat_fgfr$Fgfr2_positive <- Fgfr2_expr > 0

DimPlot(
  seurat_fgfr,
  reduction = "umap",
  group.by = "Fgfr2_positive",
  cols = c("grey80", "red"),
  pt.size = 0.5
) + ggtitle("Cells expressing Fgfr2")


##### Highlight Fgfr3 expression 
Fgfr3_expr <- GetAssayData(seurat_fgfr, layer = "data")["Fgfr3", ]
seurat_fgfr$Fgfr3_positive <- Fgfr3_expr > 0

DimPlot(
  seurat_fgfr,
  reduction = "umap",
  group.by = "Fgfr3_positive",
  cols = c("grey80", "red"),
  pt.size = 0.5
) + ggtitle("Cells expressing Fgfr3")


##### Highlight Fgfr4 expression 
Fgfr4_expr <- GetAssayData(seurat_fgfr, layer = "data")["Fgfr4", ]
seurat_fgfr$Fgfr4_positive <- Fgfr4_expr > 0

DimPlot(
  seurat_fgfr,
  reduction = "umap",
  group.by = "Fgfr4_positive",
  cols = c("grey80", "red"),
  pt.size = 0.5
) + ggtitle("Cells expressing Fgfr4")


##### Highlight co-expression of Fgfr1 and Fgfr2 #####
# Identify cells where both Fgfr1 and Fgfr2 are expressed (>0)
seurat_fgfr$Fgfr1_Fgfr2_coexpress <- (Fgfr1_expr > 0) & (Fgfr2_expr > 0)

# Quick check: how many cells co-express both
sum(seurat_fgfr$Fgfr1_Fgfr2_coexpress)
cell_types_12 <- unique(cell_types[which(seurat_fgfr$Fgfr1_Fgfr2_coexpress)])


# Plot UMAP highlighting co-expressing cells
DimPlot(
  seurat_fgfr,
  reduction = "umap",
  group.by = "Fgfr1_Fgfr2_coexpress",
  cols = c("grey80", "red"),  # grey = not co-expressing, red = co-expressing
  pt.size = 0.5
) + ggtitle("Cells co-expressing Fgfr1 and Fgfr2")



##### Highlight co-expression of Fgfr1 and Fgfr3 #####
# Identify cells where both Fgfr1 and Fgfr3 are expressed (>0)
seurat_fgfr$Fgfr1_Fgfr3_coexpress <- (Fgfr1_expr > 0) & (Fgfr3_expr > 0)

# Quick check: how many cells co-express both
sum(seurat_fgfr$Fgfr1_Fgfr3_coexpress)
cell_types_13 <- unique(cell_types[which(seurat_fgfr$Fgfr1_Fgfr3_coexpress)])

# Plot UMAP highlighting co-expressing cells
DimPlot(
  seurat_fgfr,
  reduction = "umap",
  group.by = "Fgfr1_Fgfr3_coexpress",
  cols = c("grey80", "red"),  # grey = not co-expressing, red = co-expressing
  pt.size = 0.5
) + ggtitle("Cells co-expressing Fgfr1 and Fgfr3")



##### Highlight co-expression of Fgfr1 and Fgfr4 #####
# Identify cells where both Fgfr1 and Fgfr4 are expressed (>0)
seurat_fgfr$Fgfr1_Fgfr4_coexpress <- (Fgfr1_expr > 0) & (Fgfr4_expr > 0)

# Quick check: how many cells co-express both
sum(seurat_fgfr$Fgfr1_Fgfr4_coexpress)
cell_types_14 <- unique(cell_types[which(seurat_fgfr$Fgfr1_Fgfr4_coexpress)])

# Plot UMAP highlighting co-expressing cells
DimPlot(
  seurat_fgfr,
  reduction = "umap",
  group.by = "Fgfr1_Fgfr4_coexpress",
  cols = c("grey80", "red"),  # grey = not co-expressing, red = co-expressing
  pt.size = 0.5
) + ggtitle("Cells co-expressing Fgfr1 and Fgfr4")




##### Highlight co-expression of Fgfr2 and Fgfr4 #####
# Identify cells where both Fgfr2 and Fgfr4 are expressed (>0)
seurat_fgfr$Fgfr2_Fgfr4_coexpress <- (Fgfr2_expr > 0) & (Fgfr4_expr > 0)

# Quick check: how many cells co-express both
sum(seurat_fgfr$Fgfr2_Fgfr4_coexpress)
cell_types_24 <- unique(cell_types[which(seurat_fgfr$Fgfr2_Fgfr4_coexpress)])

# Plot UMAP highlighting co-expressing cells
DimPlot(
  seurat_fgfr,
  reduction = "umap",
  group.by = "Fgfr2_Fgfr4_coexpress",
  cols = c("grey80", "red"),  # grey = not co-expressing, red = co-expressing
  pt.size = 0.5
) + ggtitle("Cells co-expressing Fgfr2 and Fgfr4")

##### Highlight co-expression of Fgfr3 and Fgfr4 #####
# Identify cells where both Fgfr3 and Fgfr4 are expressed (>0)
seurat_fgfr$Fgfr3_Fgfr4_coexpress <- (Fgfr3_expr > 0) & (Fgfr4_expr > 0)

# Quick check: how many cells co-express both
sum(seurat_fgfr$Fgfr3_Fgfr4_coexpress)
cell_types_34 <- unique(cell_types[which(seurat_fgfr$Fgfr3_Fgfr4_coexpress)])

# Plot UMAP highlighting co-expressing cells
DimPlot(
  seurat_fgfr,
  reduction = "umap",
  group.by = "Fgfr3_Fgfr4_coexpress",
  cols = c("grey80", "red"),  # grey = not co-expressing, red = co-expressing
  pt.size = 0.5
) + ggtitle("Cells co-expressing Fgfr3 and Fgfr4")



##### Highlight co-expression of Fgfr2 and Fgfr3 #####
# Identify cells where both Fgfr2 and Fgfr3 are expressed (>0)
seurat_fgfr$Fgfr2_Fgfr3_coexpress <- (Fgfr2_expr > 0) & (Fgfr3_expr > 0)

# Quick check: how many cells co-express both
sum(seurat_fgfr$Fgfr2_Fgfr3_coexpress)
cell_types_23 <- unique(cell_types[which(seurat_fgfr$Fgfr2_Fgfr3_coexpress)])

# Plot UMAP highlighting co-expressing cells
DimPlot(
  seurat_fgfr,
  reduction = "umap",
  group.by = "Fgfr2_Fgfr3_coexpress",
  cols = c("grey80", "red"),  # grey = not co-expressing, red = co-expressing
  pt.size = 0.5
) + ggtitle("Cells co-expressing Fgfr2 and Fgfr3")

# Identify cells co-expressing Fgfr1, Fgfr2, and Fgfr3
seurat_fgfr$Fgfr1_2_3_coexpress <- (Fgfr1_expr > 0) & (Fgfr2_expr > 0) & (Fgfr3_expr > 0)

# Count co-expressing cells
sum(seurat_fgfr$Fgfr1_2_3_coexpress)
cell_types_123 <- unique(cell_types[which(seurat_fgfr$Fgfr1_2_3_coexpress)])

# Plot UMAP
DimPlot(
  seurat_fgfr,
  reduction = "umap",
  group.by = "Fgfr1_2_3_coexpress",
  cols = c("grey80", "red"),
  pt.size = 0.5
) + ggtitle("Cells co-expressing Fgfr1, Fgfr2, Fgfr3")





# Identify cells co-expressing Fgfr2, Fgfr3, and Fgfr4
seurat_fgfr$Fgfr2_3_4_coexpress <- (Fgfr2_expr > 0) & (Fgfr3_expr > 0) & (Fgfr4_expr > 0)

# Count co-expressing cells
sum(seurat_fgfr$Fgfr2_3_4_coexpress)
cell_types_234 <- unique(cell_types[which(seurat_fgfr$Fgfr2_3_4_coexpress)])


# Plot UMAP
DimPlot(
  seurat_fgfr,
  reduction = "umap",
  group.by = "Fgfr2_3_4_coexpress",
  cols = c("grey80", "red"),
  pt.size = 0.5
) + ggtitle("Cells co-expressing Fgfr2, Fgfr3, Fgfr4")


# Identify cells co-expressing all four genes
seurat_fgfr$Fgfr1_2_3_4_coexpress <- (Fgfr1_expr > 0) & (Fgfr2_expr > 0) & (Fgfr3_expr > 0) & (Fgfr4_expr > 0)

# Count co-expressing cells
sum(seurat_fgfr$Fgfr1_2_3_4_coexpress)
cell_types_1234 <- unique(cell_types[which(seurat_fgfr$Fgfr1_2_3_4_coexpress)])
original_names_all4 <- cell_names[which(seurat_fgfr$Fgfr1_2_3_4_coexpress)]

# Plot UMAP
DimPlot(
  seurat_fgfr,
  reduction = "umap",
  group.by = "Fgfr1_2_3_4_coexpress",
  cols = c("grey80", "red"),
  pt.size = 0.5
) + ggtitle("Cells co-expressing Fgfr1, Fgfr2, Fgfr3, Fgfr4")

# Get cell names where all 4 FGFR genes are co-expressed


