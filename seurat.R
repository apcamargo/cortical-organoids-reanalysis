library("Seurat")
library("tximport")

file1 <- file.path("Quantification", "1_month", "alevin", "quants_mat.gz")
file2 <- file.path("Quantification", "3_month", "alevin", "quants_mat.gz")
file3 <- file.path("Quantification", "6_month", "alevin", "quants_mat.gz")
file4 <- file.path("Quantification", "10_month", "alevin", "quants_mat.gz")

txi1 <- tximport(file1, type="alevin")
so1 <- CreateSeuratObject(counts = round(txi1$counts), min.cells = 3, min.features = 200, project = "1_month")

txi2 <- tximport(file2, type="alevin")
so2 <- CreateSeuratObject(counts = round(txi2$counts), min.cells = 3, min.features = 200, project = "3_month")

txi3 <- tximport(file3, type="alevin")
so3 <- CreateSeuratObject(counts = round(txi3$counts), min.cells = 3, min.features = 200, project = "6_month")

txi4 <- tximport(file4, type="alevin")
so4 <- CreateSeuratObject(counts = round(txi4$counts), min.cells = 3, min.features = 200, project = "10_month")

so1 <- PercentageFeatureSet(so1, pattern = "^MT-", col.name = "percent.mt")
so2 <- PercentageFeatureSet(so2, pattern = "^MT-", col.name = "percent.mt")
so3 <- PercentageFeatureSet(so3, pattern = "^MT-", col.name = "percent.mt")
so4 <- PercentageFeatureSet(so4, pattern = "^MT-", col.name = "percent.mt")

so1 <- subset(so1, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 5)
so2 <- subset(so2, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 5)
so3 <- subset(so3, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 5)
so4 <- subset(so4, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 5)

so1 <- NormalizeData(so1)
so2 <- NormalizeData(so2)
so3 <- NormalizeData(so3)
so4 <- NormalizeData(so4)

so1 <- FindVariableFeatures(so1, selection.method = "vst", nfeatures = 2500)
so2 <- FindVariableFeatures(so2, selection.method = "vst", nfeatures = 2500)
so3 <- FindVariableFeatures(so3, selection.method = "vst", nfeatures = 2500)
so4 <- FindVariableFeatures(so4, selection.method = "vst", nfeatures = 2500)

object_list <- c(so1, so2, so3, so4)
anchors <- FindIntegrationAnchors(object.list = object_list)

so_combined <- IntegrateData(anchorset = anchors, dims = 1:30)

DefaultAssay(so_combined) <- "integrated"
so_combined <- ScaleData(so_combined)

so_combined <- FindNeighbors(so_combined, dims = 1:10)
so_combined <- FindClusters(so_combined, resolution = 0.5)
markers <- FindAllMarkers(so_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

so_combined <- RunPCA(so_combined)
so_combined <- RunUMAP(so_combined, dims = 1:30, umap.method = "uwot")
so_combined <- RunTSNE(so_combined)
DefaultAssay(so_combined) <- "RNA"

ExportToCellbrowser(so_combined, "cellbrowser", dataset.name="cortical.organoids", reductions=c("tsne", "umap"), cluster.field="Batch")
