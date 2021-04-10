# Analysis of digital expression data by Seurat
https://satijalab.org/seurat/
All analysis was conducted on Rstudio

## Reading of libraries
```{r}
library(needs)
needs(dplyr,Seurat,ggplot2,SeuratWrappers,SeuratData,patchwork,broom,sinaplot,monocle3,glmGamPoi)
```

## Import data
```{r}
data1=read.table(file = "shrimp1.txt.gz", header = TRUE, row.names = 1)
data2=read.table(file = "shrimp2.txt.gz", header = TRUE, row.names = 1)
data3=read.table(file = "shrimp3.txt.gz", header = TRUE, row.names = 1)
```

## Initialize the Seurat object with the raw (non-normalized data)
```{r}
Mj1=CreateSeuratObject(counts = data1, project = "Mj_hem1", min.cells = 0, min.features = 100)
Mj2=CreateSeuratObject(counts = data2, project = "Mj_hem2", min.cells = 0, min.features = 100)
Mj3=CreateSeuratObject(counts = data3, project = "Mj_hem3", min.cells = 0, min.features = 100)
```

## Percent of MT
```{r}
Mj1[["percent.mt"]] <- PercentageFeatureSet(Mj1, pattern = "^MT-")
Mj2[["percent.mt"]] <- PercentageFeatureSet(Mj2, pattern = "^MT-")
Mj3[["percent.mt"]] <- PercentageFeatureSet(Mj3, pattern = "^MT-")
```

## Plotting
```{r}
plot1 <- FeatureScatter(Mj1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Mj1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plot3 <- FeatureScatter(Mj2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(Mj2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4
plot5 <- FeatureScatter(Mj3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot6 <- FeatureScatter(Mj3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot5 + plot6
```

## Select single cell
```{r}
Mj1e <- subset(Mj1, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 &  percent.mt < 10 & percent.mt > 1)
Mj2e <- subset(Mj2, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 &  percent.mt < 10 & percent.mt > 1)
Mj3e <- subset(Mj3, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 &  percent.mt < 10 & percent.mt > 1)
```

## Merge seurat objects
```{r}
hem.combined <- merge(x = Mj1e, y =c(Mj2e, Mj3e), add.cell.ids = c("Mj1e", "Mj2e", "Mj3e"), project = "Mj_hem")
```

## Notice the cell names now have an added identifier
```{r}
head(x = hem.combined@active.ident)
table(hem.combined@meta.data$orig.ident)
```

## Split object
```{r}
list <- SplitObject(hem.combined, split.by = "orig.ident")
list <- list[c("Mj_hem1", "Mj_hem2", "Mj_hem3")]
```

## SCTransform
```{r}
list <- lapply(X = list, FUN = SCTransform, method = "glmGamPoi", vars.to.regress = "percent.mt", variable.features.n = 3000)
```

## Integration
```{r}
features <- SelectIntegrationFeatures(object.list = list, nfeatures = 3000)
list <- PrepSCTIntegration(object.list = list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT", anchor.features = features)
combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
combined.sct <- RunPCA(combined.sct, verbose = FALSE)
combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:30, n.neighbors = 50L, min.dist = 0.2, n.components = 2L, umap.method = "uwot")
DimPlot(combined.sct, reduction = "umap", split.by = "orig.ident")
```

## FindCluster
```{r}
combined.sct <- FindNeighbors(combined.sct, dims = 1:30, k.param = 50, verbose = TRUE)
combined.sct <- FindClusters(combined.sct, verbose = TRUE, algorithm = 4, resolution = 0.4)
DimPlot(combined.sct, reduction = "umap")
```

## Normalize RNA data for visualization purposes
```{r}
DefaultAssay(combined.sct) <- "RNA"
combined.sct <- NormalizeData(combined.sct, verbose = FALSE)
```
## Rename
```{r}
new.cluster.ids <- c("Hem2", "Hem6", "Hem5", "Hem4", "Hem1", "Hem3")
names(new.cluster.ids) <- levels(combined.sct)
combined.sct <- RenameIdents(combined.sct, new.cluster.ids)
levels(combined.sct) <- c("Hem1", "Hem2", "Hem3", "Hem4", "Hem5", "Hem6")
```

## FindMarkers
```{r}
Markers <- FindAllMarkers(combined.sct, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 1, test.use = "bimod")
```

## Collect gene expression
```{r}
expression <- apply(combined.sct@assays[["RNA"]]@data,1,sum)
expression <- sort(expression, decreasing = TRUE)
write.csv(expression, file="expression.csv")
```

## File saving
```{r}
write.csv(Markers, file="Markers.csv")
saveRDS(Mj1e, file = "Mj1e.rds")
saveRDS(Mj2e, file = "Mj2e.rds")
saveRDS(Mj3e, file = "Mj3e.rds")
saveRDS(hem.combined, file = "hem.combined.rds")
saveRDS(combined.sct, file = "combined.sct.rds")
```
