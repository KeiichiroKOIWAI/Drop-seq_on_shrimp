# Pseudotime analysis
## Reading of libraries
```{r}
library(needs)
needs(dplyr,Seurat,ggplot2,SeuratWrappers,SeuratData,patchwork,broom,sinaplot,monocle3,glmGamPoi)
```
## Opening data set
```{r}
library(needs)
needs(dplyr,Seurat,ggplot2,SeuratWrappers,SeuratData,patchwork,broom,sinaplot,monocle3,glmGamPoi)
combined.sct <- readRDS(file = "combined.sct.rds")
Markers <-read.csv(file = "Markers.csv", header = TRUE)
Mj1e <- readRDS(file = "Mj1e.rds")
Mj2e <- readRDS(file = "Mj2e.rds")
Mj3e <- readRDS(file = "Mj3e.rds")
```
## G2M/S phase calculation
```{r}
G2M <- scan(file="analysis/G2M_S/G2M_list.txt", what="character", sep=NULL)
S <- scan(file="analysis/G2M_S/S_list.txt", what="character", sep=NULL)

combined.sct <- CellCycleScoring(
  combined.sct,
  s.features = S,
  g2m.features = G2M
  )
DimPlot(combined.sct, group.by = "Phase", split.by = "orig.ident")
```
## Pseudotime analysis
```{r}
cds <- as.cell_data_set(combined.sct)
cds <- cluster_cells(cds, reduction_method="UMAP", k = 50)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

integrated.sub <- subset(as.Seurat(cds), monocle3_partitions == 1)
cds <- as.cell_data_set(integrated.sub)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
```
