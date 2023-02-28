
# Installation  -----------------------------------------------------------
library(devtools)
install_github("immunogenomics/presto")

#collectGarbage()

# Load Library ------------------------------------------------------------

suppressWarnings({lapply(c("Seurat",'WGCNA','ggplot2',
                           "dbplyr","cowplot","patchwork",
                           "AnnotationDbi", "GO.db", "preprocessCore", "impute", "igraph",
                           "tester"), library, character.only = T)})

library(hdWGCNA)
library(presto)
library(data.table)
library(dplyr)
# # Download seurat obj  -------------------------------------------------------
# if (FALSE) {
#   wget(url = c("https://amdproject-1.eu-central-1.linodeobjects.com/AMDIntregate.RData",
#                "https://amdproject-1.eu-central-1.linodeobjects.com/Choriod.RData"))
# }

# Standard Single-cell Processing  -----------------------------------------

load("AMDIntregate.RData")
Choriod <- seurat.integrated %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures=2000) %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  FindNeighbors(dims = 1:15) %>% 
  FindClusters(resolution = 0.02) %>% 
  RunUMAP(dims=1:10)

DimPlot(Choriod, reduction = "umap")


# Cluster Annotation ------------------------------------------------------

library(HGNChelper)

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")

# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# To prepare the gene set lets import the db

# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Eye" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)


# Lets assign the cell  types to each cluster 

es.max = sctype_score(scRNAseqData = Choriod@assays$integrated@scale.data, scaled = TRUE,
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# Merge by cluster

cL_results = do.call("rbind", lapply(unique(Choriod@meta.data$seurat_clusters),
                                     function(cl){
                                       es.max.cl = sort(rowSums(es.max[,row.names(Choriod@meta.data[Choriod@meta.data$seurat_clusters==cl, ])]),
                                                        decreasing = !0)
                                       head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(Choriod@meta.data$seurat_clusters==cl)),
                                            10)
                                       
                                     }))

sctype_scores = cL_results %>% 
  group_by(cluster) %>% 
  top_n(n=1,wt =scores)

# set low-confident (low ScType score) clusters to "unknown"


sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells /
                     4] = "Unknown"
# View result
View(sctype_scores[, 1:3])


# Lets visualize the assigned cell 

Choriod@meta.data$Cell_Identity = " "
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  Choriod@meta.data$Cell_Identity[Choriod@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(Choriod, reduction = "umap", label = T, repel = T, group.by = "Cell_Identity") + ggtitle("Choriod cell communities")
backChoriod <- Choriod
Choriod -> Choriod
# Create a vector list of the cell types

cell.types <- c("Fibroblasts","Immune cells","Muller cells","Immune cells","Muller cells",
                "Immune cells","Microglial cells","Endothelial cells")

names(cell.types) <- levels(Choriod)

Choriod <- RenameIdents(Choriod, cell.types)
# Add the cell ident to the metadata
Idents(Choriod) -> cell_types

Choriod <- Choriod
Choriod$cell_types <- cell_types
save(Choriod, file = "Choriodplus.RData")

# Custom function for visualisation

umap_theme <- function () {
  theme(axis.line = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.background = element_blank(), panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), plot.title = element_text(hjust = 0.5))
}

DimPlot(Choriod, reduction = "umap", label = T, repel = T, ) + umap_theme() + NoLegend() +
  ggtitle("Choriod cell types")



# Perform DEGs with wilcoxauROC -------------------------------------------
#load("Choriod.RData")
options(digits=2)
res <- wilcoxauc(Choriod, 'cell_types' , seurat_assay = 'RNA')

# get the topmarker 

res %>% 
  filter(logFC > log(1.2) & pct_in > 50 & padj < 0.05) %>% 
  group_by(group) %>% 
  arrange(desc(logFC), .by_group=T)  %>% 
  top_n(n=6, wt= logFC) -> res.dotplot

  filter(group == "Fibroblasts") -> res.cone

save(res.down, file = "cellspecificmarker.RData")

marker = top_markers(res, n=6, auc_min = 0.5, pct_in_min = 50)

all_markers<- marker %>%
  select(-rank) %>% 
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>%
  .[!is.na(.)]

p<- DotPlot(object = Choriod, features = all_markers)
p

theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



save(Choriod, file = "Choriod.RData")

# Visualize specific genes 
p1= FeaturePlot(Choriod, c("MEG3", "TTR", "LMO4"), ncol = 1)
p2 =VlnPlot(Choriod, features = c("MEG3", "TTR", "LMO4"), pt.size = 0)

p1 | p2 


# Visualize top genes in each cluster -------------------------------------

colors <- c("red","blue")
res.dotplot$feature <- as.factor(res.dotplot$feature)
DotPlot(object = Choriod, features = unique(res.dotplot$feature), cols = colors) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 25))+
  ggtitle("Cell Specific Marker genes")



