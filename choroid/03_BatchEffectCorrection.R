
# perform standard workflow steps to figure out if we see any batch effects --------

## Normalize the data 
merged_seurat_filtered <- NormalizeData(object = filtered_AMDP)
##Feature selection
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered)
## Scaling data
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)
## Dimensional Reduction
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)

ElbowPlot(merged_seurat_filtered)
##Cell Sub population Identification using unsupervised clustering methods [!prior information based]
## Techniques for unsupervised clustering 
## (i) k-means; (ii) hierarchical clustering; (iii) density-based clustering; and (iv) graph-based clustering 
## 
## Other cluster methods are 
## single-cell consensus clustering (SC3) (Kiselev et al., 2017) 
## shared nearest neighbor (SNN) clustering algorithm (Satija et al., 2015) ....Seurat[findNeighbors]
##SC3 is an unsupervised approach that combines multiple clustering approaches, 
## which has a high accuracy and robustness in single-cell clustering.
merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered)
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:20)


# plot
p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'phenotype') +
  ggtitle("Batch Effect")
p2 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'source_name',
              cols = c('red','green'))

grid.arrange(p1, p2, ncol = 2, nrow = 2)
#pca plot 
t1 <- DimPlot(merged_seurat_filtered, reduction = 'pca', group.by = 'phenotype')
t2 <- DimPlot(merged_seurat_filtered, reduction = 'pca', group.by = 'source_name',
              cols = c('red','green'))

## Get the list of most highly expressed genes
apply(merged_seurat_filtered@assays$RNA@data,1,mean) -> gene.expression

sort(gene.expression, decreasing = TRUE) -> gene.expression

head(gene.expression, n=10)

save(filtered_AMDP, file = 'AMD_Project.RData')



# Batch Effect Correction  ------------------------------------------------

#**************NORMALISATION & BATCH EFFECT CORRECTION*****************

# perform integration to correct for batch effects ------
#obj.list <- SplitObject(merged_seurat_filtered, split.by = 'phenotype')

obj.list <- SplitObject(filtered_AMDP, split.by = 'phenotype')

for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}
# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)

# find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features)
# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)


# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:50)


p3 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'phenotype') + 
  ggtitle("Corrected Batch effect")
p4 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'source_name',
              cols = c('red','green'))


grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)


t3=DimPlot(seurat.integrated, reduction = 'pca', group.by = 'phenotype')
t4=DimPlot(seurat.integrated, reduction = 'pca', group.by = 'source_name',
           cols = c('red','green'))

## View all the pca plots
grid.arrange(t1, t3, t2, t4, ncol = 2, nrow = 2)

p1|p3
#save

save(seurat.integrated, file = 'AMDIntregate.RData')


