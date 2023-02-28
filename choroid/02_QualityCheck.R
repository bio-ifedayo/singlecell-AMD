

# Prepare data for quality metrics ----------------------------------------

#load("choroid_cell.RData")


View(merge_max@meta.data)

# create a sample column using the rowname of the metadata
merge_max$sample <- rownames(merge_max@meta.data)

# split sample column
merge_max@meta.data <-
  separate(
    merge_max@meta.data,
    col = 'sample',
    into = c('Sample_no', 'Barcode'),
    sep = '_'
  )
## Read phenotype data to R
phen=read.delim2("https://raw.githubusercontent.com/pharmlovex/SingleCellAnalysis/main/phendata.txt", 
                 sep = ',', header = T)

phen = phen %>%
  filter(Sample.Name %in% c("GSM5676873","GSM5676875","GSM5676877",
                            "GSM5676879","GSM5676881","GSM5676883"))

# integrate phendata with merge_max
for (y in 1:nrow(merge_max@meta.data)){
  z=merge_max@meta.data$Sample_no[y]
  merge_max@meta.data$sex[y] = phen$sex[phen$Sample.Name==z]
  merge_max@meta.data$source_name[y] = phen$source_name[phen$Sample.Name==z]
  merge_max@meta.data$phenotype[y] = phen$Phenotype[phen$Sample.Name==z]
}


view(merge_max@meta.data)
# compute quality assessment metrics  for the data

# calculate mitochondrial percentage
merge_max$mitoPercent <- PercentageFeatureSet(merge_max, pattern='^MT-')
merge_max$riboPercent <- PercentageFeatureSet(merge_max, pattern="^RP[LS]")

## Get a df from the slot meta data from the seurat object 

merge_max@meta.data -> meta.data


## Add some columns to the meta.data 

# Add number of genes per UMI for each cell to metadata
meta.data$log10GenesPerUMI <- log10(meta.data$nFeature_RNA) / log10(meta.data$nCount_RNA)

# Compute percent mito ratio
#meta.data$mitoRatio <- PercentageFeatureSet(object = merge_max, pattern = "^MT-")
meta.data$mitoRatio <- meta.data$mitoPercent / 100

merge_max@meta.data <- meta.data 

# save(merge_max,
#           file=paste0(workdir,'/','merge_max.RData'))



# Clone obj to perform some manipulation 
merge_max -> data.nomalat

#
#*****23762 genes and 30587 cells in 4 samples in its sparse matrix Only nonezeros values are stored
#Compute the percentage of the largest gene 
apply(data.nomalat@assays$RNA@counts,2,max) -> data.nomalat$largest_count

apply(data.nomalat@assays$RNA@counts,2,which.max)-> data.nomalat$largest_index

rownames(data.nomalat)[data.nomalat$largest_index] -> data.nomalat$largest_gene

100 * data.nomalat$largest_count / data.nomalat$nCount_RNA -> data.nomalat$percent.Largest.Gene

data.nomalat$largest_gene -> merge_max$largest_gene
data.nomalat$percent.Largest.Gene -> merge_max$percent.Largest.Gene
rm(data.nomalat)

merge_max@meta.data -> meta.data 

# Rename columns
meta.data <- meta.data %>%
  dplyr::rename(nUMI = nCount_RNA,
                nGene = nFeature_RNA)


# Add metadata back to Seurat object
merge_max@meta.data <- meta.data


# Create .RData object to load at any time
save(merge_max, 
     file=paste0(paths,'/','merged_obj.RData'))

#################################################################
#load the merged_obj
library(ggplot2)

#load('merged_obj.RData')

merge_max@meta.data -> meta.data
## Exploring the quality of the data ..........................................


# Visualize the number of cell counts per sample
meta.data %>% 
  ggplot(aes(x=phenotype, fill=phenotype)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 20),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none") +
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 25)) +
  geom_text(aes(label= ..count..), stat = "count", vjust =2.0, size = 10)+
  ggtitle("Number of Choriod cells")
# Visualize the number UMIs/transcripts per cell
meta.data %>% 
  ggplot(aes(color=phenotype, x=nUMI, fill= phenotype)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 23000)



# Visualize the distribution of genes detected per cell via histogram
meta.data %>% 
  ggplot(aes(color=phenotype, x=nGene, fill= phenotype)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)


# filtering# Visualize the distribution of genes detected per cell via boxplot
meta.data %>% 
  ggplot(aes(x=phenotype, y=log10(nGene), fill=phenotype)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoPercent)) + 
  geom_point() + 
  scale_colour_gradient(low = "grey", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept =1000) +
  geom_hline(yintercept = 500) +
  facet_wrap(~phenotype)


meta.data %>% 
  ggplot(aes(color=phenotype, x=mitoPercent, fill=phenotype)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic()


# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, color = phenotype, fill=phenotype)) +
  geom_density(alpha = 0.2) +
  theme_classic()

VlnPlot(merge_max, features=c("nUMI","mitoPercent", "riboPercent","percent.Largest.Gene"))



# Cell level filtering ----------------------------------------------------

 
merge_filtered <- subset(merge_max, subset = nUMI > 1000 &
                           nUMI < 20000 &
                           nGene > 500 &
                           mitoPercent < 25 &
                           percent.Largest.Gene < 25 &
                           riboPercent < 45 )

# View quality after 
VlnPlot(merge_filtered, features=c("nUMI","mitoPercent", "riboPercent","percent.Largest.Gene"))


dim(merge_filtered)

#merged_seurat




# gene level filtering ----------------------------------------------------

 


# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = merge_filtered, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 50

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_AMDP <- CreateSeuratObject(filtered_counts, meta.data = merge_filtered@meta.data)




# Visualize the outcome  --------------------------------------------------

 
VlnPlot(filtered_AMDP, features=c("nUMI","mitoPercent", "riboPercent","percent.Largest.Gene"))

#filtered_AMDP@meta.data->meta

#meta$source_name[grep('Chorod', meta$source_name)]='Choroid'

#merge_max@meta.data <- meta.data

filtered_AMDP@meta.data %>% 
  ggplot(aes(x=phenotype, fill=phenotype)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 20),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none") +
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 25)) +
  geom_text(aes(label= ..count..), stat = "count", vjust =2.0, size = 10)+
  ggtitle("Number of Retina cells")


