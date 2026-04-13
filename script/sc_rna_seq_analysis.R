#
#TITLE: scRNA-seq Analysis 
#
#BY: Liona Vu

################ 1.0 Load libraries and EDA  #################### 

library(Seurat)
library(sctransform)
library(dplyr)
library(ggplot2)
#BiocManager::install('glmGamPoi')
library(glmGamPoi)
library(future)
library(SingleR)
library(celldex)
library(openxlsx)
library(MAST)
library(cowplot)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library(patchwork)
library(aplot)

#run this first!!!!!
options(future.globals.maxSize = 128 * 1024^3)

#increase vector memory limit to 200GB'
mem.maxVSize(vsize = 200000)

# Run the following when on HPC cluster but do not run before sctransform(), will end up in a zombie-like state and the script will be stuck
#plan("multicore", workers = 8) 

#For interactive session when debugging on HPC
#salloc --time=1:00:00 --mem=185G --cpus-per-task=8 --account=def-cottenie

#load file
seurat_file <-readRDS(file = "seurat_ass4.rds")

#EDA
head(seurat_file, n = 20)
dim(seurat_file)
#25129 features, 156572 cells
table(seurat_file@meta.data$orig.ident)
table(seurat_file@meta.data$biosample_id)


################ 2.0 Quality Control  #################### 

#gets the number of mitochondrial genes starting with mt (mice genes starts with lowercase)
seurat_file[["percent.mt"]] <- PercentageFeatureSet(seurat_file, pattern = "^mt-")
seurat_file$percent.mt

# Visualize QC metrics as a violin plot
VlnPlot(seurat_file, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3, pt.size = 0) 
ggsave("figures/QC_violin_plots.png", height = 10, width = 18)

#plot scatter to see feature relationships closely
plot_1 <- FeatureScatter(seurat_file, feature1 = "nCount_RNA", feature2= "percent.mt")
plot_1
ggsave("figures/feature_scatter_rnavsmt.png", height = 10, width = 10)

plot_2 <- FeatureScatter(seurat_file, feature1 = "nCount_RNA", feature2= "nFeature_RNA")
plot_2
ggsave("figures/feature_scatter_rna_count_feature.png", height = 10, width = 10)

#Based on the Vlnplots and scatter plots, will remove reads that have less than 500 genes, greater than 15% mtDNA
subset_file <- subset(seurat_file, subset = nFeature_RNA > 500 & percent.mt < 15)
saveRDS(subset_file, file = "output/01_subset_file.RDS", compress = FALSE)
dim(subset_file)

#Plot Vlnplot after subsetting
VlnPlot(subset_file, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave("figures/QC_violin_plots_filtered.png", height = 10, width = 18)

################ 3.0 Data transformation and PCA ##############################

#Basic SCTansform normalization of the data
#subset_file <- readRDS(file = "output/01_subset_file.RDS")
subset_file_sctransformed <- SCTransform(subset_file, verbose = TRUE, method = "glmGamPoi")

saveRDS(subset_file_sctransformed, file = "output/02_subset_file_sctransformed.RDS", compress = FALSE)

#The other option Normalize data
#subset_file_trans <- NormalizeData(subset_file, normalization.method = "LogNormalize", scale.factor = 10000)
#Identification og highly feature variables
#subset_file_trans <- FindVariableFeatures(subset_file_trans, selection.method = "vst", nfeatures = 2000)

#see the top 10 features 
#top10 <- head(VariableFeatures(subset_file_trans), n =10)

#Plot variable features plot with the top 10 features
#plot_vf <- LabelPoints(plot = VariableFeaturePlot(subset_file_trans), points = top10, repel = TRUE)
#ggsave(plot_vf, file = "figures/variable_features_plot.png", height = 10, width = 10)

#all_genes <- rownames(subset_file_trans)

#Scale data
#subset_file_scaled <- ScaleData(subset_file_trans, features = all_genes)

#write RDS file
#saveRDS(subset_file_scaled, file = "output/subset_file_scaled.RDS", compress = FALSE)
#subset_file_scaled <- readRDS("subset_file_scaled.rds")

#run a PCA to produce principal components that can be used to cluster our cells
pca_cluster <- RunPCA(subset_file_sctransformed, features = VariableFeatures(object = subset_file_sctransformed))
saveRDS(pca_cluster, file = "output/03_pca.RDS", compress = FALSE) #set to FALSE, because it will take a longer time
      
#Read in PCA RDS file
#pca_cluster <- readRDS("output/03_pca_cluster.RDS")

#Find best number of dimensions
elbow_plot <- ElbowPlot(pca_cluster, ndims = 50)
#From above, will try out PC at different intervals 10, 20, 30, 40, and 50 as per Seurat vignette's recommendation)

ggsave("figures/elbow_plot.png", height =10, width =10)

################ 4.0 Systematic Analysis for best PCA and Resolution ##############################

#Downstream analysis, make function to test out different PCA, neighbours, clustering, and make umap.
best_pca <- function(pca_result, number) {
  number <- as.numeric(number)
  
  neighbours <- FindNeighbors(pca_result, dims = 1:number)
  
  clustering <- FindClusters(neighbours, resolution = 0.5)
  
  umap <- RunUMAP(clustering, dims = 1:number)
  saveRDS(umap, file = paste0("output/umap_", number, ".RDS"))
  
  dimplot <- DimPlot(umap, reduction = "umap", label = TRUE)
  
  return(dimplot)
}

#Run for 10, 15, 20, 30, 40, and 50 PCA to compare and contrast
pca_10 <- best_pca(pca_cluster, 10)
ggsave(pca_10, file = "figures/dimplot_pca_10.png", height = 10, width = 10)

pca_20 <- best_pca(pca_cluster, 20)
ggsave(pca_20, file = "figures/dimplot_pca_20.png", height = 10, width = 10)

pca_30 <- best_pca(pca_cluster, 30)
ggsave(pca_30, file = "figures/dimplot_pca_30.png", height = 10, width = 10)

pca_36 <- best_pca(pca_cluster, 15)
ggsave(pca_36, file = "figures/dimplot_pca_36.png", height = 10, width = 10)

pca_40 <- best_pca(pca_cluster, 40)
ggsave(pca_40, file = "figures/dimplot_pca_40.png", height = 10, width = 10)

pca_50 <- best_pca(pca_cluster, 50)
ggsave(pca_50, file = "figures/dimplot_pca_50.png", height = 10, width = 10)

#Combine PCA
pca_combined <- grid.arrange(pca_10, pca_15, pca_20, pca_30, pca_40, pca_50, ncol = 3)
ggsave(pca_combined, file = "figures/pca_combined.png", height = 15, width = 20)

#Based on PCA above, will chose 36 as the ideal PCA for analysis
#Determine best resolution for PCA 36 at 0.6, 0.8, 1.0
pca_36 <- FindNeighbors(pca_cluster, dims = 1:36)

pca_36_0.6 <- FindClusters(pca_36, dims = 1:36, resolution = 0.6)
pca_36_0.6 <- RunUMAP(pca_36_0.6, dims = 1:36)
pca_36_umap_0.6 <- DimPlot(pca_36_0.6, reduction = "umap", label = TRUE)
ggsave(pca_36_umap_0.6, file = "figures/pca_36_0.6_umap.png", height = 15, width = 15)

pca_36_0.8 <- FindClusters(pca_36, dims = 1:36, resolution = 0.8)
pca_36_0.8 <- RunUMAP(pca_36_0.8, dims = 1:36)
pca_36_umap_0.8 <- DimPlot(pca_36_0.8, reduction = "umap", label = TRUE)
ggsave(pca_36_umap_0.8, file = "figures/pca_36_0.8_umap.png", height = 15, width = 15)

pca_36_1.0 <- FindClusters(pca_36, dims = 1:36, resolution = 1.0)
pca_36_1.0 <- RunUMAP(pca_36_1.0, dims = 1:36)
pca_36_umap_1.0 <- DimPlot(pca_36_1.0, reduction = "umap", label = TRUE)
ggsave(pca_36_umap_1.0, file = "figures/pca_36_1.0_umap.png", height = 15, width = 15)

#Based on the resolution above, resolution of 0.5 or 0.6 are good, with 0.6 having a bit more overclustering, will use res 0.5 instead.

#Make final plot for UMAP for PCA36
DimPlot(umap_pca_36, reduction = "umap", label = TRUE, repel = TRUE) +
  labs(x = "UMAP 1", y = "UMAP 2", title = "UMAP Clustering Plot")
ggsave(file = "figures/pca_36_0.5_umap.png", height = 15, width = 15)

#umap_pca_36 <- readRDS("output/umap_36.RDS")
#plot3 <- DimPlot(umap_pca_36, reduction = "umap", group.by = "organ_custom")
#ggsave(plot3, file = "figures/before_integration_organ.png", height = 10, width = 10)
#plot4 <- DimPlot(umap_pca_36, reduction = "umap", group.by = "time")
#ggsave(plot4, file = "figures/before_integration_time.png", height = 10, width = 10)

#clustering by organ? might try to integrate
#umap_pca_36[["RNA"]] <- split(umap_pca_36[["RNA"]], f = umap_pca_36$organ_custom)
#umap_pca_36 <- SCTransform(umap_pca_36, method = "glmGamPoi", verbose = TRUE)
#umap_pca_36 <- RunPCA(umap_pca_36, verbose = TRUE)

# Integrate
#umap_pca_36 <- IntegrateLayers(object = umap_pca_36, 
              #                 method = HarmonyIntegration, 
             #                  orig.reduction = "pca", 
              #                 new.reduction = "harmony", 
              #                 verbose = FALSE)

#umap_pca_36[["RNA"]] <- JoinLayers(umap_pca_36[["RNA"]])

# Re-cluster and viz
#umap_pca_36 <- FindNeighbors(umap_pca_36, reduction = "harmony", dims = 1:36)
#umap_pca_36 <- FindClusters(umap_pca_36, resolution = 0.5)
#umap_pca_36 <- RunUMAP(umap_pca_36, dims = 1:36, reduction = "harmony")

#saveRDS(umap_pca_36, "output/04_UMAP_integrated.RDS", compress = FALSE)

#plot <- DimPlot(umap_pca_36, reduction = "umap", group.by = "organ_custom")
#ggsave(plot, file = "figures/UMAP_integrated_organ_custom.png", height = 10, width =10)

#plot2 <- DimPlot(umap_pca_36, reduction = "umap")
#ggsave(plot2, file = "figures/UMAP_integrated.png", height = 10, width =10)

#Integration did not seem to work as intended, still separated by organ so might continue with original.

################ 5.0 Cell Type Annotation ########################
umap_pca_36 <- readRDS("output/umap_36.RDS")

# Annotation tutorial from: https://www.youtube.com/watch?v=7RuPGaWcY0Y
# Get reference dataset
# ref <- celldex::MouseRNAseqData()
ref <- celldex::ImmGenData()
View(as.data.frame(colData(ref)))

# Run single R
umap_pca_counts <- GetAssayData(umap_pca_36, layer = 'counts')

pred <- SingleR(test = umap_pca_counts,
                ref = ref, 
                labels = ref$label.main)

rownames(umap_pca_36@meta.data)
rownames(pred)

#Add cells to seurat object based on matching rownames
umap_pca_36$singleR_labels <-pred$labels[match(rownames(umap_pca_36@meta.data), rownames(pred))]
table(umap_pca_36@meta.data$singleR_labels)

#plot UMAP with annotated cells
DimPlot(umap_pca_36, reduction = "umap", group.by = "singleR_labels", label = TRUE, repel = TRUE) +
  ggtitle("Cell Annotation") +
  labs(x = "UMAP 1", y = "UMAP 2")
ggsave("figures/UMAP_cell_annotation.png", height = 20, width = 20)

#Annotation Diagnostics
# scores within cells
plotScoreHeatmap(pred)
ggsave("figures/heatmap_diagnostics.png", height = 20, width = 20)

#delta values
plotDeltaDistribution(pred)

################ 6.0 Feature Plots ############################################ 

# Tutorial link: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

#W try find all DE markers for all clusters
all_markers <- FindAllMarkers(umap_pca_36, only.pos = TRUE, logfc.threshold = 0.1,
                              test.use = "MAST", 
                              max.cells.per.ident = 500, #Downsampling to speed up computational time
                              random.seed = 123) #Set seed for reproducibility

#filter by cluster on avg logfold change greater than 1 and adj p values less than 0.05
filtered_all_markers <- all_markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1 & p_val_adj < 0.05)

#save to excel sheet in case
write.xlsx(filtered_all_markers, file = "output/all_markers.xlsx")

#Write function to filter by cluster,number of genes and plot feature plot
plot_feat_plot <- function(dataframe, cluster_n, genes, umap) {
  cluster_n = as.numeric(cluster_n)
  genes = as.numeric(genes)
  
  top_genes_cluster <- dataframe %>%
    filter(cluster == cluster_n) %>%
    head(n =genes) %>%
    pull()
  
 plot <- FeaturePlot(umap, features = top_genes_cluster) &
    labs(x = "UMAP 1", y = "UMAP 2")
 
  return(plot)
}

#Feature Plot for cluster 0
plot_feat_plot(filtered_all_markers, 0, 4, umap_pca_36)

#Feature Plot for cluster 2
plot_feat_plot(filtered_all_markers, 2, 4, umap_pca_36)
ggsave(file = "04_featurePlot_cluster02.png", height = 10, width = 10)
#immune cells that express complement genes, microglia??

#Feature Plot for cluster 11
#plot_feat_plot(filtered_all_markers, 11, 6, umap_pca_36)
#myeloid cells, neutrophil activation? 
#s100a9, Il1b, s100a8, Tyrobp, Srgn, Cxcr2, S100a11, Cebpb, Macl1

#plot_feat_plot(filtered_all_markers, 12, 6, umap_pca_36)
#epithelial cells

#Curious about heatmap, choose top 5 genes
#filtered_all_markers %>%
 # group_by(cluster) %>%
  #slice_head(n = 5) %>%
  #ungroup() -> top10

#DoHeatmap(umap_pca_36, features = top10$gene) + NoLegend()


################ 7.0 Pseudobulking and DEG ############################################
#Picked cluster 2 for analysis
# Need to do pseudobulking to correct for false postives
# From before, there seems to be mouse id missing from metadata, needs to be removed
table(umap_pca_36@meta.data$mouse_id) #56018 missing!

#remove samples without mouse IDs
umap_pca_36_filt <- subset(umap_pca_36, subset = mouse_id != "")

#Sanity check
table(umap_pca_36_filt@meta.data$mouse_id) # 0 now
table(umap_pca_36_filt@meta.data$seurat_clusters) 
table(umap_pca_36_filt@meta.data$time) 

# Aggregation of data 
pseudo_umap_36 <- AggregateExpression(umap_pca_36_filt, assays = "RNA", return.seurat = T, group.by = c("mouse_id", "time", "organ_custom", "seurat_clusters"))
pseudo_umap_36@meta.data$organ_custom

#pseudo_umap_36$cluster.time <- paste(pseudo_umap_36$seurat_clusters, pseudo_umap_36$time, sep = "_")
#Idents(pseudo_umap_36) <- "cluster.time"

#subset by tissue types
pseudo_umap_36_OM <- subset(pseudo_umap_36, subset = organ_custom == "OM")
pseudo_umap_36_RM <- subset(pseudo_umap_36, subset = organ_custom == "RM")
pseudo_umap_36_LNG <- subset(pseudo_umap_36, subset = organ_custom == "LNG")

#Subset by cluster
pseudo_umap_36_OM_2<- subset(pseudo_umap_36_OM, subset = seurat_clusters == "2")
pseudo_umap_36_RM_2<- subset(pseudo_umap_36_RM, subset = seurat_clusters == "2")
pseudo_umap_36_LNG_2<- subset(pseudo_umap_36_LNG, subset = seurat_clusters == "2")

pseudo_umap_36_OM_2$time <- pseudo_umap_36_OM_2$time
pseudo_umap_36_RM_2$time <- pseudo_umap_36_RM_2$time
pseudo_umap_36_LNG_2$time <- pseudo_umap_36_LNG_2$time

Idents(pseudo_umap_36_OM_2) <- "time"

test2 <- FindMarkers(pseudo_umap_36_OM_2, ident.1 = "Naive", 
                   ident.2 = "D05",
                   logfc.threshold = 0.1, #Lower threshold for GSEA
                   test.use = "DESeq2") 

#create vector of timepoints to be looped over
tp <- c("D02", "D05", "D08", "D14")

#Initiate an empty list
list_OM_comparison <- list()

#Might do a loop for all posible comparisons for each tissue
for (i in tp) {
  deg <- FindMarkers(pseudo_umap_36_OM_2, ident.1 = i, # condition
                                   ident.2 = "Naive", #baseline
                                   logfc.threshold = 0.1, #Lower threshold for GSEA
                                   test.use = "DESeq2") #better for pseudobulk
  
  #Append to list
  list_OM_comparison[[i]] <- deg
}

head(list_OM_comparison)


#For RM
Idents(pseudo_umap_36_RM_2) <- "time"

list_RM_comparison <- list()

for (i in tp) {
  deg <- FindMarkers(pseudo_umap_36_RM_2, ident.1 = i, 
                            ident.2 = "Naive",
                            logfc.threshold = 0.1, #Lower threshold for GSEA
                            test.use = "DESeq2") #better for pseudobulk
  
  list_RM_comparison[[i]] <- deg

}

#For LNG
Idents(pseudo_umap_36_LNG_2) <- "time"
list_LNG_comparison <- list()

for (i in tp) {
  deg <- FindMarkers(pseudo_umap_36_LNG_2, ident.1 = i, 
                            ident.2 = "Naive",
                            logfc.threshold = 0.1, #Lower threshold for GSEA
                            test.use = "DESeq2") #better for pseudobulk
  
  list_LNG_comparison[[i]] <- deg
}

#find numver of significant up and down genes
#no_deg <- function (list_input) {

#make loop to count # of up and down regulated, change the list 
for (i in tp) {
  filtered <- list_RM_comparison[[i]] %>%
    filter(p_val_adj < 0.05) #filter for sig adj adj values
  print(nrow(filtered)) #gets total # of genes
    for (n in names(filtered)) {
      up <- sum(filtered[n]$avg_log2FC>0)
      print(paste0("Up ", up))
      down <- sum(filtered[n]$avg_log2FC <0)
      print(paste0("Down ", down))
    }
}  
#LNG, 74, 63, not sure why the 0s appeared but at lest we got the #DEGs correct

#Grab genes, most significant genes are at the top, change day to have dif timepoints
top_genes_LNG <- head(rownames(list_LNG_comparison[["D02"]]), n = 4) 

#plotting DEG as violin plots
#Reorder the samples, relevel for x axis so naive starts first
levels(pseudo_umap_36_LNG_2) <- c("Naive", tp)
vlnplot_LNG <- VlnPlot(pseudo_umap_36_LNG_2, features = top_genes_LNG, ncol = 2)

# Add title
vlnplot_LNG + plot_annotation(title = "A. Top differentially expressed genes from Naive to D02 in LNG samples",
                              theme = theme(plot.title = element_text(face = "bold")))

#Save
ggsave2(file= "figures/Vlnplot_Naive_D02_LNG.png", height = 10, width = 10)

#For OM
top_genes_OM <- head(rownames(list_OM_comparison$D02), n = 4)

levels(pseudo_umap_36_OM_2) <- c("Naive", tp)
vlnplot_OM <- VlnPlot(pseudo_umap_36_OM_2, features = top_genes_OM, ncol = 2)

#plot Violin plot
vlnplot_OM + plot_annotation(title = "A. Top differentially expressed genes from Naive to D02 in OM samples",
                theme = theme(plot.title = element_text(face = "bold")))


ggsave2(file= "figures/Vlnplot_Naive_D02_OM.png", height = 10, width = 10)

# For RM
top_genes_RM <- head(rownames(list_RM_comparison$D02), n = 4)

levels(pseudo_umap_36_RM_2) <- c("Naive", tp)
vlnplot_RM <- VlnPlot(pseudo_umap_36_RM_2, features = top_genes_RM, ncol = 2)

vlnplot_RM + plot_annotation(title = "A. Top differentially expressed genes from Naive to D02 in RM samples",
                             theme = theme(plot.title = element_text(face = "bold")))

ggsave2(file= "figures/Vlnplot_Naive_D02_RM.png", height = 10, width = 10)

################ 8.0 GSEA Analysis ##################################################

#For all tissues types, create function for data wrangling and GSEA
gsea_analysis <- function(input) {
gene_list_df <- as.data.frame(input)
gene_list_df <- na.omit(gene_list_df)

gene_list <- gene_list_df$avg_log2FC

gene_list_df$gene_names <- rownames(gene_list_df)

#Name the list of vectors
names(gene_list) <- gene_list_df$gene_names

#Make named sorted from greater to least, gseGO does not take in dataframes!
gene_list <- sort(gene_list, decreasing = TRUE)

#check keytypes
#keytypes(org.Mm.eg.db) #Use ALIAS

#perform GSEA
gsea_results <- gseGO(geneList = gene_list,
                      ont = "BP",
                      OrgDb = org.Mm.eg.db, 
                      keyType = "ALIAS",
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      verbose =  TRUE, 
                      eps = 0, #P values less than 1.0 e-10
                      by = "fgsea",
                      seed = TRUE)

#Simplify to remove redundant genes
gsea_clust <- clusterProfiler::simplify(x = gsea_results, 
                                            cutoff = 0.7,
                                            by = "p.adjust",
                                            select_fun = min)
return(gsea_clust)
}

#Function to plot
plot_gsea <- function (df_input, day, tissue) {
  gsea <- gsea_analysis(df_input)
  
  gsea <- gseaplot2(gsea, geneSetID = c(1,2,3), 
                                 color = c("#E495A5", "#86B875", "#7DB0DD"),
                                 pvalue_table = FALSE, 
                                 title = paste0("Naive to ", day,  " Influenza A Infection in ", tissue), 
                                 base_size = 15)
  ggsave2(gsea, file = paste0("figures/08_gseaplot_naive_", day, "_", tissue,".png"), height = 10, width = 12)
  
  return(gsea)
}

#Apply function
set.seed(420)
gsea_d02_LNG <- plot_gsea(list_LNG_comparison[["D02"]], "D02", "LNG")
gsea_d02_OM <- plot_gsea(list_OM_comparison[["D02"]], "D02", "OM")
gsea_d02_RM <- plot_gsea(list_RM_comparison[["D02"]], "D02", "RM")

gsea_d14_LNG <- plot_gsea(list_LNG_comparison[["D14"]], "D14", "LNG")
gsea_d14_OM <- plot_gsea(list_OM_comparison[["D14"]], "D14", "OM")
gsea_d14_RM <- plot_gsea(list_RM_comparison[["D14"]], "D14", "RM")

#Combine and save plots
plot_list(gglist = list(gsea_d02_LNG, gsea_d02_OM, gsea_d02_RM,
                        gsea_d14_LNG, gsea_d14_OM, gsea_d14_RM), 
          byrow = T, ncol = 2, tag_levels = "A")

ggsave2(file = "figures/08_combined_gsea_plots.png", height = 25, width = 20)


