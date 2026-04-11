# Single-cell-RNA-sequencing
A workflow outlining the steps to analyze single-cell RNA sequencing data on mice samples infected with influenza A at various timepoints d02, d05, d08, and day14.

## Introduction
Influenza A virus (IAV) is a highly contagious pathogen that often causes the seasonal flu in humans. Once a person is infected through contact with respiratory droplets, the virus starts its infection through the nasal mucosa, the tissue that lines the nasal cavity. Consequently, influenza A begins to hijack the molecular machinery of epithelial cells, which is a primary component of the nasal mucosa, and starts to replicate. Epithelial cell are further subdivided into different types of cells that include goblet cells, brush cells, small granule cells, and basal cells (Petruson, Hansson and Karlsson, 1984a). In addition, olfactory neurons, responsible for smell, are also present (Petruson, Hansson and Karlsson, 1984b). As a result of infection, the viral replication triggers a localized immune response inside the nasal cavity that includes cell damage and apoptosis, recruitment of immune cells, and causes inflammation (Zhu et al., 2022). Immune cells such as nasal submucosal dendritic cells, monocytes, and neutrophils are recruited to the site of infection to deal with the viral infection. During an influenza infection, the initial innate immune system initially starts defending against the pathogen followed by the adaptive immune response. However, the molecular mechanisms, cellular mechanism, and transcriptional dynamics underlying the immune response when faced with an influenza infection remains elusive. Therefore, the overall goal of scRNA-seq analysis is to measure the entire transcriptome, the entire set of RNA, for each individual cell in a sample/tissue. In this analysis, scRNAseq analysis will be applied to *Mus musculus (M. musculus)* that were infected with influenza A. *M. musculus* models are well studied in the field in the context of influenza A research where the immune response is well documented and comparable to human immune responses (Thangavel and Bouvier, 2014). Specifically, this analysis will characterize the cellular heterogeneity of diverse cell types within the nasal mucosa and understand what genes are regulating each cell during an influenza infection. In addition, this analysis aims to identify transcriptional changes and differences across timepoints 0 (no infection or naive), 02, 05, 08, and 14 post infection in immune response in response to an influenza infection from a published dataset (Kazer et al., 2024).  Thus, understanding the biology of IAV infection can open new avenues on more effective on preventative or treatment methods such as the development of vaccinations or anti-viral drugs. 
	
Recent advancements in single cell RNA sequencing (scRNA-seq) have revolutionized transcriptomic studies and our understanding of cell types. Specifically, it allows scientists to study the transcriptional activities of a cell at a certain timepoint. While traditional bulk RNA sequencing has also been widely used, it homogenizes samples by averaging transcripts across all cell types thus, making it unable to detect elusive cell types or masking important biological differences between individual cells. Therefore, scRNAseq would be the most appropriate tool since the murine nasal mucosa contains heterogeneous cells such as goblet, brush, and olfactory neurons that may be missed if bulk RNA seq was employed instead. Notably, scRNA-seq has the advantage in gaining a deeper understanding of samples by capturing cellular heterogeneity and unmask rare elusive cell types. The standard workflow includes aligning raw reads to reference genome, quality control filtering, data normalization, feature variable selection, dimensionality reduction, clustering, cell annotation, and differential gene expression. Many bioinformatics tools exist to analyze scRNA-seq data. For example, Scanpy (Wolf, Angerer and Theis, 2018) is a Python-based toolkit while Seurat (Butler et al., 2018) is employed in RStudio or R. Both toolkits are widely used to implement scRNA-seq workflows. Benchmarking studies have reported that there are considerably differences between the two toolkits when defaults parameters are employed for differential gene expression analysis. For example, a study reported that using the PBMC 10k cell dataset, Seurat reported more significant marker genes than Scanpy, with Seurat reporting 16,866 genes and Scanpy reported 11,007 genes (Rich et al., 2024). Seurat has been reported to have better plotting functions with the ggplot2 package whereas Scanpy employs seaborn or matplotlib which are also good but not as extensive as ggplot2. However, Seurat may struggle with computationally demanding datasets whereas Scanpy excels instead. For example, employing the FindClusters() function in Seurat has reported to crash often while Scanpy’s Leiden implementation produced the same results (Rich et al., 2024). In addition, choosing which toolkit often depends on the user’s familiarity programming language. All in all, both toolkits have their own strengths and weaknesses, with the decision often depending on the dataset and features that the user wants to employ. Therefore, based on the comparison, Seurat will be employed.

When performing dataset normalization, there are many methods available. For example, Seurat has innate methods such as the log normalization method where the feature counts of each cell are divided by the total counts for the cell. Subsequently, the result is multiplied by the scale factor (default is 10,000) then natural log transformed. An alternate method called SCtransform uses a regularized negative binomial regression to normalize and stabilize the variance of scRNA-seq data (Hafemeister and Satija, 2019). Notably, it has been reported that the SCtransform normalization is superior to other statistical models such as the default normalization method or the Poisson regression where the latter  prefers sparser datasets and shallower sequencing depth (Choudhary and Satija, 2022). However, as sequencing technologies improve, the Poisson regression does not properly capture all the variances in scRNA sequencing data. In contrast, the negative binomial regression performs better at accounting for technical factors such as sequencing depth and unique molecular identifiers (UMI) differences between each cell (Hafemeister and Satija, 2019). Therefore, SCtransform will be used for data normalization.
	
The type of statistical methods also affects the results when performing differential expressed gene analysis (DEG) on scRNA-seq data. A benchmarking study has shown that a parametric hurdle model method, Model-based Analysis of Single-cell Transcriptomics (MAST) performs better for scRNA data than the non-parametric Wilcoxon rank sum test which is Seurat’s default testing method. Specifically, MAST had one of the best performances when it came to batch effects and moderate sequencing depth whereas the Wilcoxon test only exceled at shallow sequencing depth. In addition, the MAST method had better F1 scores ranging from 0.5 – 0.6 whereas the Wilcoxon test has a F1 score of less than 0.5 (Nguyen et al., 2023). Other than MAST, pseudobulk with DESeq2 method can also be applied for scRNAseq to reduce the amount of false positive rates. In fact, it is one of the most-well established statistical methods for scRNAseq. Therefore, both MAST and pseudobulk DESeq2 will be employed.
	
## Methodology
Based on previous comparisons and justifications, Seurat will be employed with SCtransform for normlization, statistical tests MAST, and pseudobulk DESeq2 used to determine the top significant genes for feature plots, and DEG, respectfully.  

### Quality Control, Clustering, and Annotation
The Seurat object containing the dataset and metadata was obtained from the link (https://aacgenomicspublic.blob.core.windows.net/public/seurat_ass4.rds). Briefly, mice were infected with the influenza A virus, sacrificed, and three separate nasal tissue samples: respiratory (RM), olfactory (OM), and lateral nasal gland (LNG) were harvested at different timepoints at day 0 (naïve), day 2, day 5, day 8, and day 14 for a total n = 3 biological replicates for each timepoint (Kazer et al., 2024). Both desktop RStudio (v 4.5.1) and the R (v4.5.0) module on the High-Performance Computing (HPC) provided by The Digital Research Alliance of Canada (Nibi cluster), were employed for this analysis, where the HPC was usef for more computationally intensive functions. 

The Seurat (v5.4.0) (Butler et al., 2018) package was employed for scRNA-seq analysis. Quality control was performed using the subset() function to remove cells that contains greater than 15% mitochondrial reads and less than 500 detected genes. Subsequently, the entire dataset was normalized using the SCTransform() function with the model formula y ~ log_umi and method = "glmGamPoi" from the glmGamPoi package (v1.22.0) (Ahlmann-Eltze and Huber, 2021) parameter to increase computational efficiency. The SCTransform() has been reported to be more accurate at normalizing scRNAseq data than the Normalize() function and encompassed finding variable features and data scaling (Choudhary and Satija, 2022).

Dimensionality reduction was performed using principal component analysis (PCA). An elbow plot was employed to determine the best PCA for dimensionality reduction. Different PCA values were tested as recommended by the Seurat’s vignette and a range of 10-50 PCAs were tested with PCA 36 chosen as the best PCA due to the over clustering observed at higher PCA dimensions. The function FindNeighbors() used the shared nearest-neighbour (SNN) with k.param = 20 (nearest neighbours). and Louvain algorithm was utilized for clustering. Clustering resolutions of 0.5 – 1.0 was used to test out different resolutions, and a resolution of 0.5 was utilized as the ideal resolution of single cells and to prevent over clustering. Finally, Uniform Manifold Approximation and Projection (UMAP) embedding was performed with RunUMAP() function with dims = 1:36. UMAP embedding was chosen over other methods such as t-Distributed Stochastic Neighbor Embedding (t-SNE) because UMAP is better at preserving global structure. 

### Cell Annotation and Feature Plots
	
Automated cell annotation was performed using the SingleR (v2.12.0) package with the database ImmGenData() from the celldex package (v1.20.0) (Aran et al., 2019). This database was used due to its well-known annotation of murine immune cell types and due to the appropriateness of the database in modelling infectious disease biology for this analysis. The parameters label.main, the more general label, was used to obtain broader and more robust cell type classification. Results were plotted as a cell annotated UMAP. Diagnostic plots were also employed to ensure that cells were correctly annotated. All markers for the dataset were determined with the FindAllMarkers() function from the Seurat package (Butler et al., 2018) with the following parameters: logfc.threshold = 0.1 to detect subtle gene changes, test.use = "MAST" as it handles scRNAseq datasets better by handling potential dropout events, max.cells.per.ident = 500 for downsampling to speed up computational time and random.seed = 123 to ensure reproducibility. Downsampling of 500 is often sufficient and robust to detect markers without losing statistical power. Subsequently, the average logfold of > 1 and adjusted p values < 0.05 were filtered to retain statistically significant genes. Top genes for cluster 2 as identified by the clustering steps before were plotted as feature plots. 

### DEG and GSEA analysis
Before performing DEG analysis, 56, 018 cells with missing mouse IDs were excluded. To address the high false positives rates associated with scRNAseq DEG analysis, pseudobulking was performed to improve statistical rigor with the AggregateExpression() and grouping cells by time, mouse id, the type of tissue and the Seurat cluster. Samples were subsetted by tissues and pairwise comparisons between day 0 (naïve) and post-infection timepoints of day 02, day 05, day 08, and day 14, were performed for a total of 12 comparisons. The FindMarkers() function was employed with parameters logfc.threshold = 0.1, a lower threshold required for further Gene Set Enrichment Analysis (GSEA) and test.use = "DESeq2" due to DESeq2’s better handling of bulk datasets. Plots were generated for all pairwise comparisons but only naïve and day 02 and day 14 will be shown to model early and late influenza A infection.

GSEA was performed with clusterProfiler package (v4.18.4) (Yu et al., 2012) using annotations from the *M. musculus* genome database. GSEA was employed since it is a more robust method than over-representation analysis (ORA) and it considers all ranked distributions of genes rather only top identified DEGs for ORA, making it a more statistically powerful method. Statistical analysis of GSEA was employed with the following parameters: p-values were adjusted using the Benjamini-Hochberg correction to control for false discovery rates (FDR), with adjusted p-values < 0.05 considered significant for gene and fast gene set enrichment analysis (fgsea) (Korotkevich, Sukhov and Sergushichev, 2019), a method used to identify enriched and significant pathways. To remove redundant GO terms, the simplify() from the clusterProfiler package was applied. For GSEA, genes were ranked according to the log2foldchange values obtained from Seurat’s FindMarker() function for each pairwise comparison. GO annotations were retrieved from the *M. musculus* database using Bioconductor’s annotation package, org.Mm.eg.db with a focus on biological process (BP) ontology. Results were plotted as GSEA enrichment plots.

## Results
### UMAP embedding reveals 39 distinct clusters and annotation of immune cells
To visualize the high-dimensional transcriptomic landscape, UMAP embedding with 36 PCA was performed. This analysis revealed 39 distinct clusters with the Louvain algorithm employed (Figure 1). Clusters 0 and 1 represented the clusters with the highest number of cells. The number of neighbours (20) and PCA dimensions of 36 were chosen to prioritize global structure of the dataset, and to avoid overclustering. The identified UMAP clusters served as the foundational clusters and will be used for further downstream analysis such as DEG and GSEA. 

Next, automatic cell annotation was employed using the Immgen database from the celldex package. Distinct cluster of cell types representing specialized immune cells such as neutrophils, B cells and macrophages were observed (Figure 2). Interestingly, microglia and macrophages were identified as separate cell types but as 1 cluster in Figure 1, highlighting the similarity in gene expression between the two cells. Indeed, macrophages are a specialized type of macrophage that are resident immune cells of the central nervous system. Notably, minimal over clustering was observed between different clusters. To rule out the potential of false positive cell annotations, a diagnostic plot was employed to assess the accuracy of the cell annotation process (Figure 3). Overall, by visual assessment, no major cell types were misidentified. For example, there are blocks that represent the correct cell types for both B cells and epithelial cells. In summary, clustering of cells and cell types were observed during the UMAP embedding and cell annotation.

<br>

<img width="4500" height="4500" alt="01_UMAP_pca_36_0 5_res" src="https://github.com/user-attachments/assets/d14ddb5c-29e2-49da-b647-b64ba621e379" />

**Figure 1. UMAP visualization of scRNA seq nasal mucosal dataset of *M. musculus* infected with influenza A displaying clustering.** 39 clusters were determined with 20 nearest neighbours and 36 PCAs. Each point is coloured accordingly by their clusters as shown on the legend.

<br>

<img width="6000" height="6000" alt="02_UMAP_cell_annotation" src="https://github.com/user-attachments/assets/89b828da-e958-4ac8-ac7c-b0221c42fcf2" />

**Figure 2. Cell type annotation of *M. musculus* nasal mucosa.** Briefly, the Immgen database from the celldex R package was used to annotate the entire dataset. Major cell populations identified included endothelial cells, and various types of immune cells such as B cells, T cells, macrophages and microglia. Cell types are denoted by different colours by the legend on the left.

<br>

<img width="1815" height="1381" alt="03_heatmap_diagnostics" src="https://github.com/user-attachments/assets/0f808af4-a8cf-407e-b7a6-878b4a99d3ab" />

**Figure 3. Diagnostic heatmap plot to validate proper cell annotations across identified cell clusters.** Overall, the cell clusters indicate correct identification. Distinct block expression ensures that the marker genes are specific to the identified cluster and ensures robust-cell type identification. The rows are labeled by cell type while the scores indicate which clusters were identified. Legend on the right indicates the expression as yellow indicates high expression and dark blue indicates low expression. Cell clusters are colour-coded matching the scores on the heatmap. 

<br>

### Feature plots reveal macrophage/microglial specific genes for cluster 02

Next, feature plots with the top 4 most significant genes with adjusted p value < 0.05 and log fold change > 0.1 were plotted on UMAP for cluster 02. Cluster 02 was chosen due to its macrophages/microglia annotation and due to their known roles during immune responses. Thefore, it may reveal relevant biology in the context of influenza A infection. The top four determined genes were *C1qa, C1qb, C1qc,* and *Cst3.* Notably, the first three genes, also known as Complement C1q A-chain, Complement C1q B-chain, and Complement C1q C-chain, are part of the complement pathway where they encode the three polypeptide chains A, B, and C that make up the C1q glycoprotein (Petry et al., 1996). *Cst3*, also known as Cystatin 3, is a gene highly expressed in microglia (Masuda et al., 2019). In fact, these feature plots results are consistent with the cell annotation from Figure 2 where this cluster was annotated as microglia or macrophages and thus validates the previous findings. While the genes in the compliment C1q pathway are mostly localized to cluster 2, *Cst3* expression is also expressed in other clusters suggesting that *Cst3* is not a specific enough marker for microglia and may play different roles in different cells. All in all, the feature plots presented genes that are characteristic of microglia and macrophages. 

<br>

<img width="1240" height="1108" alt="04_feature_plot_cluster_02" src="https://github.com/user-attachments/assets/a1819e0f-55f2-4d8f-9658-52c8a3417c1d" />

**Figure 4. Feature plots of the topmost significant genes in cluster 2.** Top genes included genes that play a major role in the classical C1q complement pathway and *Cst3*, a marker for a type of macrophage, microglia. Genes were filtered on adjusted p value < 0.05 and log fold change of > 1. Legend represents gene expression, with the darker blue colour indicating higher expression.

<br>

### Top DEGs during an influenza infection are related to ribosomal processes and stress response
DEG analysis was performed for all different pairwise comparisons between naïve and post influenza infection timepoints. Log fold threshold of > 0.1 was chosen to detect subtle DEGs. DEG analysis for all 3 different tissue types during the early D0 to D02 timepoints showed mostly genes related to ribosomal processes such as *Rpl39, Rp41, Rps21, Rsp16* (Figure 5 - 7). Interestingly, *Plek* which codes for the Pleckstrin protein, was one the most differentially expressed in RM samples which showed higher expression in day 02 in RM tissues (Figure 7). Notably, Plek has been reported to be highly upregulated in inflammatory diseases and is a regulator of inflammation (Alim et al., 2022). 

To rule out that ribosomal genes processes are mostly dominating the DEG analysis and are reflecting true biology, a separate analysis for the timepoint at day 14 and for comparison between early to late viral infection. Little to no ribosomal genes were detected as the top DEGs (Figure 5 - 7). In fact, the later timepoint had DEGs such as *Gadd45b, Nfkbiz, Lars2*, and *Ubb*, the first two genes which had higher expression in day 14 samples while *Ubb* had higher expression in the day 0 RM samples (Figure 5 B, 7B).  Specifically, *Gadd45b* codes for Growth Arrest And DNA Damage Inducible Beta and *Nfkbiz* codes for Nuclear Factor Kappa B Inhibitor Zeta. Notably, *Gadd45b* is known to regulate cellular stress response when an environmental stressor is present by mediating the NF-kappa pathway B (Papa et al., 2004; Salerno et al., 2012). This is in line with the biology of viral infections as cells are often in a stressed state after fighting and recovering from an acute infection/stressor such as influenza A.

Total number of DEGs for all comparisons with adjusted p values < 0.05 are summarized with the number of both downregulated and upregulated genes determined (Table 1). Interestingly, there were more DEGs in the RM tissue samples where there were more downregulated genes observed compared to upregulated. This may suggest a more dynamic transcriptional change in the RM tissues as compared to other tissue types. 

<br>

<img width="2080" height="1040" alt="05_Vlnplot_Naive_D02_LNG_combined" src="https://github.com/user-attachments/assets/17b192a4-80de-443a-8eed-8ba791e3f503" />

**Figure 5. Violin plots of the topmost significant genes in cluster 2 of LNG tissues.** Top DEGs expressed in A. day 02 samples and B. day 14 samples. Day 0 (naive) was used as the baseline for DEG. Genes were filtered on adjusted p value < 0.05 and log fold change of > 0.1.  X-axis represents the timepoints whereas the y-axis labels expression level. 

<br>

<img width="2080" height="1040" alt="05_Vlnplot_day02_day12_OM_Combined" src="https://github.com/user-attachments/assets/434800b2-5d96-4a72-99c8-f9105a7b61a8" />

**Figure 6. Violin plots of the topmost significant genes in cluster 2 of OM tissue samples.** Top DEGs expressed in A. day 02 samples and B. day 14 samples. Day 0 was used as the baseline for DEG. Genes were filtered adjusted p value < 0.05 and log fold change of > 0.1.  X-axis represents the timepoints whereas the y-axis labels expression level. 

<br>

<img width="2177" height="1090" alt="05_Vlnplot_D02_D14_RM_combined" src="https://github.com/user-attachments/assets/2b2e8394-d490-42f4-ad85-6b2647c8fd66" />

**Figure 7. Violin plots of the topmost significant genes in cluster 2 of RM tissue samples.** Top DEGs expressed in A. day 02 samples and B. day 14 samples. Day 0 was used as the baseline for DEG. Genes were filtered adjusted p value < 0.05 and log fold change of > 0.1.  X-axis represents the timepoints whereas the y-axis labels expression level.

<br>

| Tissue | Comparison | Total DEGs | Upregulated | Downregulated |
|--------|--------|--------|--------|--------|
|LNG     |D02	|137	|74		|63
|LNG     |D05	|107	|74		|33
|LNG	|D08	|110	|104	|6
|LNG	|D14	|67		|55		|12
|OM		|D02	|77		|70		|7
|OM		|D05	|92		|31		|61
|OM		|D08	|76		|70		|6
|OM		|D14	|42		|31		|11
|RM		|D02	|426	|89		|337
|RM		|D05	|364	|31		|333
|RM		|D08	|371	|57		|314
|RM		|D14	|486	|43		|443

**Table 1. Summary of the number of DEGs for each comparison for each tissue type.** All comparison was performed to day0 (naïve) as a baseline. A log fold change > 0.1 was used. Only significant genes with adjusted p valued < 0.05 were included in the table.

<br>

### GSEA show biological processes related to biosynthesis and immune response 
GSEA analysis with a focus on biological process was performed using the DEG identified from the FindMarkers() function and was employed between DEG determined between naïve and day 02 and day 14 for all three tissue types. Day 02 and day14 would represent early and late immune responses, respectively. A p value cutoff of 0.05 was used to filter for significant genes corrected with the Benjamini-Hochberg correction applied. GSEA revealed two main biological themes during early influenza infection: 1. biosynthesis and translation, and 2. immune response. Notably, all tissue samples had gene set IDs related to either metabolic processes or translation (Figure 8 A-C). In LNG and OM tissues but not RM, gene set IDs terms revealed biological processes related to immune response such as humoral response, cellular response to biotic stimulus, and cellular response to molecule of bacterial origin (Figure 8 A, B). The latter term is somewhat unexpected given that influenza A is a virus, not a bacterium. In day 14, immune processes are present in the context of late infection and healing such as cellular detoxication and B cell activation in OM and RM tissue samples (Figure 8E, F). All in all, GSEA revealed processes related to cells responding to an influenza infection.
<br>

<img width="6000" height="7500" alt="08_combined_gsea_plots" src="https://github.com/user-attachments/assets/2700ebb0-e7f1-4bc1-b40f-716d32320f4e" />

**Figure 8.  Gene set enrichment analysis (GSEA) enrichment plots between naïve and day02 post infection of Cluster 02 for all three tissue samples.** Enrichment plots for different tissue samples at A-C. day 02  and D-E. day 14. Curves show enrichment score for top 3 gene sets identified for each comparison. Each tick represents a gene that falls within the ranked list. The lower panel denotes the ranking metric of how each gene is ranked. Only the top 3 significant GO biological process gene sets were plotted. P value cutoff < 0.05 with the Benjamini-Hochberg correction.

<br>

## Discussion
scRNAseq has revolutionized the field of gene expression analysis and has allowed us to gain a deeper understanding of the transcriptome at the single cell level. In the event of an influenza infection, cells of the innate immune system such as macrophages are first recruited to the site of infection followed by the adaptive immune response. In this analysis, scRNA seq analysis was applied to three murine tissues in the nasal mucosa that were infected with IAV at 4 different timepoints, with day 02 modeling early infection and day 14 as the late infectious period.

### UMAP embedding and cell annotation of cluster 02
39 clusters were determined from UMAP embedding (Figure 1). From those 39 clusters, cell types were annotated (Figure 2). Notably, UMAP has been shown to preserve global structure better than t-SNE (t-distributed stochastic neighbor embedding) therefore comparisons between clusters can be inferred. For example, the macrophages and microglia cell annotation, denoted together as cluster 02, were observed to be in proximity to each other suggesting that they are similar biologically. In fact, microglia are a specialized type of macrophage that mainly reside in the central nervous system. In addition, the presence of microglia here is not surprising given that one of the tissue types include in the analysis is the olfactory tissue which is part of the brain. 

Feature plots further validated the identity of the macrophage and microglia cluster where genes that play a role in the classical C1q complement pathway and *Cst3* were expressed in this cluster (Figure 4). Notably, the complement system plays an important role in the innate immune system for the initial step of recognizing pathogens such as influenza viruses. Consequently, macrophages will often express C1q receptors that recognize these C1q proteins and promote phagocytosis of macrophages to engulf viral particles (Bohlson et al., 2014). Therefore, the DEGs of the complement pathway make sense in this context due to immune response to IAV. *Cst3* is a cysteine proteinase inhibitor that is shown to play a role in neurodegenerative diseases (Masuda et al., 2019). This finding is consistent with other scRNA-seq analysis where *Cst3* was also determined as a marker of microglia in both mice and human microglia (Masuda et al., 2019). However, the role of *Cst3* in the context of influenza A infection has been limited with most of the literature focused as a marker for microglia or on its role in neurodegenerative diseases such as Alzheimer’s (Keren-Shaul et al., 2017).

### DEGs in early vs late immune response
DEGs was performed for all pairwise timepoint comparisons with the naïve (day 0) acting as the baseline. In the day 02 timepoint, which represents the early stages of infection, ribosomal genes such as *Rsp21, Rsp27,* and *Rsp16* were differentially expressed for all three tissue types (Figure 5-7). An explanation for this observation is that viruses cannot replicate on its own and must take advantage of the host’s ribosomal machinery to replicate its own genetic material. Therefore, it is possible that macrophages are infected with IAV which ramps up ribosomal activity to produce more influenza A viral copies. It has been reported that the certain influenza viruses such as H1N1 and H5N1, can infect human macrophages and increase inflammatory response through Type I interferon (Lee et al., 2009).

During the later timepoint of day 14, genes associated with stress were differentially expressed such as *Ubb, Gadd45b* and *Nfkbiz*. Notably, *Gadd45b* is known to mediate the NF-kappa pathway B under environmental stress (Papa et al., 2004; Salerno et al., 2012). This is consistent with the biology of viral infections as cells are often in a stressed state after fighting and recovering from an acute infection/stressor such as influenza A. Interestingly, *Ubb*, known as Ubiquitin B, is a gene that codes for the Ubiquitin B protein which targets proteins for degradation. An explanation of the higher expression of Ubb in the earlier timepoints is caused by the influenza virus where the virus hijacks the degradation pathway to remove all non-viral proteins such as the host interferon response to prevent an immune response to prioritize its own replication and promoting their survival. In contrast, the host cell might attempt to ubiquitinate the assembled viral proteins to removal the viral proteins and protect itself (Park et al., 2022). Future studies are warranted to distinguish between the two mechanisms under the context of viral infection.

### GSEA 
GSEA was performed where it revealed biological processes related to biosynthesis of molecules and immune responses (Figure 8). An explanation of the biosynthesis of molecules terms is that the virus needs to replicate and thus needs to ramp up the production of molecules to produce more viral particles. Literature has reported that PA-X and nonstructural protein 1 (NS1) are proteins that belongs to IAV and aids with IAV to replicate. Notably, PAX-1 can target host mRNA for degradation, and block host mRNA processing (Hu, Ma and Liu, 2018) while NS1 ensures that the viral mRNA is exported to be translated by the ribosomes (Pott, Kryvenko and Vadász, 2026). Notably, other types of IAV such as H1N1 employs the PAX-1 mediated method  (Hu, Ma and Liu, 2018). Consequently, as a response to viral replication, immune responses are activated. Interestingly, the term B cell activation was observed in the later timepoint in RM tissue samples (Figure 8 F). This was a surprising finding despite analysis being conducted in cluster 2 which is characterized by macrophages or microglia, both being part of the innate immune response. Nevertheless, this is consistent with current known biology where cells in the adaptive immune response such as B cells, play a role during later viral infection. Notably, B cells produce antibodies specific to the virus that help macrophages better detect pathogens to clear up the infection and establish long-term immunity. In summary, GSEA suggest processes related to ribosomal translation and immune responses. 

### Future Work
In this analysis, only automated cell annotation was employed. Future steps should also use manual annotation in conjecture with automatic annotation to increase the confidence that the right types of cells are classified. In addition, it would be interesting to perform a timepoint trajectory analysis to unveil the progression of influenza A and to observed which genes are differentially expressed. In addition, future steps should implement a multiomics approach where other levels of biology such as the proteome and the metabolome. Having a cohesive multi-approach analysis on multiple levels of biology will entail a robust understanding of influenza infection and its rewiring effects on the nasal tissues.

## References
Ahlmann-Eltze, C. and Huber, W. (2021) “glmGamPoi: fitting Gamma-Poisson generalized linear models on single cell count data,” Bioinformatics, 36(24), pp. 5701–5702. Available at: https://doi.org/10.1093/bioinformatics/btaa1009.

Alim, M.A. et al. (2022) “Pleckstrin Levels Are Increased in Patients with Chronic Periodontitis and Regulated via the MAP Kinase-p38α Signaling Pathway in Gingival Fibroblasts,” Frontiers in Immunology, Volume 12-2021. Available at: https://doi.org/10.3389/fimmu.2021.801096.

Aran, D. et al. (2019) “Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage,” Nature Immunology, 20(2), pp. 163–172. Available at: https://doi.org/10.1038/s41590-018-0276-y.

Bohlson, S.S. et al. (2014) “Complement, C1q, and C1q-Related Molecules Regulate Macrophage Polarization,” Frontiers in Immunology, Volume 5-2014. Available at: https://doi.org/10.3389/fimmu.2014.00402.

Butler, A. et al. (2018) “Integrating single-cell transcriptomic data across different conditions, technologies, and species,” Nature Biotechnology, 36(5), pp. 411–420. Available at: https://doi.org/10.1038/nbt.4096.

Choudhary, S. and Satija, R. (2022) “Comparison and evaluation of statistical error models for scRNA-seq,” Genome Biology, 23(1), p. 27. Available at: https://doi.org/10.1186/s13059-021-02584-9.

Hafemeister, C. and Satija, R. (2019) “Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression,” Genome Biology, 20(1), p. 296. Available at: https://doi.org/10.1186/s13059-019-1874-1.

Hu, J., Ma, C. and Liu, X. (2018) “PA-X: a key regulator of influenza A virus pathogenicity and host immune responses,” Medical Microbiology and Immunology, 207(5), pp. 255–269. Available at: https://doi.org/10.1007/s00430-018-0548-z.

Kazer, S.W. et al. (2024) “Primary nasal influenza infection rewires tissue-scale memory response dynamics,” Immunity, 57(8), pp. 1955-1974.e8. Available at: https://doi.org/10.1016/j.immuni.2024.06.005.

Keren-Shaul, H. et al. (2017) “A Unique Microglia Type Associated with Restricting Development of Alzheimer’s Disease,” Cell, 169(7), pp. 1276-1290.e17. Available at: https://doi.org/10.1016/j.cell.2017.05.018.

Korotkevich, G., Sukhov, V. and Sergushichev, A. (2019) “Fast gene set enrichment analysis,” bioRxiv, p. 060012. Available at: https://doi.org/10.1101/060012.

Lee, S.M.Y. et al. (2009) “Systems-Level Comparison of Host-Responses Elicited by Avian H5N1 and Seasonal H1N1 Influenza Viruses in Primary Human Macrophages,” PLOS ONE, 4(12), p. e8072. Available at: https://doi.org/10.1371/journal.pone.0008072.

Masuda, T. et al. (2019) “Spatial and temporal heterogeneity of mouse and human microglia at single-cell resolution,” Nature, 566(7744), pp. 388–392. Available at: https://doi.org/10.1038/s41586-019-0924-x.

Nguyen, H.C.T. et al. (2023) “Benchmarking integration of single-cell differential expression,” Nature Communications, 14(1), p. 1570. Available at: https://doi.org/10.1038/s41467-023-37126-3.

Papa, S. et al. (2004) “Gadd45β mediates the NF-κB suppression of JNK signalling by targeting MKK7/JNKK2,” Nature Cell Biology, 6(2), pp. 146–153. Available at: https://doi.org/10.1038/ncb1093.

Park, E.-S. et al. (2022) “The Roles of Ubiquitination in Pathogenesis of Influenza Virus Infection,” International Journal of Molecular Sciences, 23(9), p. 4593. Available at: https://doi.org/10.3390/ijms23094593.

Petruson, B., Hansson, H.-A. and Karlsson, G. (1984b) “Structural and Functional Aspects of Cells in the Nasal Mucociliary System,” Archives of Otolaryngology, 110(9), pp. 576–581. Available at: https://doi.org/10.1001/archotol.1984.00800350018006.

Pott, J., Kryvenko, V. and Vadász, I. (2026) “The ribosomal landscape in influenza A virus infection: from molecular mechanisms to clinical relevance,” European Respiratory Review, 35(179), p. 250049. Available at: https://doi.org/10.1183/16000617.0049-2025.

Rich, J.M. et al. (2024) “The impact of package selection and versioning on single-cell RNA-seq analysis,” bioRxiv, p. 2024.04.04.588111. Available at: https://doi.org/10.1101/2024.04.04.588111.

Salerno, D.M. et al. (2012) “Gadd45a and Gadd45b modulate innate immune functions of granulocytes and macrophages by differential regulation of p38 and JNK signaling,” Journal of Cellular Physiology, 227(11), pp. 3613–3620. Available at: https://doi.org/10.1002/jcp.24067.

Thangavel, R.R. and Bouvier, N.M. (2014) “Animal models for influenza virus pathogenesis, transmission, and immunology,” Journal of Immunological Methods, 410, pp. 60–79. Available at: https://doi.org/10.1016/j.jim.2014.03.023.

Wolf, F.A., Angerer, P. and Theis, F.J. (2018) “SCANPY: large-scale single-cell gene expression data analysis,” Genome Biology, 19(1), p. 15. Available at: https://doi.org/10.1186/s13059-017-1382-0.

Yu, G. et al. (2012) “clusterProfiler: an R Package for Comparing Biological Themes Among Gene Clusters,” OMICS: A Journal of Integrative Biology, 16(5), pp. 284–287. Available at: https://doi.org/10.1089/omi.2011.0118.

Zhu, F. et al. (2022) “H1N1 Influenza Virus-Infected Nasal Mucosal Epithelial Progenitor Cells Promote Dendritic Cell Recruitment and Maturation,” Frontiers in Immunology, Volume 13-2022. Available at: https://doi.org/10.3389/fimmu.2022.879575.

