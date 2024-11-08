###################################
##### GEO accession:GSE151974 #####
###################################
library(Seurat)
library(SeuratObject)
library(clustree)
library(dplyr)
library(tidyverse)
library(patchwork)
options(stringsAsFactors = FALSE)

### Data importing ###
raw_count <- read.csv("GSE151974_raw_umi_matrix_postfilter.csv", row.names = 1) # Download: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151974/suppl/GSE151974%5Fcell%5Fmetadata%5Fpostfilter.csv.gz
metadata <- read.csv("GSE151974_cell_metadata_postfilter.csv", row.names = 1) # Download: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151974/suppl/GSE151974%5Fraw%5Fumi%5Fmatrix%5Fpostfilter.csv.gz

sce <- CreateSeuratObject(counts = raw_count,
                          meta.data = metadata)


### P7 seperation ###
Idents(sce) <- "CellType"
unique(sce$CellType)
subset <- subset(sce, idents = c("Cap-a", 
                                 "Art", 
                                 "Cap", 
                                 "Vein", 
                                 "Col13a1+ fibroblast", 
                                 "Col14a1+ fibroblast", 
                                 "Myofibroblast",
                                 "SMC"))
Idents(subset) <- "Sample"
P7_subset <- subset(subset, idents = c("P7_Hyperoxia", "P7_Normoxia"))


### ProcessExper ###
source("./processExper.R")
P7_integrated = processExper(P7_subset, 
                             sct.method = T, 
                             reduction = "cca")


### Cell clustering and dimensional reduction ###
P7_integrated <- RunPCA(P7_integrated)
dims = 1:40
P7_integrated <- RunUMAP(P7_integrated,
                         dims = dims,
                         reduction.name = "umap") %>%
                 FindNeighbors(dims = dims) %>%
                 FindClusters(resolution = c(seq(0, 1, .1)))

clustree(P7_integrated) # Select suitable resolution

P7_integrated <- FindClusters(P7_integrated, resolution = 0.5)


### Cell annotation ###
p <- DotPlot(P7_integrated, 
             assay = "SCT",
             group.by = "seurat_clusters",
             features = c("Gpihbp1", "Kit", # gCap
                          "Car4", "Kdr", # aCap
                          "Cxcl12", "Pcsk5", # Art
                          "Vegfc", "Prss23", # Vein
                          "Pecam1", "Eng", "Cd34", "Cdh5", # Gen ECs
                          "Col1a1", "Col1a2", "Col3a1", "Fn1", "Tagln", "Acta2", "Myl9", "Myh11", # Mesenchyme
                          "Tgfbi","Wnt5a" #Myofibroblast
                         ) 
            ) +
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10),
        panel.background = element_rect(color = "black"))+
  coord_flip()

p1 <- DimPlot(P7_integrated, 
              reduction = "umap", 
              group.by = "seurat_clusters", 
              label = T)
wrap_plots(p + p1) # Reference

cluster_ids <- c("gCap",          #cluster 0
                 "Fibroblast",    #cluster 1
                 "Fibroblast",    #cluster 2
                 "aCap",          #cluster 3
                 "Myofibroblast", #cluster 4
                 "Myofibroblast", #cluster 5
                 "Fibroblast",    #cluster 6
                 "EndoMT",        #cluster 7
                 "Fibroblast",    #cluster 8
                 "Fibroblast",    #cluster 9
                 "Vein",          #cluster 10
                 "Art",           #cluster 11
                 "SMC",           #cluster 12
                 "Fibroblast",    #cluster 13
                 "Myofibroblast", #cluster 14
                 "Fibroblast"     #cluster 15
                )

names(cluster_ids) <- levels(P7_integrated)
P7_integrated <- RenameIdents(P7_integrated, cluster_ids)
P7_integrated$Celltype_fine <- Idents(P7_integrated)

Idents(P7_integrated) <- "seurat_clusters"
cluster_ids1 <- c("Endothelium",   #cluster 0
                  "Fibroblast",    #cluster 1
                  "Fibroblast",    #cluster 2
                  "Endothelium",   #cluster 3
                  "Myofibroblast", #cluster 4
                  "Myofibroblast", #cluster 5
                  "Fibroblast",    #cluster 6
                  "EndoMT",        #cluster 7
                  "Fibroblast",    #cluster 8
                  "Fibroblast",    #cluster 9
                  "Endothelium",   #cluster 10
                  "Endothelium",   #cluster 11
                  "SMC",           #cluster 12
                  "Fibroblast",    #cluster 13
                  "Myofibroblast", #cluster 14
                  "Fibroblast"     #cluster 15
)

names(cluster_ids1) <- levels(P7_integrated)
P7_integrated <- RenameIdents(P7_integrated, cluster_ids1)
P7_integrated$Celltype_main <- Idents(P7_integrated)


### Cell atlas visualization ###
source("./cellatlas_umap.R")
cluster_atlas <- cellatlas_umap(P7_integrated, 
                                idents = "seurat_clusters",
                                levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"), 
                                hull_alpha = 0.1, 
                                hull_size = 0.5, 
                                hull_lty = 2, 
                                hull_delta = 0.8, 
                                dot_color = c("#f38989", "#f9a341", "#f48521", "#ef6a45", "#549da3", "#96cb8f", "#f9b769", "#ec5051", "#a4cde1", "#67a4cc","#a48cbe","#e32422","#b2db87","#4dae47","#5c9e43","#b79973"),  # A list of color codes
                                dot_size = 0.5, 
                                dot_alpha = 1, 
                                label_size = 4, 
                                label_color = F 
)
cluster_atlas

celltype_atlas <- cellatlas_umap(P7_integrated, 
                                 idents = "Celltype_fine",
                                 levels = c("gCap",
                                            "aCap",
                                            "Art",
                                            "Vein",
                                            "EndoMT",
                                            "Fibroblast",
                                            "Myofibroblast",
                                            "SMC"), 
                                 hull_alpha = 0.1, 
                                 hull_size = 0.5, 
                                 hull_lty = 2,
                                 hull_delta = 0.8,
                                 dot_color = c("#ea5c6f", "#f7905a", "#e187cb", "#fb948d", "#e2b159", "#ebed6f", "#b2db87", "#7ee7bb"),  
                                 dot_size = 0.5,
                                 dot_alpha = 1, 
                                 label_size = 4, 
                                 label_color = F
                                 )
celltype_atlas

oxygen_atlas <- cellatlas_umap(P7_integrated, 
                               idents = "Oxygen",
                               levels = c("Normoxia", "Hyperoxia"),
                               hull_alpha = 0,
                               hull_size = 0, 
                               hull_lty = 2, 
                               hull_delta = 0, 
                               dot_color = c("#4DBBD5FF", "#E64B35FF"),  
                               dot_size = 0.5, 
                               dot_alpha = 1, 
                               label_size = 4, 
                               label_color = T
                               )
oxygen_atlas


### Dotplot of cell markers ###
source("./sce_dotplot.R")
marker = c("Gpihbp1", "Kit", # gCap
           "Car4", "Kdr",  # aCap
           "Cxcl12", "Pcsk5", # Art
           "Vegfc", "Prss23", # Vein
           "Pecam1", "Eng", "Cd34", "Cdh5", # General Endothelium
           "Col1a1", "Col1a2", "Col3a1", "Fn1", # Fibroblast
           "Tagln", "Acta2", "Myl9", "Myh11", # SMC
           "Tgfbi", "Wnt5a" # Myofibroblast
          )

P7_dotplot <- sce_dotplot(P7_integrated,
                          assay = "SCT",
                          idents = "Celltype_fine",
                          markers = marker,
                          levels = c("gCap",
                                     "aCap",
                                     "Art",
                                     "Vein",
                                     "EndoMT",
                                     "Fibroblast",
                                     "Myofibroblast",
                                     "SMC"
                                    ) 
                          )

P7_dotplot


### Dotplot of metabolic genes ###
marker1 = c("Adh5", "Aldh3a2", "Gapdh", "Gpi1", "Ldha", "Pkm", "Slc2a1", # Glycolysis
            "G6pdx", "Pgd", "Taldo1", "Tkt", # Pentose phosphate pathway
            "Dhfr", "Mthfd1", "Mthfd2", "Shmt1", "Shmt2", "Slc25a32", # One carbon metabolism
            "Acat1", "Aldh7a1", "Gcdh",  # Amino acid metabolism
            "Cpt1a", "Cpt1c", "Cpt2", # Carnitine shuttle
            "Acaa2", "Acadm", "Acads", "Acadsb", "Acadvl", "Acox3", "Echs1", "Eci1", "Eci2", "Hadh", "Hadha","Hadhb", # β-oxidation
            "Acly", "Acsl1", "Acsl4", "Acaca", "Degs1", "Elovl6", "Fasn", "Scd1", "Slc25a1", "Thrsp", # Fatty acid synthesis
            "Apoe", "Lipa", "Lipe", "Lpl", "Mttp", "Pnpla2", # Lipid metabolism and transport
            "Cyc1", "Ndufs1", "Ndufs2", "Ndufs3", "Pdha1", "Sdha" # Oxidative phosphorylation
            )


P7_EndoMT_dotplot <- sce_dotplot(subset(P7_integrated, idents = c("EndoMT")),
                                 assay = "SCT",
                                 idents = "Oxygen",
                                 markers = marker1,
                                 levels = c("Normoxia", "Hyperoxia"),
                                 title = "EndoMT",
                                 title_size = 10
                                )

P7_EndoMT_dotplot

P7_EndoMT_integrated <- subset(P7_integrated, Celltype_fine == "EndoMT")
P7_EndoMT_Hyperoxia_dotplot <- sce_dotplot(subset(P7_EndoMT_integrated, Oxygen == "Hyperoxia"),
                                           assay = "SCT",
                                           idents = "Sex",
                                           markers = marker1,
                                           levels = c("Female", "Male"),
                                           title = "EndoMT (Hyperoxia)",
                                           title_size = 10
                                          )

P7_EndoMT_Hyperoxia_dotplot

### Vlnplot of EndoMT markers ###
source("./Split_Vln_stacked.R")
levels = c("gCap",
           "aCap",
           "Art",
           "Vein",
           "EndoMT",
           "Fibroblast",
           "Myofibroblast",
           "SMC")
P7_integrated$Celltype_fine <- factor(P7_integrated$Celltype_fine, levels = levels) 
Idents(P7_integrated) <- "Celltype_fine"

Split_Vln_stacked(P7_integrated,
                  assay = "SCT",
                  feature = c("Acta2", "Myl9", "Tagln",
                              "Cdh5", "Eng", "Pecam1"),
                  split.plot = F,
                  pt.size = 0, 
                  face = "bold",
                  text_x_size = 12,
                  text_y_size = 10,
                  title_y_size = 14,
                  cols = c("#ea5c6f", "#f7905a", "#e187cb", "#fb948d", "#e2b159", "#ebed6f", "#b2db87", "#7ee7bb"),
                  test = F)


### Trajectory analysis ###
library(slingshot)

P7_integrated1 <- as.SingleCellExperiment(P7_integrated, assay = "SCT")

sds <- slingshot(P7_integrated1, 
                 clusterLabels = "Celltype_main",
                 reducedDim = "UMAP" ,
                 start.clus = "Endothelium"
)

source("./pseudotime_umap.R")
pseu1 <- pseudotime_umap(P7_integrated,
                         pseudotime_obj = sds,
                         idents = "Celltype_fine",  
                         sub_idents = "slingPseudotime_1", 
                         dot_color = c("#440154FF", "#FDE725FF"),
                         dot_size = 0.5, 
                         dot_alpha = 1, 
                         label_size = 4
)
pseu1

pseu2 <- pseudotime_umap(P7_integrated,
                         pseudotime_obj = sds,
                         idents = "Celltype_fine",  
                         sub_idents = "slingPseudotime_2", 
                         dot_color = c("#440154FF", "#FDE725FF"),
                         dot_size = 0.5, 
                         dot_alpha = 1, 
                         label_size = 4
)
pseu2

### Volcano plot ###
source("./FindMarker_genes.R")
source("./sce_volcanoplot.R")

P7_integrated <- PrepSCTFindMarkers(P7_integrated)

Idents(P7_integrated) <- "Celltype_fine"
DEG <- FindMarker_genes(object = subset(P7_integrated, Oxygen == "Normoxia"), 
                        assay = "SCT",
                        clusters = c("gCap", "aCap", "Art", "Vein", "EndoMT", "Fibroblast", "Myofibroblast", "SMC"),
                        comparison = c("Sex", "Female", "Male"),
                        logfc.threshold = 0,  
                        min.cells.group = 1
                       )  

DEG1 <- FindMarker_genes(object = subset(P7_integrated, Oxygen == "Hyperoxia"), 
                         assay = "SCT",
                         clusters = c("gCap", "aCap", "Art", "Vein", "EndoMT", "Fibroblast", "Myofibroblast", "SMC"),
                         comparison = c("Sex", "Female", "Male"),
                         logfc.threshold = 0,  
                         min.cells.group = 1
                        )  

levels <- c("gCap", "aCap", "Art", "Vein", "EndoMT", "Fibroblast", "Myofibroblast", "SMC")
volcano1 <- sce_volcanoplot(DEG,
                            levels = levels, 
                            title = "Male vs Female (Nox)", 
                            group_col =  c("#e74a32","#0da9ce"),
                            cluster_col = c("#ea5c6f", "#f7905a", "#e187cb", "#fb948d", "#e2b159", "#ebed6f", "#b2db87", "#7ee7bb"),
                            ptype = "pvalue"
                            )

volcano2 <- sce_volcanoplot(DEG1,
                            levels = levels, 
                            title = "Male vs Female (Hox)", 
                            group_col =  c("#e74a32","#0da9ce"),
                            cluster_col = c("#ea5c6f", "#f7905a", "#e187cb", "#fb948d", "#e2b159", "#ebed6f", "#b2db87", "#7ee7bb"),
                            ptype = "pvalue"
                            )
wrap_plot(volcano1/volcano2)

### GSEA analysis ###
library(msigdbr)
library(tibble)
library(patchwork)

c2_mmu <- msigdbr(species = "Mus musculus", category = c("C2")) %>% as.data.frame %>% dplyr::select(gs_name,entrez_gene,gene_symbol) %>% as.data.frame
c5_mmu <- msigdbr(species = "Mus musculus", category = c("C5")) %>% as.data.frame %>% dplyr::select(gs_name,entrez_gene,gene_symbol) %>% as.data.frame

c2_KEGG_mmu_pathway <- c2_mmu[grep(c("KEGG"), c2_mmu$gs_name),]
c2_PID_mmu_pathway <- c2_mmu[grep(c("PID_"), c2_mmu$gs_name),]
c2_BIOCARTA_mmu_pathway <- c2_mmu[grep(c("BIOCARTA"), c2_mmu$gs_name),]
c2_REACTOME_mmu_pathway <- c2_mmu[grep(c("REACTOME"), c2_mmu$gs_name),]
c2_WIKIPATHWAYS_mmu_pathway <- c2_mmu[grep(c("WIKIPATHWAYS"), c2_mmu$gs_name),]
c5_GOBP_mmu_pathway <- c5_mmu[grep(c("GOBP"), c5_mmu$gs_name),]

c2_KEGG_mmu_fgsea_sets <- c2_KEGG_mmu_pathway %>% split(x = .$gene_symbol, f = .$gs_name)
c2_PID_mmu_fgsea_sets <- c2_PID_mmu_pathway %>% split(x = .$gene_symbol, f = .$gs_name)
c2_BIOCARTA_mmu_fgsea_sets <- c2_BIOCARTA_mmu_pathway %>% split(x = .$gene_symbol, f = .$gs_name)
c2_REACTOME_mmu_fgsea_sets <- c2_REACTOME_mmu_pathway %>% split(x = .$gene_symbol, f = .$gs_name)
c2_WIKIPATHWAYS_mmu_fgsea_sets <- c2_WIKIPATHWAYS_mmu_pathway %>% split(x = .$gene_symbol, f = .$gs_name)
c5_GOBP_mmu_fgsea_sets <- c5_GOBP_mmu_pathway %>% split(x = .$gene_symbol, f = .$gs_name)

mmu_fgsea_sets <- c(c2_KEGG_mmu_fgsea_sets, 
                    c2_PID_mmu_fgsea_sets,
                    c2_BIOCARTA_mmu_fgsea_sets,
                    c2_REACTOME_mmu_fgsea_sets,
                    c2_WIKIPATHWAYS_mmu_fgsea_sets,
                    c5_GOBP_mmu_fgsea_sets
)

Idents(P7_integrated) <- "Celltype_fine"
DEG2 <- FindMarker_genes(dataset = P7_integrated, 
                         assay = "SCT",
                         clusters = c("gCap", "aCap", "Art", "Vein", "EndoMT", "Fibroblast", "SMC", "Myofibroblast"),
                         comparison = c("Oxygen", "Normoxia", "Hyperoxia"),
                         logfc.threshold = 0,  
                         min.cells.group = 1
                        )  

source("./sce_GSEA.R")
gsea_res <- sce_GSEA(DEG2, pathway = mmu_fgsea_sets)
gsea_res1 <- do.call(rbind, gsea_res) %>% t()
#write.csv(gsea_res1, file = "./P7_Hyperoxia_vs_Normoxia_MSigDB_C2CP_GSEA_results.csv", row.names = F)

DEG3 <- FindMarker_genes(object = subset(P7_integrated, Oxygen == "Normoxia"), 
                         assay = "SCT",
                         clusters = c("EndoMT"),
                         comparison = c("Sex", "Female", "Male"),
                         logfc.threshold = 0,  
                         min.cells.group = 1
                        )  

gsea_res2 <- sce_GSEA(DEG3, pathway = mmu_fgsea_sets)
gsea_res3 <- do.call(rbind, gsea_res2) %>% t()
#write.csv(gsea_res3, file = "./P7_EndoMT_Normoxia_Male_vs_Female_MSigDB_C2CP_GSEA_results.csv", row.names = F)

DEG4 <- FindMarker_genes(object = subset(P7_integrated, Oxygen == "Hyperoxia"), 
                         assay = "SCT",
                         clusters = c("EndoMT"),
                         comparison = c("Sex", "Female", "Male"),
                         logfc.threshold = 0,  
                         min.cells.group = 1
                         )  

gsea_res4 <- sce_GSEA(DEG4, pathway = mmu_fgsea_sets)
gsea_res5 <- do.call(rbind, gsea_res4) %>% t()
#write.csv(gsea_res5, file = "./P7_EndoMT_Hyperoxia_Male_vs_Female_MSigDB_C2CP_GSEA_results.csv", row.names = F)



### GSEA sigPathway visualization ###
source("./sce_GSEAbarplot.R")
sigPathway <- read.csv("./P7_EndoMT_Hyperoxia_vs_Normoxia_GSEA_sigPathway.csv", header = T)

levels = c("Carbohydrate metabolism",
           "Energy metabolism",
           "Lipid metabolism",
           "Signal transduction",
           "Signaling molecules and interaction",
           "Cell growth and death",
           "Immune system"
)
bar <- sce_GSEAbarplot(sigPathway, 
                       levels = levels[1:7], 
                       maxPathway = c("Carbohydrate catabolic process",
                                      "ADP metabolic process",
                                      "Membrane lipid metabolic process",
                                      "Positive regulation of Erk1 and Erk2 cascade",
                                      "Cytokine mediated signaling pathway",
                                      "p53 signaling pathway",
                                      "Interleukin-4 and interleukin-13 signaling"
                       ),
                       category_color = "grey10",
                       pathway_color = c("Carbohydrate metabolism" = "#66C2A5", 
                                         "Energy metabolism" = "#8DA0CB", 
                                         "Lipid metabolism" = "#A6D854", 
                                         "Signal transduction" = "#6a3d9a",
                                         "Signaling molecules and interaction" = "#E5C494",
                                         "Cell growth and death" = "#FC8D62",
                                         "Immune system" = "#E78AC3"),
                       title = "Hyperoxia vs. Normoxia (EndoMT)",
                       title_size = rel(1),
                       num_size = 0,
                       text_x_size = rel(1),
                       text_y_size = rel(1),
                       xlim = c(0, 2.3))

bar

sigPathway1 <- read.csv("./P7_EndoMT_Normoxia_Male_vs_Female_MSigDB_C2CP_GSEA_sigPathway.csv", header = T)

levels <- c("Signal transduction",
            "Immune system")

bar1 <- sce_GSEAbarplot(sigPathway1, 
                        levels = levels[1:2], 
                        maxPathway = c("Tgfb receptor signaling in EMT",
                                       "Positive regulation of IL4 production"),
                        title = "Male vs Female (EndoMT)",
                        category_color = "grey10",
                        pathway_color = c("Signal transduction" = "#6a3d9a",
                                          "Immune system" = "#E78AC3"),
                        title_size = rel(1),
                        num_color = "black",
                        num_size = 0,
                        text_x_size = rel(1),
                        text_y_size = rel(1),
                        xlim = c(-2.7, 3.2))

bar1

sigPathway2 <- read.csv("./P7_EndoMT_Hyperoxia_Male_vs_Female_MSigDB_C2CP_GSEA_sigPathway.csv", header = T)

levels <- c("Signal transduction",
            "Immune system")

bar2 <- sce_GSEAbarplot(sigPathway2,
                       levels = levels[1:2],
                       maxPathway = c("SMAD2/3/4 transcriptional activity",
                                      "Positive regulation of IL6 production"),
                       title = "Male vs Female (EndoMT)",
                       category_color = "grey10",
                       pathway_color = c("Signal transduction" = "#6a3d9a",
                                         "Immune system" = "#E78AC3"),
                       title_size = rel(1),
                       num_color = "black",
                       num_size = 0,
                       text_x_size = rel(1),
                       text_y_size = rel(1),
                       xlim = c(-2.5, 2.5))

bar2

### P14 seperation ###
P14_subset <- subset(subset, idents = c("P14_Hyperoxia", "P14_Normoxia"))

### ProcessExper ###
P14_integrated = processExper(P14_subset, 
                              sct.method = T, 
                              reduction = "cca")


### Cell clustering and dimensional reduction ###
P14_integrated <- RunPCA(P14_integrated)
dims = 1:40
P14_integrated <- RunUMAP(P14_integrated,
                          dims = dims,
                          reduction.name = "umap") %>%
                  FindNeighbors(dims = dims) %>%
                  FindClusters(resolution = c(seq(0, 1, .1)))

clustree(P14_integrated) # Select suitable resolution

P14_integrated <- FindClusters(P14_integrated, resolution = 0.5)

### Cell annotation ###
p2 <- DotPlot(P14_integrated, 
              assay = "SCT",
              group.by = "seurat_clusters",
              features = c("Gpihbp1", "Kit", # gCap
                           "Car4", "Kdr", # aCap
                           "Cxcl12", "Pcsk5", # Art
                           "Vegfc", "Prss23", # Vein
                           "Pecam1", "Eng", "Cd34", "Cdh5", # Gen ECs
                           "Col1a1", "Col1a2", "Col3a1", "Fn1", "Tagln", "Acta2", "Myl9", "Myh11", # Mesenchyme
                           "Tgfbi","Wnt5a" #Myofibroblast
                           ) 
             ) +
      theme(axis.title = element_blank(),
            axis.line = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(size = 10),
            panel.background = element_rect(color = "black"))+
      coord_flip()

p3 <- DimPlot(P14_integrated, 
              reduction = "umap", 
              group.by = "seurat_clusters", 
              label = T)
wrap_plots(p2 + p3) # Reference



cluster_ids2 <- c("Fibroblast",    #cluster 0
                  "Myofibroblast", #cluster 1
                  "gCap",          #cluster 2
                  "aCap",          #cluster 3
                  "EndoMT",        #cluster 4
                  "Fibroblast",    #cluster 5
                  "Fibroblast",    #cluster 6
                  "Fibroblast",    #cluster 7
                  "Art",           #cluster 8
                  "Vein",          #cluster 9
                  "SMC",           #cluster 10
                  "Fibroblast",    #cluster 11
                  "Fibroblast"     #cluster 12
)

names(cluster_ids2) <- levels(P14_integrated)
P14_integrated <- RenameIdents(P14_integrated, cluster_ids2)
P14_integrated$Celltype_fine <- Idents(P14_integrated)

Idents(P14_integrated) <- "seurat_clusters"
cluster_ids3 <- c("Fibroblast",    #cluster 0
                  "Myofibroblast", #cluster 1
                  "Endothelial",   #cluster 2
                  "Endothelial",   #cluster 3
                  "EndoMT",        #cluster 4
                  "Fibroblast",    #cluster 5
                  "Fibroblast",    #cluster 6
                  "Fibroblast",    #cluster 7
                  "Endothelial",   #cluster 8
                  "Endothelial",   #cluster 9
                  "SMC",           #cluster 10
                  "Fibroblast",    #cluster 11
                  "Fibroblast"     #cluster 12
)

names(cluster_ids3) <- levels(P14_integrated)
P14_integrated <- RenameIdents(P14_integrated, cluster_ids3)
P14_integrated$Celltype_main <- Idents(P14_integrated)

### Cell proportion ###
table(GSE151974_subset_P7_integrated1$Celltype_main, GSE151974_subset_P7_integrated1$Oxygen)
table(GSE151974_subset_P7_integrated1$Celltype_main, GSE151974_subset_P7_integrated1$Oxygen)




###################################
##### GEO accession:GSE211356 #####
###################################

### Data importing ###
fs = list.files(pattern = "GSM.") # Download: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE211356&format=file
samples = str_split(fs, "_", simplify = T)[, 1]

lapply(unique(samples), function(x){
  y = fs[grepl(x, fs)]
  folder = paste0(str_split(y[1], "_", simplify = T)[, 1])
  dir.create(folder, recursive = T)

  file.rename(paste0(y[1]), file.path(folder, "barcodes.tsv.gz"))
  file.rename(paste0(y[2]), file.path(folder, "features.tsv.gz"))
  file.rename(paste0(y[3]), file.path(folder, "matrix.mtx.gz"))
})
fs_name <- c("GSM6467309", "GSM6467310", "GSM6467311", "GSM6467312")
sce <- lapply(fs_name, function(x){
              a = Read10X(x)
              sce <- CreateSeuratObject(a)
              })


sce1 <- merge(x = sce[[1]],
              y = c(sce[[2]], sce[[3]], sce[[4]]),
              add.cell.ids = c("GSM6467309", "GSM6467310", "GSM6467311", "GSM6467312"),
              merge.data = TRUE)

### Quality control ###
source("./RunQC.R")
sce1 <- RunQC(sce1,
              org = "mus", 
              LowerFeatureCutoff = 200, 
              UpperFeatureCutoff = "MAD", 
              UpperMitoCutoff = 5, 
              Hb = F,
              doubletdetection = T, 
              decontXCutoff = 0.2,
              dir = "./")


orig_ident <- sapply(strsplit(rownames(sce1@meta.data), "_"), function(x) x[1])
sce1 <- AddMetaData(sce1, 
                    metadata = orig_ident,
                    col.name = "orig.ident")

sce1$Oxygen <- c(rep("Hyperoxia", 8633),
                 rep("Hyperoxia", 6926),
                 rep("Normoxia", 7979),
                 rep("Normoxia", 10566))

sce1$Sex <- c(rep("Female", 8633),
              rep("Male", 6926),
              rep("Female", 7979),
              rep("Male", 10566))

sce1@meta.data <- sce1@meta.data[-4]

### ProcessExper ###
sce_integrated = processExper(sce1, 
                              org = "mus",
                              cyclescoring = F,
                              sct.method = T,
                              reduction = "cca")

### Cell clustering and dimensional reduction ###
sce_integrated <- RunPCA(sce_integrated)
dims = 1:40
sce_integrated <- RunUMAP(sce_integrated,
                          dims = dims,
                          reduction.name = "umap") %>%
                  FindNeighbors(dims = dims) %>%
                  FindClusters(resolution = c(seq(0, 1, .1)))

clustree(sce_integrated) # Select suitable resolution

sce_integrated <- FindClusters(sce_integrated, resolution = 0.3)

### Cell annotation ###
marker = c("Epcam", # Epithelium
           "Pecam1", # Endothelium
           "Col1a1", # Mesenchyme
           "Ptprc" # Immune
           )

p <- DotPlot(sce_integrated, 
             assay = "SCT",
             features = marker
            ) +
     theme(axis.title = element_blank(),
           axis.line = element_blank(),
           axis.ticks.x = element_blank(),
           axis.text.x = element_text(size = 10),
           panel.background = element_rect(color = "black"),
           legend.position = "none") +
     coord_flip()

p1 <- DimPlot(sce_integrated, 
              reduction = "umap", 
              group.by = "seurat_clusters", 
              label = T)
wrap_plots(p + p1) # Reference


cluster_ids <- c("Endothelium", #cluster 0
                 "Epithelium",  #cluster 1
                 "Immune",      #cluster 2
                 "Epithelium",  #cluster 3
                 "Endothelium", #cluster 4
                 "Immune",      #cluster 5
                 "Immune",      #cluster 6
                 "Immune",      #cluster 7
                 "Immune",      #cluster 8
                 "Endothelium", #cluster 9
                 "Immune",      #cluster 10
                 "Mesenchyme",  #cluster 11
                 "Mesenchyme",  #cluster 12
                 "Endothelium", #cluster 13
                 "Immune",      #cluster 14
                 "Immune",      #cluster 15
                 "Mesenchyme",  #cluster 16
                 "Endothelium", #cluster 17
                 "Immune",      #cluster 18
                 "Immune",      #cluster 19
                 "Other",      #cluster 20
                 "Endothelium", #cluster 21
                 "Epithelium",  #cluster 22
                 "Immune",      #cluster 23
                 "Other"       #cluster 24
                )

names(cluster_ids) <- levels(sce_integrated)
sce_integrated <- RenameIdents(sce_integrated, cluster_ids)
sce_integrated$Celltype_main <- Idents(sce_integrated)

### ECs and Mesenchyme seperation ###
Idents(sce_integrated) <- "Celltype_main"
subset <- subset(sce_integrated, idents = c("Endothelium", "Mesenchyme"))

subset1 <- CreateSeuratObject(counts = subset@assays$RNA@counts,
                              meta.data = subset@meta.data)


###ProcessExper
subset1 <- processExper(subset1, 
                        org = "mus",
                        cyclescoring = F,
                        sct.method = T,
                        reduction = "cca")


### Cell clustering and dimensional reduction ###
subset1 <- RunPCA(subset1)
dims = 1:40
subset1 <- RunUMAP(subset1,
                   dims = dims,
                   reduction.name = "umap") %>%
           FindNeighbors(dims = dims) %>%
           FindClusters(resolution = c(seq(0, 1, .1)))

clustree(subset1) #Find suitable resolution value

subset1 <- FindClusters(subset1, resolution = 0.7)


marker1 = c("Gpihbp1", "Kit", #gCap
            "Car4", "Kdr", #aCap
            "Cxcl12", "Pcsk5", #Art
            "Vegfc", "Prss23", #Vein
            "Ccl21a", "Mmrn1", # Lymph
            "Pecam1", "Eng", "Cd34", "Cdh5", #Gen Endothelium
            "Col1a1", "Col1a2", "Col3a1", "Fn1", # Fibroblast
            "Tagln", "Acta2", "Myl9", "Myh11", # SMC
            "Tgfbi", "Wnt5a" # Myofibroblast
            )

p2 <- DotPlot(subset1, 
              assay = "SCT",
              features = marker1
             ) +
      theme(axis.title = element_blank(),
            axis.line = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(size = 10),
            panel.background = element_rect(color = "black")
           ) +
      coord_flip()

p3 <- DimPlot(subset1, 
              reduction = "umap", 
              group.by = "seurat_clusters", 
              label = T)
wrap_plots(p2 + p3) # Reference

cluster_ids1 <- c("gCap",          #cluster 0
                  "aCap",          #cluster 1
                  "Art",           #cluster 2
                  "gCap",          #cluster 3
                  "Myofibroblast", #cluster 4
                  "gCap",          #cluster 5
                  "gCap",          #cluster 6
                  "gCap",          #cluster 7
                  "Art",           #cluster 8
                  "gCap",          #cluster 9
                  "aCap",          #cluster 10
                  "Fibroblast",    #cluster 11
                  "Vein",          #cluster 12
                  "Fibroblast",    #cluster 13
                  "Fibroblast",    #cluster 14
                  "aCap",          #cluster 15
                  "EndoMT",        #cluster 16
                  "SMC",           #cluster 17
                  "Lymph",         #cluster 18
                  "aCap",          #cluster 19
                  "Fibroblast",    #cluster 20
                  "Fibroblast"     #cluster 21
                 )

names(cluster_ids) <- levels(subset1)
subset1 <- RenameIdents(subset1, cluster_ids)
subset1$Celltype_fine <- Idents(subset1)

### Vlnplot of EndoMT markers ###
Idents(subset1) <- "Celltype_fine"
subset2 <- subset(subset1, idents = c("gCap", "aCap", "Art", "Vein", "EndoMT", "Fibroblast", "Myofibroblast", "SMC"))
                  
levels = c("gCap",
           "aCap",
           "Art",
           "Vein",
           "EndoMT",
           "Fibroblast",
           "Myofibroblast",
           "SMC")
                     
subset2$Celltype_fine <- factor(subset2$Celltype_fine, levels = levels)

Idents(subset2) <- "Celltype_fine"

Split_Vln_stacked(subset2,
                  assay = "SCT",
                  feature = c("Acta2", "Myl9", "Tagln",
                              "Cdh5", "Eng", "Pecam1"),
                  pt.size = 0, 
                  text_x_size = 12,
                  text_y_size = 10,
                  title_y_size = 14,
                  cols = c("#ea5c6f", "#f7905a", "#e187cb", "#fb948d", "#e2b159", "#ebed6f", "#b2db87", "#7ee7bb")
                 )
