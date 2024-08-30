###################################
##### GEO accession:GSE151974 #####
###################################
library(Seurat)
library(SeuratObject)
library(clustree)
library(patchwork)
options(stringsAsFactors = FALSE)

### Data importing ###
raw_count <- read.csv('GSE151974_raw_umi_matrix_postfilter.csv', row.names = 1) # https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151974/suppl/GSE151974%5Fcell%5Fmetadata%5Fpostfilter.csv.gz
metadata <- read.csv('GSE151974_cell_metadata_postfilter.csv', row.names = 1) # https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151974/suppl/GSE151974%5Fraw%5Fumi%5Fmatrix%5Fpostfilter.csv.gz

sce <- CreateSeuratObject(counts = raw_count,
                          meta.data = metadata)


### P7 seperation ###
Idents(sce) <- 'CellType'
unique(sce$CellType)
subset <- subset(sce, idents = c('Cap-a', 'Art', 'Cap', 'Vein', 'Col13a1+ fibroblast', 'Col14a1+ fibroblast', 'Myofibroblast','SMC'))
Idents(subset) <- 'Sample'
P7_subset <- subset(subset, idents = c('P7_Hyperoxia', 'P7_Normoxia'))


### ProcessExper ###
source('./processExper.R')
P7_integrated = processExper(P7_subset, 
                             sct.method = T, 
                             reduction = 'cca')


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
             assay = 'SCT',
             group.by = 'seurat_clusters',
             features = c('Gpihbp1', 'Kit', # gCap
                          'Car4', 'Kdr', # aCap
                          'Cxcl12', 'Pcsk5', # Art
                          'Vegfc', 'Prss23', # Vein
                          'Pecam1', 'Eng', 'Cd34', 'Cdh5', # Gen ECs
                          'Col1a1', 'Col1a2', 'Col3a1', 'Fn1', 'Tagln', 'Acta2', 'Myl9', 'Myh11', # Mesenchyme
                          'Tgfbi','Wnt5a' #Myofibroblast
                          ) 
             ) +
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10),
        panel.background = element_rect(color = 'black'))+
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
                 'Myofibroblast', #cluster 4
                 'Myofibroblast', #cluster 5
                 "Fibroblast",    #cluster 6
                 'EndoMT',        #cluster 7
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
P7_integrated$Celltype <- Idents(P7_integrated)

Idents(P7_integrated) <- 'seurat_clusters'
cluster_ids1 <- c("Endothelium",   #cluster 0
                  "Fibroblast",    #cluster 1
                  "Fibroblast",    #cluster 2
                  "Endothelium",   #cluster 3
                  'Myofibroblast', #cluster 4
                  'Myofibroblast', #cluster 5
                  "Fibroblast",    #cluster 6
                  'EndoMT',        #cluster 7
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
source('./cellatlas_umap.R')
atlas <- cellatlas_umap(P7_integrated, 
                        idents = 'Celltype',
                        levels = c('gCap',
                                   'aCap',
                                   'Art',
                                   'Vein',
                                   'EndoMT',
                                   'Fibroblast',
                                   'Myofibroblast',
                                   'SMC'), 
                        hull_alpha = 0.1, 
                        hull_size = 0.5, 
                        hull_lty = 2,
                        hull_delta = 0.8,
                        dot_color = c("#ea5c6f","#f7905a","#e187cb","#fb948d","#e2b159","#ebed6f","#b2db87","#7ee7bb"),  
                        dot_size = 0.5,
                        dot_alpha = 1, 
                        label_size = 4, 
                        label_color = F
                       )
atlas

oxygen_atlas <- cellatlas_umap(P7_integrated, 
                               idents = 'Oxygen',
                               levels = c('Normoxia','Hyperoxia'),
                               hull_alpha = 0,
                               hull_size = 0, 
                               hull_lty = 2, 
                               hull_delta = 0, 
                               dot_color = c("#4DBBD5FF",'#E64B35FF'),  
                               dot_size = 0.5, 
                               dot_alpha = 1, 
                               label_size = 4, 
                               label_color = T
                              )
oxygen_atlas


### Dotplot of cell markers ###
source('./sce_dotplot.R')
marker = c('Gpihbp1', 'Kit', # gCap
           'Car4', 'Kdr',  # aCap
           'Cxcl12', 'Pcsk5', # Art
           'Vegfc', 'Prss23', # Vein
           'Pecam1', 'Eng', 'Cd34', 'Cdh5', # General Endothelium
           'Col1a1', 'Col1a2', 'Col3a1', 'Fn1', 'Tagln', 'Acta2', 'Myl9', 'Myh11', # Mesenchyme
           'Tgfbi', 'Wnt5a' # Myofibroblast
          )

P7_dotplot <- sce_dotplot(P7_integrated,
                          assay = 'SCT',
                          idents = 'Celltype',
                          markers = marker,
                          levels = c('gCap',
                                     'aCap',
                                     'Art',
                                     'Vein',
                                     'EndoMT',
                                     'Fibroblast',
                                     'Myofibroblast',
                                     'SMC'
                                    ) 
                         )

P7_dotplot


### Dotplot of metabolism-related genes ###
marker1 = c('Adh5', 'Aldh3a2', 'Gapdh', 'Gpi1', 'Ldha', 'Pkm', 'Slc2a1', # Glycolysis
            'G6pdx', 'Pgd', 'Taldo1', 'Tkt', # Pentose phosphate pathway
            'Dhfr', 'Mthfd1', 'Mthfd2', 'Shmt1', 'Shmt2', 'Slc25a32', # One carbon metabolism
            'Acat1', 'Aldh7a1', 'Gcdh',  # Amino acid metabolism
            'Cpt1a', 'Cpt1c', 'Cpt2', # Carnitine shuttle
            'Acaa2', 'Acadm', 'Acads', 'Acadsb', 'Acadvl', 'Acox3', 'Echs1', 'Eci1', 'Eci2', 'Hadh', 'Hadha','Hadhb', # β-oxidation
            'Acsl1', 'Acsl4', 'Acaca', 'Degs1', # Fatty acid synthesis
            'Apoe', 'Lipa', 'Lipe', 'Lpl', 'Mttp', 'Pnpla2', # Lipid metabolism and transport
            'Cyc1', 'Ndufs1', 'Ndufs2', 'Ndufs3', 'Pdha1', 'Sdha' # Oxidative phosphorylation
           )

P7_Ec_dotplot <- sce_dotplot(subset(P7_integrated, idents = c('gCap','aCap','Art','Vein','EndoMT')),
                             assay = 'SCT',
                             idents = 'Oxygen',
                             markers = marker1,
                             levels = c('Normoxia', 'Hyperoxia'),
                             title = 'Endothelium',
                             title_size = 10
                            )

P7_EndoMT_dotplot <- sce_dotplot(subset(P7_integrated, idents = c('EndoMT')),
                                 assay = 'SCT',
                                 idents = 'Oxygen',
                                 markers = marker1,
                                 levels = c('Normoxia', 'Hyperoxia'),
                                 title = 'EndoMT',
                                 title_size = 10
                                )

wrap_plots(P7_Ec_dotplot + P7_EndoMT_dotplot)


### Vlnplot of EndoMT markers ###
levels = c('gCap',
           'aCap',
           'Art',
           'Vein',
           'EndoMT',
           'Fibroblast',
           'Myofibroblast',
           'SMC')
P7_integrated$Celltype <- factor(P7_integrated$Celltype, levels = levels) 
Idents(P7_integrated) <- 'Celltype'

Split_Vln_stacked(P7_integrated,
                  assay = 'SCT',
                  feature = c('Acta2','Myl9','Tagln',
                              'Cdh5','Eng','Pecam1'),
                  split.plot = F,
                  pt.size = 0, 
                  size = 10,
                  cols = c("#ea5c6f","#f7905a","#e187cb","#fb948d","#e2b159","#ebed6f","#b2db87","#7ee7bb"),
                  test = F)


### Trajectory analysis ###
library(slingshot)

P7_integrated1 <- as.SingleCellExperiment(P7_integrated, assay = 'SCT')

sds <- slingshot(P7_integrated1, 
                 clusterLabels = "Celltype_main",
                 reducedDim = 'UMAP' ,
                 start.clus = "Endothelium"
                )

source('./pseudotime_umap.R')
pseu1 <- pseudotime_umap(P7_integrated,
                         pseudotime_obj = sds,
                         idents = 'Celltype',  
                         sub_idents = 'slingPseudotime_1', 
                         dot_color = c("#440154FF", "#FDE725FF"),
                         dot_size = 0.5, 
                         dot_alpha = 1, 
                         label_size = 4
                        )
pseu1

pseu2 <- pseudotime_umap(P7_integrated,
                         pseudotime_obj = sds,
                         idents = 'Celltype',  
                         sub_idents = 'slingPseudotime_2', 
                         dot_color = c("#440154FF", "#FDE725FF"),
                         dot_size = 0.5, 
                         dot_alpha = 1, 
                         label_size = 4
                        )
pseu2


### GSEA analysis ###
library(msigdbr)
library(tibble)
library(patchwork)

c2_mmu <- msigdbr(species = "Mus musculus", category = c("C2")) %>% as.data.frame %>% dplyr::select(gs_name,entrez_gene,gene_symbol) %>% as.data.frame
c5_mmu <- msigdbr(species = "Mus musculus", category = c("C5")) %>% as.data.frame %>% dplyr::select(gs_name,entrez_gene,gene_symbol) %>% as.data.frame

c2_KEGG_mmu_pathway <- c2_mmu[grep(c('KEGG'), c2_mmu$gs_name),]
c2_PID_mmu_pathway <- c2_mmu[grep(c('PID_'), c2_mmu$gs_name),]
c2_BIOCARTA_mmu_pathway <- c2_mmu[grep(c('BIOCARTA'), c2_mmu$gs_name),]
c2_REACTOME_mmu_pathway <- c2_mmu[grep(c('REACTOME'), c2_mmu$gs_name),]
c2_WIKIPATHWAYS_mmu_pathway <- c2_mmu[grep(c('WIKIPATHWAYS'), c2_mmu$gs_name),]
c5_GOBP_mmu_pathway <- c5_mmu[grep(c('GOBP'), c5_mmu$gs_name),]

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

source("./FindMarker_genes.R")
P7_integrated <- PrepSCTFindMarkers(P7_integrated)
DEGs <- FindMarker_genes(dataset = P7_integrated, 
                         assay = 'SCT',
                         clusters = c('gCap','aCap','Art','Vein','EndoMT','Fibroblast','SMC','Myofibroblast'),
                         comparison = c("Oxygen", "Normoxia", "Hyperoxia"),
                         logfc.threshold = 0,  
                         min.cells.group = 1
                        )  

source("./sce_GSEA.R")
gsea_res <- sce_GSEA(DEGs, pathway = mmu_fgsea_sets)

gsea_res1 <- gsea_res
gsea_res1 <- do.call(rbind, gsea_res1) %>% t()
#write.csv(gsea_res1, file = "./P7_MSigDB_C2CP_GSEA_results.csv", row.names = F)


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
                       maxPathway = c('Carbohydrate catabolic process',
                                      'ADP metabolic process',
                                      'Membrane lipid metabolic process',
                                      'Positive regulation of Erk1 and Erk2 cascade',
                                      'Cytokine mediated signaling pathway',
                                      'p53 signaling pathway',
                                      'Interleukin-4 and interleukin-13 signaling'
                                     ),
                       category_color = 'grey10',
                       pathway_color = c("Carbohydrate metabolism" = "#66C2A5", 
                                         "Energy metabolism" = "#8DA0CB", 
                                         "Lipid metabolism" = "#A6D854", 
                                         "Signal transduction" = "#6a3d9a",
                                         "Signaling molecules and interaction" = '#E5C494',
                                         "Cell growth and death" = '#FC8D62',
                                         "Immune system" = '#E78AC3'),
                       title = 'Hyperoxia vs. Normoxia (EndoMT)',
                       title_size = rel(1),
                       num_size = 2.5,
                       text_x_size = rel(1),
                       text_y_size = rel(1),
                       xlim = c(0, 2.3))

bar
