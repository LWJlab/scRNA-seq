###################################
##### GEO accession:GSE151974 #####
###################################
library(Seurat)
library(SeuratObject)
library(clustree)
library(patchwork)
options(stringsAsFactors=FALSE)

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
            group.by = 'Celltype',
            features = c('Gpihbp1', 'Kit', # gCap
                         'Car4', 'Kdr', # aCap
                         'Cxcl12', 'Pcsk5', # Art
                         'Vegfc', 'Prss23', # Vein
                         'Pecam1', 'Eng', 'Cd34', 'Cdh5', # Gen ECs
                         'Col1a1', 'Col1a2', 'Col3a1', 'Fn1', 'Tagln', 'Acta2', 'Myl9', 'Myh11', # Mesenchyme
                         'Tgfbi','Wnt5a' #Myofibroblast
                         ) 
            ) +
    theme(axis.title =element_blank(),
          axis.line = element_blank(),
          axis.ticks.x =element_blank(),
          axis.text.x = element_text(size=10),
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
           'Pecam1','Eng', 'Cd34', 'Cdh5', # General Endothelium
           'Col1a1', 'Col1a2', 'Col3a1', 'Fn1', 'Tagln', 'Acta2', 'Myl9','Myh11', # Mesenchyme
           'Tgfbi','Wnt5a' # Myofibroblast
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
                                     'Myfibroblast',
                                     'SMC'
                                    ), 
                          )

P7_dotplot


