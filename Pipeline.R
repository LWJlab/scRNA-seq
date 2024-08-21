###################################
############ GSE151974 ############
###################################
library(Seurat)
library(SeuratObject)

options(stringsAsFactors=FALSE)

### Data importing ###
raw_count <- read.csv('GSE151974_raw_umi_matrix_postfilter.csv', row.names = 1)
metadata <- read.csv('GSE151974_cell_metadata_postfilter.csv', row.names = 1)

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
GSE151974_subset_P7_integrated = processExper(P7_subset, 
                                              org = 'mus',
                                              cyclescoring = F,
                                              sct.method = T,
                                              reduction = 'cca')
