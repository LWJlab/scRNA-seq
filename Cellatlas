###################################
############ GSE151974 ############
###################################
library(Seurat)
library(SeuratObject)

options(stringsAsFactors=FALSE)

### Data loading ###
GSE151974_raw_count <- read.csv('GSE151974_raw_umi_matrix_postfilter.csv', row.names = 1)
dim(GSE151974_raw_count)

GSE151974_metadata <- read.csv('GSE151974_cell_metadata_postfilter.csv', row.names = 1)
head(GSE151974_metadata)
dim(GSE151974_metadata)
table(GSE151974_metadata$CellType, GSE151974_metadata$Oxygen)

GSE151974_seurat <- CreateSeuratObject(counts = GSE151974_raw_count,
                                       meta.data = GSE151974_metadata)

head(GSE151974_seurat)

Idents(GSE151974_seurat) <- 'CellType'

unique(GSE151974_seurat$CellType)

GSE151974_subset <- subset(GSE151974_seurat, idents = c('Cap-a', 'Art', 'Cap', 'Vein', 
                                                        'Col13a1+ fibroblast', 'Col14a1+ fibroblast', 'Myofibroblast',
                                                        'SMC'))
GSE151974_seurat

###### P7 ######
Idents(GSE151974_subset) <- 'Sample'
GSE151974_subset_P7 <- subset(GSE151974_subset, idents = c('P7_Hyperoxia', 'P7_Normoxia'))
unique(GSE151974_subset_P7$orig.ident)
Idents(GSE151974_subset_P7)

###ProcessExper
setwd('/data2/fanjie/Yao_H/')
source('/data2/fanjie/Yao_H/processExper.R')

GSE151974_subset_P7_integrated = processExper(GSE151974_subset_P7, 
                                              org = 'mus',
                                              cyclescoring = F,
                                              sct.method = T,
                                              reduction = 'cca')
