message("Loading libraries")
library(Signac)
library(dplyr)
library(ggplot2)
library(SeuratDisk)
library(biovizBase)
library(ensembldb)
library(data.table)
library(Matrix)
library(genomation)
library(feather)
library(SCopeLoomR)
library(devtools)
load_all("/staging/leuven/stg_00002/lcb/ijanic/software/seurat/")

message("Loading objects")
obj.multi = readRDS("/staging/leuven/stg_00002/lcb/ijanic/projects/FCA_atac/integration/matrices/mutiome_obj_v4.rds")
obj.rna <- readRDS("/scratch/leuven/338/vsc33893/bridge/objects/rna.rds")

#Create ATAC
path <- "/scratch/leuven/338/vsc33893/bridge/inputs/brain_counts.feather"
atac_query_counts = read_feather(path)
obj.atac = CreateSeuratObject(counts = atac_query_counts, assay="ATAC")


# Drop first dimension for ATAC reduction
dims.atac <- 2:449
dims.rna <- 1:61
DefaultAssay(obj.multi) <-  "RNA"
DefaultAssay(obj.rna) <- "RNA"
obj.rna.ext <- PrepareBridgeReference(reference = obj.rna,
                                      bridge = obj.multi, 
                                      reference.reduction = "pca",
                                      reference.dims = dims.rna,
                                      normalization.method = "LogNormalize",
                                      slsi.reduction.dims = 450,
                                      transfer.anchors.k.anchor = 5,
                                      transfer.anchors.n.trees = 50
)


bridge.anchor <- FindBridgeTransferAnchors(extended.reference = obj.rna.ext, 
                                           query = obj.atac,
                                           reduction = "lsiproject",
                                           dims = dims.atac
)

obj.atac <- MapQuery(anchorset = bridge.anchor, 
                     reference = obj.rna, 
                     query = obj.atac, 
                     refdata = list(
                       l1 = "annot_1"),
                     reduction.model = "umap"
)

message("Saving the matrix")
path <- "/scratch/leuven/338/vsc33893/bridge/results/brain_optimized_slsi.feather"
obj.atac@meta.data$rownames = row.names(obj.atac@meta.data)
write_feather(obj.atac@meta.data, path)