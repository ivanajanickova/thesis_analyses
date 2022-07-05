message("Loading Libraries")
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
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



##########################################################
##### Preparing RNA object ###############################
##########################################################

message("Loading RNA head objects")
Convert("/staging/leuven/stg_00002/lcb/ijanic/projects/FCA_atac/integration/matrices/s_fca_biohub_head_10x.h5ad", "/staging/leuven/stg_00002/lcb/ijanic/projects/FCA_atac/integration/matrices/rna_head.h5Seurat")
rna_obj <- LoadH5Seurat("/staging/leuven/stg_00002/lcb/ijanic/projects/FCA_atac/integration/matrices/rna_head.h5Seurat")
rna_loom <-  open_loom("/staging/leuven/stg_00002/lcb/ijanic/projects/FCA_atac/integration/matrices/rna_head_loom.loom", mode = "r")

clustering = get_clusterings_with_name(rna_loom)
rna_obj = RenameCells(rna_obj, rownames(clustering))
clustering_columns = c('-1', '-2', '-3', '-4', '-5')
for(clust in clustering_columns){
  name = paste0("annot_", strsplit(clust, '-')[[1]][2])
  rna_obj = AddMetaData(rna_obj, clustering$clust)
}
rna_obj = AddMetaData(rna_obj, clustering$`-1`, col.name = "annot_1")
Idents(rna_obj) = clustering$`-1`
rna_obj@assays$SCT = rna_obj@assays$RNA

######## SAVING #########
rna_obj = SCTransform(object = rna_obj) %>% RunPCA() %>% RunUMAP(dims = 1:50, return.model = TRUE) 
saveRDS(rna_obj, file = "/staging/leuven/stg_00002/lcb/ijanic/projects/FCA_atac/integration/matrices/rna_annotated.rds")



##########################################################
##### Preparing Multiome object ##########################
##########################################################

message("Loading Head Multiome Data")
file_path = "/staging/leuven/stg_00002/lcb/lcb_projects/FDM/multiome/NextSeq2000_20210719__NextSeq2000_20210720__NextSeq2000_20210723/FDM__2cf092__10x_Multiome_cDNA_lib_fly_head_1/outs/filtered_feature_bc_matrix.h5"

input.1 = Read10X_h5(file_path)


# change the colname formatting so it matches ATAC
change_colnames = function(count_mat, sample_name) {
  new_colnames = c()
  for(colname in colnames(count_mat)){
    new = paste0(colname, sample_name)
    new_colnames = c(new_colnames, new)
  }
  colnames(count_mat) = new_colnames
  count_mat
}


rna_counts1 = input.1$`Gene Expression`

### Adding sample string to the RNA barcodes
rna_counts1 = change_colnames(rna_counts1, "___FDM__2cf092__10x_Multiome_cDNA_lib_fly_head_1")
message("RNA counts format")
print(head(rna_counts1))


atac_counts <- read_feather("/scratch/leuven/338/vsc33893/bridge/inputs/multiome_atac_v4.feather")
print("ATAC counts format")
print(head(atac_counts))
### Filter non-intersecting cells
# Keep only cell present in both atac_counts and rna_counts dataframe
atac_colnames = colnames(atac_counts)
rna_colnames = colnames(rna_counts1)

common_cells = intersect(atac_colnames, rna_colnames)
message("The lenght of the common cells is....")
print(length(common_cells))

atac_counts = atac_counts %>% select(common_cells)
rna_counts1 = as.data.frame(rna_counts1) %>% select(intersect(common_cells, colnames(rna_counts1)))

### Create Seurat objs from the RNA count matrices and merge the objects
# Create Seurat object
obj.multi = CreateSeuratObject(counts = rna_counts1)

# Get % of mitochondrial genes
obj.multi[["percent.mt"]] <- PercentageFeatureSet(obj.multi, pattern = "^MT-")



### Add ATAC assay to the multiome object


atac_assay = CreateAssayObject(counts = atac_counts)
obj.multi[["ATAC"]] <- atac_assay


DefaultAssay(obj.multi) <- "RNA"
obj.multi = NormalizeData(obj.multi, normalization.method = "LogNormalize")
`
# normalize multiome ATAC
DefaultAssay(obj.multi) <- "ATAC"
obj.multi <- FindTopFeatures(obj.multi, min.cutoff = "q0")

print("Saving")
saveRDS(obj.multi, "/staging/leuven/stg_00002/lcb/ijanic/projects/FCA_atac/integration/matrices/mutiome_obj_v4.rds")
print("Saved")