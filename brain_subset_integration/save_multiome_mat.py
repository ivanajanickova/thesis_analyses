print("Loading libraries")
import os
from pycisTopic.cistopic_class import *
from pycisTopic.qc import *
from pycisTopic.lda_models import *
from pycisTopic.clust_vis import *
from pycisTopic.diff_features import *
from pycisTopic.pseudobulk_peak_calling import *
from pycisTopic.iterative_peak_calling import *
from pycisTopic.topic_binarization import *
from pycisTopic.gene_activity import *
from loomxpy.loomxpy import SCopeLoom
import regex as re
import pickle
import sys
import pyarrow.feather as feather
print(sys.version)

print("Loading atac")
with open("/staging/leuven/stg_00002/lcb/ijanic/projects/subclustering/atlas_with_brain_annots.pkl", "rb") as f:
    atlas_obj = pickle.load(f)  
print("FCA-atac atlas loaded")

mo_cells = atlas_obj.cell_data.loc[atlas_obj.cell_data.Sample == "FDM__2cf092__10x_Multiome_cDNA_lib_fly_head_1",:].index
print(f"We have {len(mo_cells)} multiome head cells.")
mo_cells_set = {cell for cell in mo_cells}

mo_df = atlas_obj.cell_data.loc[mo_cells,:]


print("Calculating data matrix")
data_mat =atlas_obj.selected_model.cell_topic_harmony
# Scale the data
data_mat = pd.DataFrame(sklearn.preprocessing.StandardScaler().fit_transform(
            data_mat), index=data_mat.index.to_list(), columns=data_mat.columns)


mo_df = pd.DataFrame()

print("Selecting Brain samples")
for col_name in data_mat.columns:
    if col_name in mo_cells_set:
        mo_df[col_name] = data_mat[[col_name]]
mo_df['rownames'] = mo_df.index


print("saving")
feather.write_feather(mo_df,  "/scratch/leuven/338/vsc33893/bridge/inputs/multiome_atac_v4.feather")
print("saved")