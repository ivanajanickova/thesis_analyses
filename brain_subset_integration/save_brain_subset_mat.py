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

brain_cells = [atlas_obj.cell_data.index[i] for i in range(0, len(atlas_obj.cell_data)) if type(atlas_obj.cell_data["brain_annot"][i]) == str]
print(f"We have {len(brain_cells)} cells from the fly brain project.")

brain_df = atlas_obj.cell_data[~atlas_obj.cell_data["brain_annot"].isnull()]

brain_cells_set = set()
for i in range(0, len(brain_df)):
    brain_cells_set.add(brain_df.iloc[i, :].name)
        
print(f"Count of the selected cells {len(brain_cells_set)}")


print("Calculating data matrix")
data_mat =atlas_obj.selected_model.cell_topic_harmony
# Scale the data
data_mat = pd.DataFrame(sklearn.preprocessing.StandardScaler().fit_transform(
            data_mat), index=data_mat.index.to_list(), columns=data_mat.columns)


brain_df = pd.DataFrame()

print("Selecting Brain samples")
for col_name in data_mat.columns:
    if col_name in brain_cells_set:
        brain_df[col_name] = data_mat[[col_name]]
brain_df['rownames'] = brain_df.index

annot_df = pd.DataFrame()
annot_df['brain_annot'] = atlas_obj.cell_data.loc[brain_cells_set, "brain_annot"]
annot_df['cells'] =  atlas_obj.cell_data.loc[brain_cells_set, :].index

print("Filtering brain samples")

print("saving")
feather.write_feather(brain_df,  "/scratch/leuven/338/vsc33893/bridge/inputs/brain_counts.feather")
feather.write_feather(annot_df,  "/scratch/leuven/338/vsc33893/bridge/inputs/brain_annots.feather")
print("saved")
