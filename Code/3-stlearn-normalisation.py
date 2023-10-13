# load python libraries
import stlearn as st
from pathlib import Path
import csv
import pandas as pd
import feather

# specify path to usual 'outs' folder, the result of 10X Genomics Space Ranger
path_to_visium_outs_folder = "stlearn-data/outs"

# specify path to spot tiles (the intermediate result of image pre-processing)
path_to_tiles_folder = "stlearn-data/tmp/tiles"

# specify the intermediate output path
path_to_output_folder = "stlearn-data/results"

# specify path to csv files exported from R containing gene and spot names
path_to_gene_names = "names_of_genes.csv"
path_to_spot_names = "names_of_spots.csv"

# specify path to output feather file
# will be uploaded to R later
path_to_stlearn_output_feather_file = "stlearn_normalise_matrix.feather"

BASE_PATH = Path(path_to_visium_outs_folder)
TILE_PATH = Path(path_to_tiles_folder)
TILE_PATH.mkdir(parents=True, exist_ok=True)
OUT_PATH = Path(path_to_output_folder)
OUT_PATH.mkdir(parents=True, exist_ok=True)

# find genes that weren't discarded during filtering
gene_list = []
with open(path_to_gene_names, newline='') as genes_file:
    reader = csv.reader(genes_file, delimiter=',')
    for row in reader:
        gene_list.append(row[0])
gene_list = gene_list[1:]

# find spots that weren't discarded during filtering
spot_list = []
with open(path_to_spot_names, newline='') as spots_file:
    reader2 = csv.reader(spots_file, delimiter=',')
    for row in reader2:
        spot_list.append(row[0])
spot_list = spot_list[1:]

# load Visium data
data = st.Read10X(BASE_PATH, count_file='raw_feature_bc_matrix.h5')
data.layers["raw"] = data.X

# add columns for filtering genes and spots
data.var['keep_gene'] = data.var['gene_ids'].isin(gene_list)
data.obs['keep_spots'] = data.obs.index.isin(spot_list)

# get filtered data
filtered_data = data[data.obs['keep_spots'], data.var['keep_gene']].copy()

# run first normalisation step
st.pp.normalize_total(filtered_data)
st.pp.log1p(filtered_data)

# run pre-processing for spot image
st.pp.tiling(filtered_data, TILE_PATH)

# extract high-level features from tile images
st.pp.extract_feature(filtered_data)

# runpPrincipal component analysis for gene expression data
st.em.run_pca(filtered_data,n_comps=50)

# run second normalisation step
data_SME = filtered_data.copy()
st.spatial.SME.SME_normalize(data_SME, use_data="raw")
data_SME.X = data_SME.obsm['raw_SME_normalized']

# save normalised data to feather file
stlearn_normalise_matrix = pd.DataFrame(data_SME.X)
stlearn_normalise_matrix = stlearn_normalise_matrix.transpose()
feather.write_dataframe(stlearn_normalise_matrix, path_to_stlearn_output_feather_file)

