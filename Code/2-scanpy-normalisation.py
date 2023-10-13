# load python libraries
import scanpy as sc
import anndata
import pandas as pd
import feather

# specify path to h5ad file exported from R
path_to_h5ad_file = "rds_file.h5ad"

# specify path to output feather file
# will be uploaded to R later
path_to_scanpy_output_feather_file = "scanpy_normalise_matrix.feather"

# read h5ad file exported from R
rds_file = anndata.read_h5ad(path_to_h5ad_file)

# run scanpy normalisation algorithm
sc.pp.normalize_total(rds_file, inplace=True)

# save normalised data to feather file
scanpy_normalise_matrix = rds_file.to_df()
scanpy_normalise_matrix = pd.DataFrame(scanpy_normalise_matrix)
scanpy_normalise_matrix = scanpy_normalise_matrix.transpose()
feather.write_dataframe(scanpy_normalise_matrix, path_to_scanpy_output_feather_file)

