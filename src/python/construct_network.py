import pandas as pd
import numpy as np
from tqdm import tqdm

group = "L"
ugrn_path = "data/unweighted_network/5k/ugrn.csv"
chrom_rp_path = f"data/tmp/outputs/{group}/chrom_rp.csv"
peak_rp_path = f"data/tmp/outputs/{group}/peak_rp.csv"
output_file = f"data/tmp/outputs/{group}/network.csv"

ugrn = pd.read_csv(ugrn_path, names=["tf", "gene"])
# print(ugrn.shape)  # Edgelist: (10229161, 2)

# tf_list = ugrn["tf"].unique()
# gene_list = ugrn["gene"].unique()

# Import TF_RP
peak_rp = pd.read_csv(peak_rp_path, index_col=0)
# peak_rp = peak_rp.loc[peak_rp.index.isin(gene_list), peak_rp.columns.isin(tf_list)]
# print(peak_rp.shape)  # Matrix: (18833 genes, 1713 TFs)

# Import H3K27ac_RP
chrom_rp = pd.read_csv(chrom_rp_path, index_col=0, names=["h3k27ac"])
chrom_rp = chrom_rp.loc[peak_rp.index.to_list()]  # Ensure indices are the same
# print(chrom_rp.shape)  # Column vector: (18833 genes, 1)


# Compute weights
def normalize(peak_rp):
    # Step 1: Log-transform the raw RP values
    log_peak_rp = np.log(peak_rp + 1)  # Add 1 to avoid log(0)

    # Step 2: Calculate the mean of the log-transformed RP values for each column (TF)
    mean_log_peak_rp = log_peak_rp.mean(axis=0)

    # Step 3: Subtract the mean log-transformed RP value of each column from the log-transformed RP values in that column
    normalized_log_peak_rp = log_peak_rp.subtract(mean_log_peak_rp)

    # Step 4: Exponentiate the result to obtain the normalized RP values
    # normalized_peak_rp = (
    #     np.exp(normalized_log_peak_rp) - 1
    # )  # Subtract 1 to revert the +1 added during log-transformation
    return normalized_log_peak_rp


normalized_peak_rp = normalize(peak_rp)
normalized_peak_rp = np.maximum(normalized_peak_rp, 0)
weights = normalized_peak_rp.to_numpy() * chrom_rp.to_numpy()
weights = pd.DataFrame(weights, columns=peak_rp.columns, index=peak_rp.index)
# print(weights.shape)  # Matrix: (18833 genes, 1713 TFs)

# Update weights
# Convert the weights DataFrame to a MultiIndex Series
weights_series = weights.stack()
weights_series.index.names = ["gene", "tf"]

# Set the index of ugrn to be a MultiIndex with 'gene' and 'tf'
ugrn = ugrn.set_index(["gene", "tf"])

# Use the `update` method to update the 'weight' column using the values from the weights_series
ugrn["weight"] = 0  # Initialize the 'weight' column
ugrn["weight"].update(weights_series)

# Reset the index
ugrn = ugrn.reset_index()

ugrn = ugrn[["tf", "gene", "weight"]]

ugrn.to_csv(output_file, index=False)
