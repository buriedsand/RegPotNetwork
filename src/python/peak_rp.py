import sys
import numpy as np
import pandas as pd
from tqdm import tqdm

input_file = sys.argv[1]
output_file = sys.argv[2]
refseq_file = sys.argv[3]

L = 100_000  # Size of genomic inverval on either side of the TSS
delta = 10_000  # Decay distance


def calculate_single_gene_rp(sample_values, t_k, L=100000, mu=1.0):
    R_jk = 0
    for row in sample_values.itertuples():
        # Create a NumPy array of nucleotide positions
        i_values = np.arange(getattr(row, "start"), getattr(row, "end"))

        # Compute d, w_i, and s_ij values
        d_values = np.abs(i_values - t_k) / L
        w_i_values = 2 * np.exp(-mu * d_values) / (1 + np.exp(-mu * d_values))

        s_ij_values = getattr(row, "s_ij")

        # Compute R_jk value
        R_jk += np.sum(w_i_values * s_ij_values)

    return R_jk


def calculate_peak_rp(binary_values, refseq_genes, L=100000, mu=1.0):
    peak_rp = dict()

    for gene in tqdm(refseq_genes.itertuples(), total=refseq_genes.shape[0]):
        gene_id = getattr(gene, "gene_id")
        gene_chr = getattr(gene, "chr")
        t_k = getattr(gene, "tss")

        interval_start = max(t_k - L, 0)
        interval_end = t_k + L
        interval_values = binary_values[
            (binary_values["chr"] == gene_chr)
            & (binary_values["end"] >= interval_start)
            & (binary_values["start"] <= interval_end)
        ]

        R_jk = calculate_single_gene_rp(interval_values, t_k, L, mu)

        peak_rp[gene_id] = R_jk

    return peak_rp


binary_values = pd.read_csv(
    input_file,
    sep="\t",
    names=["chr", "start", "end", "s_ij"],
)

refseq_genes = pd.read_csv(
    refseq_file,
    sep="\t",
    names=["gene_id", "chr", "strand", "tss"],
)

mu = -np.log(L / (3 * delta))

peak_rp = calculate_peak_rp(binary_values, refseq_genes, L=L, mu=mu)

with open(output_file, "w") as f:
    for key in peak_rp.keys():
        f.write("%s,%s\n" % (key, peak_rp[key]))
