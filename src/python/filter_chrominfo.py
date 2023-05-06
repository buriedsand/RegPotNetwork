import sys
import pandas as pd

input_file = sys.argv[1]
output_file = sys.argv[2]

# Read the chromInfo file
df = pd.read_csv(input_file, sep="\t", header=None, names=["chr", "length"])

# Filter out unwanted contigs
filtered_df = df[
    df["chr"].str.endswith(
        (
            "chr1",
            "chr2",
            "chr3",
            "chr4",
            "chr5",
            "chr6",
            "chr7",
            "chr8",
            "chr9",
            "chr10",
            "chr11",
            "chr12",
            "chr13",
            "chr14",
            "chr15",
            "chr16",
            "chr17",
            "chr18",
            "chr19",
            "chr20",
            "chr21",
            "chr22",
            "chrX",
            "chrY",
            # "chrM",
        )
    )
]

# Write the filtered file
filtered_df.to_csv(output_file, sep="\t", header=False, index=False)
