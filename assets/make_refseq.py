import os
import pandas as pd

ASSEMBLY = "hg19"

url = (
    f"https://hgdownload.soe.ucsc.edu/goldenPath/{ASSEMBLY}/database/ncbiRefSeq.txt.gz"
)
refseq_file = f"assets/refseq_genes/{ASSEMBLY}_refseq_genes.txt.gz"

os.makedirs(os.path.dirname(refseq_file), exist_ok=True)
os.system(f"wget -O {refseq_file} {url}")

output_file = f"assets/refseq_genes/{ASSEMBLY}_refseq_genes.tsv"

# Read the input file
df = pd.read_csv(
    refseq_file,
    sep="\t",
    header=None,
    names=[
        "bin",
        "name",
        "chr",
        "strand",
        "txStart",
        "txEnd",
        "cdsStart",
        "cdsEnd",
        "exonCount",
        "exonStarts",
        "exonEnds",
        "score",
        "name2",
        "cdsStartStat",
        "cdsEndStat",
        "exonFrames",
    ],
)

# Calculate the TSS for each entry
df["tss"] = df.apply(
    lambda row: row["txStart"] if row["strand"] == "+" else row["txEnd"], axis=1
)

# Calculate the length of each isoform
df["length"] = df["txEnd"] - df["txStart"]

# Sort the data frame by gene symbol and length in descending order
df = df.sort_values(by=["name2", "length"], ascending=[True, False])

# Drop duplicates based on gene symbol and keep the first occurrence (longest isoform)
df = df.drop_duplicates(subset="name2", keep="first")

# Write the output file containing gene symbol, chromosome, strand, and TSS
df[["name2", "chr", "strand", "tss"]].to_csv(
    output_file, sep="\t", header=None, index=False
)

os.system(f"rm {refseq_file}")
