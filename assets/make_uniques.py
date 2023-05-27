input_file = "data/unweighted_network/10k/ugrn.csv"
output_tf_file = "assets/unique_tfs.txt"
output_genes_file = "assets/unique_genes.txt"

unique_tf = set()
unique_genes = set()

with open(input_file, "r") as f:
    for line in f:
        tf, gene = line.strip().split(",")
        unique_tf.add(tf)
        unique_genes.add(gene)

with open(output_tf_file, "w") as f:
    for tf in unique_tf:
        f.write(tf + "\n")

with open(output_genes_file, "w") as f:
    for gene in unique_genes:
        f.write(gene + "\n")

print("Unique TFs saved to:", output_tf_file)
print("Unique genes saved to:", output_genes_file)
