import sys

input_files = sys.argv[1:-1]
output_file = sys.argv[-1]

with open(output_file, "w") as output:
    for input_file in input_files:
        tf_name = input_file.split("/")[-1].split(".")[0]
        with open(input_file, "r") as f:
            # Skip header line
            next(f)
            for line in f:
                target_gene = line.strip().split("\t")[0]
                if target_gene:
                    output.write(f"{tf_name},{target_gene}\n")
