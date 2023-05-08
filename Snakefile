# Main Snakefile

include: "uterine.smk"

with open("assets/factors.txt", "r") as f:
    lines = [line.strip() for line in f]
    TF_LIST = [line.split(",")[0] for line in lines][:3]

# Define the target rule
rule all:
    input:
        expand("data/tmp/outputs/{group}/peak_rp.csv", group=("M", "L")),
        expand("data/tmp/outputs/{group}/chrom_rp.csv", group=("M", "L"))

rule compile_cpp:
    input: "src/cpp/peak_rp.cpp"
    output: "src/cpp/peak_rp"
    shell: "g++ -o {output} {input} -O3 -std=c++17 -pthread"

rule download_tf_chipseq_data:
    output:
        temp("data/tf_chipseq/{tf}.bed")
    params:
        url="https://chip-atlas.dbcls.jp/data/hg19/assembled/Oth.ALL.50.{tf}.AllCell.bed"
    shell:
        "wget -O {output} {params.url}"

rule preprocess_tf_chipseq_data:
    input:
        "data/tf_chipseq/{tf}.bed"
    output:
        protected("data/tf_chipseq_filtered/{tf}.bed")
    params:
        blacklist="assets/blacklists/hg19_blacklist.bed"
    shell:
        "bedtools intersect -a {input} -b {params.blacklist} -v -wa -sorted > {output}"

## Process 2 groups
rule intersect_tf_chipseq_h3k27ac:
    input:
        h3k27ac="inputs/h3k27ac_consensus/{group}.bed",
        tf_chipseq="data/tf_chipseq_filtered/{tf}.bed"
    output:
        temp("data/tmp/tf_chipseq_intersect/{group}/{tf}.bed")
    shell:
        "bedtools intersect -a {input.tf_chipseq} -b {input.h3k27ac} -wa -sorted > {output}"

rule merge_tf_chipseq_intersect:
    input:
        "data/tmp/tf_chipseq_intersect/{group}/{tf}.bed"
    output:
        temp("data/tmp/tf_chipseq_merged/{group}/{tf}.bed")
    shell:
        "bedtools merge -i {input} > {output}"

rule binarize_chipseq_binding:
    input:
        "data/tmp/tf_chipseq_merged/{group}/{tf}.bed"
    output:
        "data/tmp/tf_chipseq_binarized/{group}/{tf}.bed"
    shell:
        "awk 'BEGIN {{FS=OFS=\"\t\"}} !seen[$1, $2, $3]++ {{print $1, $2, $3, 1}}' {input} > {output}"


rule calculate_peak_rp:
    input:
        "data/tmp/tf_chipseq_binarized/{group}/{tf}.bed"
    output:
        "data/tmp/tf_peak_rp/{group}/{tf}.csv"
    params:
        refseq_genes="assets/refseq_genes/hg19_refseq_genes.tsv"
    shell:
        "src/cpp/peak_rp {input} {output} {params.refseq_genes}"

rule aggregate_peak_rp:
    input:
        lambda wildcards: expand(f"data/tmp/tf_peak_rp/{wildcards.group}/{{tf}}.csv", tf=TF_LIST)
    output: "data/tmp/outputs/{group}/peak_rp.csv"
    shell: "python src/python/aggregate_files.py {input} {output}"

rule calculate_chrom_rp:
    input:
        "inputs/h3k27ac_average_signal/{group}.bed"
    output:
        "data/tmp/outputs/{group}/chrom_rp.csv"
    params:
        refseq_genes="assets/refseq_genes/hg19_refseq_genes.tsv"
    shell:
        "src/cpp/peak_rp {input} {output} {params.refseq_genes}"