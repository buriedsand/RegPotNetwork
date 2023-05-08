# Main Snakefile

include: "uterine.smk"

with open("assets/factors.txt", "r") as f:
    lines = [line.strip() for line in f]
    TF_LIST = [line.split(",")[0] for line in lines][:3]

# Define the target rule
rule all:
    input:
        expand("data/tmp/tf_chipseq_binarized/M/{tf}_binarized.bed", tf=TF_LIST),
        expand("data/tmp/tf_chipseq_binarized/L/{tf}_binarized.bed", tf=TF_LIST)
        # expand("data/output/tf_peak_rp/{tf}.txt", tf=TF_LIST),
        # "data/output/chrom_rp.txt"

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
        tf_chipseq="data/tmp/tf_chipseq_filtered/{tf}.bed"
    output:
        temp("data/tmp/tf_chipseq_intersect/{group}/{tf}.bed")
    shell:
        "bedtools intersect -a {input.tf_chipseq} -b {input.h3k27ac} > {output}"

rule merge_tf_chipseq_intersect:
    input:
        "data/tmp/tf_chipseq_intersect/{group}/{tf}.bed"
    output:
        temp("data/tmp/tf_chipseq_merged/{group}/{tf}.bed")
    shell:
        "bedtools merge -i {input} > {output}"

rule binarize_chipseq_binding:
    input:
        "data/tmp/tf_chipseq_merged/{group}/{tf}_merged.bed"
    output:
        "data/tmp/tf_chipseq_binarized/{group}/{tf}_binarized.bed"
    shell:
        "awk 'BEGIN {{FS=OFS=\"\t\"}} !seen[$1, $2, $3]++ {{print $1, $2, $3, 1}}' {input} > {output}"


# rule calculate_peak_rp:
#     input:
#         "data/tmp/tf_chipseq_intersect/{tf}.bed"
#     output:
#         "data/output/tf_peak_rp/{tf}.txt"
#     shell:
#         "python src/python/calculate_peak_rp.py {input} {output}"

# rule calculate_chrom_rp:
#     input:
#         "data/tmp/h3k27ac_average_signal.bed"
#     output:
#         "data/output/chrom_rp.txt"
#     shell:
#         "python src/python/calculate_chrom_rp.py {input} {output}"