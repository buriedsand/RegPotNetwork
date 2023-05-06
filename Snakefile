with open("assets/factors.txt", "r") as f:
    lines = [line.strip() for line in f]
    TF_LIST = [line.split(",")[0] for line in lines]
    TF_LIST = TF_LIST[:3]

# Define the target rule
rule all:
    input:
        expand("data/tf_chipseq_filtered/{tf}.bed", tf=TF_LIST)
        # expand("data/output/tf_peak_rp/{tf}.txt", tf=TF_LIST),
        # "data/output/chrom_rp.txt"

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

# rule intersect_tf_chipseq_h3k27ac:
#     input:
#         h3k27ac="data/tmp/h3k27ac_consensus.bed",
#         tf_chipseq="data/tmp/tf_chipseq_filtered/{tf}.bed"
#     output:
#         "data/tmp/tf_chipseq_intersect/{tf}.bed"
#     shell:
#         "bedtools intersect -a {input.tf_chipseq} -b {input.h3k27ac} > {output}"

# rule calculate_peak_rp:
#     input:
#         "data/tmp/tf_chipseq_intersect/{tf}.bed"
#     output:
#         "data/output/tf_peak_rp/{tf}.txt"
#     shell:
#         "python src/python/calculate_peak_rp.py {input} {output}"

# rule calculate_average_h3k27ac_signal:
#     input:
#         "data/tmp/h3k27ac_consensus.bed"
#     output:
#         "data/tmp/h3k27ac_average_signal.bed"
#     shell:
#         "python src/python/calculate_average_signal.py {input} {output}"

# rule calculate_chrom_rp:
#     input:
#         "data/tmp/h3k27ac_average_signal.bed"
#     output:
#         "data/output/chrom_rp.txt"
#     shell:
#         "python src/python/calculate_chrom_rp.py {input} {output}"


###

# rule download_h3k27ac_data:
#     output:
#         "data/tmp/h3k27ac_data.bed"
#     shell:
#         "wget -O {output} <URL>"

# rule preprocess_h3k27ac_data:
#     input:
#         "data/tmp/h3k27ac_data.bed"
#     output:
#         "data/tmp/h3k27ac_filtered.bed"
#     shell:
#         "python src/python/filter_problematic_regions.py {input} {output}"

# rule merge_h3k27ac_peaks:
#     input:
#         "data/tmp/h3k27ac_filtered.bed"
#     output:
#         "data/tmp/h3k27ac_consensus.bed"
#     shell:
#         "python src/python/merge_peaks.py {input} {output}"