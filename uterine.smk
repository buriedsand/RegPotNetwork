with open("inputs/samples.txt", "r") as f:
    lines1 = [line.strip() for line in f]
    SAMPIDS = [line.split(",")[0] for line in lines1]

def get_sample_file(wildcards):
    idx = SAMPIDS.index(wildcards[0])
    return lines1[idx].split(",")[2]

rule preprocess_h3k27ac_data:
    input:
        "inputs/h3k27ac_data"
    output:
        "inputs/h3k27ac_preprocessed/{sampid}.bed"
    params:
        sample_file=get_sample_file,
        blacklist="assets/blacklists/hg19_blacklist.bed"
    shell:
        "bedtools intersect -a inputs/h3k27ac_data/{params.sample_file} -b {params.blacklist} -v -wa -sorted > {output}"

rule merge_h3k27ac_peaks:
    input:
        files=expand("inputs/h3k27ac_preprocessed/{sampid}.bed", sampid=SAMPIDS),
        sample_sheet="inputs/mumerge.input"
    output:
        "inputs/h3k27ac_merged"
    shell:
        "mkdir {output} && python mumerge/mumerge/mumerge.py -i {input.sample_sheet} -o {output}/ -s -v"

rule consense_h3k27ac_peaks:
    input:
        "inputs/h3k27ac_merged"
    output:
        "inputs/h3k27ac_consensus"
    shell:
        "python inputs/make_consensus.py"

rule calculate_average_h3k27ac_signal:
    input:
        "inputs/h3k27ac_consensus"
    output:
        "inputs/h3k27ac_average_signal"
    shell:
        "bash inputs/calculate_average_signal.sh"