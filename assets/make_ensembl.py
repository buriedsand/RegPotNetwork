import requests
from tqdm import tqdm

genes_file = "assets/unique_genes.txt"
output_file = "assets/refseq_genes/ensembl_genes.tsv"


def get_tss(gene_list, batch_size=500):
    url = "https://www.ensembl.org/biomart/martservice?query="

    def create_xml_query(genes):
        return f"""<?xml version="1.0" encoding="UTF-8"?>
        <!DOCTYPE Query>
        <Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="0" count="" datasetConfigVersion="0.6">
            <Dataset name="hsapiens_gene_ensembl" interface="default">
                <Filter name="hgnc_symbol" value="{','.join(genes)}"/>
                <Attribute name="hgnc_symbol"/>
                <Attribute name="chromosome_name"/>
                <Attribute name="strand"/>
                <Attribute name="transcription_start_site"/>
                <Attribute name="transcript_length"/>
                <Attribute name="transcript_biotype"/>
            </Dataset>
        </Query>"""

    tss_data = []

    for i in tqdm(range(0, len(gene_list), batch_size)):
        genes_batch = gene_list[i : i + batch_size]
        xml_query = create_xml_query(genes_batch)

        response = requests.get(url + xml_query)
        response.raise_for_status()

        data = response.text.strip().split("\n")
        tss_data.extend([row.split("\t") for row in data])

    return tss_data


def find_longest_transcripts(tss_data):
    longest_transcripts = {}
    for row in tss_data:
        gene, chr_name, strand, tss, length, biotype = row
        length = int(length)
        strand = "+" if int(strand) == 1 else "-"

        if biotype != "protein_coding":
            continue

        if gene not in longest_transcripts or length > longest_transcripts[gene][3]:
            longest_transcripts[gene] = (chr_name, strand, tss, length)

    return longest_transcripts


# Load the gene list from the file
with open(genes_file, "r") as file:
    gene_list = [line.strip() for line in file.readlines()]

# Get the TSS data
tss_data = get_tss(gene_list)

# Find the longest transcripts
longest_transcripts = find_longest_transcripts(tss_data)

# Print the results
with open(output_file, "w") as output:
    for gene, (chr_name, strand, tss, length) in longest_transcripts.items():
        output.write(f"{gene}\tchr{chr_name}\t{strand}\t{tss}\n")
