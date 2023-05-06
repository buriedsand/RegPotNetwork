# Regulatory Potential-Derived Gene Regulatory Networks

**Brandon Lukas**

***Department of Biomedical Engineering, University of Illinois at Chicago***

> IEEE-EMBS International Conference on Biomedical and Health Informatics (BHI) 2023

## Introduction


Transcription factors (TFs) are DNA-binding proteins that regulate the expression of target genes in a cell type-specific manner. TFs govern a wide range of cellular processes, such as cell growth, differentiation, and apoptosis, and are occasionally implicated in diseases, including cancer, diabetes, and neurodegenerative disorders. By enhancing our understanding of the complex molecular mechanisms underlying both normal and pathological biology, TF studies can contribute to the development of targeted therapies and personalized medicine. Although chromatin immunoprecipitation followed by sequencing (ChIP-seq) offers the most direct evidence of TF binding, examining binding profiles for every possible TF and cell type combination is experimentally infeasible. To this end, numerous computational approaches have been developed to identify potentially relevant TFs using alternative experimental data sources, such as gene expression profiles or epigenomic data, in conjunction with large-scale databases containing TF binding information.

Motif analysis and deep learning approaches like FactorNet have been developed to predict context-specific TF binding sites. However, motif analysis depends on availability of known motifs, and motif discovery is computationally expensive. Moreover, supervised learning models like FactorNet are impractical for imputing TF binding sites on a large scale due to computational and data availability limitations.

Gene regulatory networks (GRNs) mathematically represent the regulatory relationships between TFs and their target genes. Numerous network-based statistical inference methods, such as AUCell, GSEA, ULM, and WSUM, leverage GRNs and gene expression data to infer TF activity. The accuracy of the predicted activity from these methods immensely depends on the quality of the GRN.

The most widely-used curated GRN is the DoRothEA network. While the connections in this network are supported by different layers of evidence, the network does not reflect all context-specific interactions in particular tissue or cell types of interest, thereby resulting in incomplete or inaccurate predictions.

## Methods

The regulatory potential (RP) model quantifies the regulatory influence of transcriptional regulators, including TFs and chromatin regulators (CRs), on a gene with respect to distance from the gene's transcription start site (TSS). The RP model assumes the regulatory influence of a transcriptional regulator decays exponentially with distance from the TSS.

$$L=\begin{bmatrix}

r_{1,1}   & r_{1,2}   & r_{1,3}   & \cdots    & r_{1,m} \\
r_{2,1}   & r_{2,2}   & r_{2,3}   & \cdots    & r_{2,m} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
r_{n,1}   & r_{n,2}   & r_{n,3}   & \cdots    & r_{n,m} \\

\end{bmatrix}$$

where $r_{i, j}$ is the regulatory potential of transcription factor $i$ on gene $j$.

### Regulatory Potential Calculation

To calculate the regulatory potential of a transcription factor, we first ...

_(continue describing the methods in more detail, including any preprocessing of data, calculations, and algorithms used)_

## Results

In this study, we applied the RP model to construct context-specific GRNs for different tissue and cell types. We compared the performance of our context-specific GRNs to that of the DoRothEA network using various network-based statistical inference methods. Our results demonstrated that ...

_(present your results, including any comparisons, improvements, or interesting findings)_

## Discussion

Our findings show that constructing context-specific GRNs based on the RP model can significantly improve the accuracy of TF activity predictions in specific tissue and cell types. This highlights the importance of considering the context when inferring GRNs and can help to identify novel regulatory mechanisms in complex diseases.

_(discuss any limitations, future work, and implications of your study)_

## Conclusion

