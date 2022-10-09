# dna10x

This is a data processing pipeline for scDNA-seq data generating by tagmentation generating using the 10x Genomics scATAC-seq kit. It has the following dependencies:

- python v3.9.13
- pysam v0.19.1 
- bwa v0.7.17
- samtools v1.16.1
- cellranger-atac v2.0.0
- cutadapt v4.1

Try the following to install a conda environment:
```
conda create -n cutadapt -c bioconda -c conda-forge cutadapt python=3.9 bwa pysam samtools numpy
```
