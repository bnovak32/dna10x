# dna10x

This is a data processing pipeline for scDNA-seq data generating by tagmentation generating using the 10x Genomics scATAC-seq or Multiome kit.  It has been tested on sequencing data generated using the Illumina NextSeq 550, Illumina NovaSeq 6000, and Element Aviti. It has the following dependencies:

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

If you already have a fastq, then cellranger-atac is not required, but otherwise the pipeline optionally includes a step where the mkfastq command is used to generate fastq files.

Here's an example command:

```
python dna10x.py --bcl DIRECTORY_WITH_BCL_DATA --samplesheet SAMPLESHEET.csv -d OUTPUT_DIRECTORY --barcodes CELL_BARCODE_STANDARD_LIST.txt -t N_THREADS --reference GENOME.fa -i 1000 -p BARCODE_START_CYCLE -rc -c -ad CTGTCTCTTATACACATCT -sc
```

where CELL_BARCODE_STANDARD_LIST.txt is a one-column table of, for example, 10x Genomics cell barcode sequences.
