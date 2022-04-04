#! /usr/bin/python
import argparse
from demux_atac import demux
from os.path import exists

def parse_user_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r1','--read1-fastq',required=True,help='Path to read 1 fastq file.')
    parser.add_argument('-r2','--read2-fastq',required=True,help='Path to read 2 fastq file.')
    parser.add_argument('-r3','--read3-fastq',required=True,help='Path to read 3 fastq file.')
    parser.add_argument('-b','--barcodes',required=True,help='Path to 10x barcode whitelist.') 
    parser.add_argument('-r','--revcomp',required=True,help='Reverse complement cell barcodes.')
    return parser

def revcomp(seq):
    sdict={'A':'T','T':'A','G':'C','C':'G'}
    newseq = [sdict[s] for s in seq]
    return ''.join(newseq[::-1])


parser = parse_user_input()
ui = parser.parse_args()

if exists(ui.barcodes):
    if ui.revcomp=='1':
        bcset = set([revcomp(line.split()[0]) for line in open(ui.barcodes)])
    else:
        bcset = set([line.split()[0] for line in open(ui.barcodes)])
else:
    print("Error: can't find barcode whitelist file.")
    exit()

demux(ui.read1_fastq,ui.read2_fastq,ui.read3_fastq,bcset)

