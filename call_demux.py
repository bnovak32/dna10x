#! /usr/bin/python
import argparse
from demux_atac import demux
from os.path import exists

def parse_user_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r1','--read1-fastq',required=True,help='Path to read 1 fastq file.')
    parser.add_argument('-r2','--read2-fastq',required=True,help='Path to read 2 fastq file.')
    parser.add_argument('-r3','--read3-fastq',required=True,help='Path to read 3 fastq file.')
    return parser

parser = parse_user_input()
ui = parser.parse_args()

demux(ui.read1_fastq,ui.read2_fastq,ui.read3_fastq)

