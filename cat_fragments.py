#! /usr/bin/python
import argparse

def parse_user_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input-sample',required=True,help='Sample or run name (e.g. PTO025).')
    parser.add_argument('-r','--reference',required=True,help='Path to reference genome.')
    parser.add_argument('-o','--output',required=True,help='Path to output fragment file.')
    return parser

parser=parse_user_input()
ui=parser.parse_args()

chrs=sorted([line.split()[0][1::] for line in open(ui.reference) if line[0]=='>'])

run=ui.input_sample
with open(ui.output,'w') as g:
    for i,ch in enumerate(chrs):
        infile = run+'/outs/fastq_path/'+run+'/'+run+'/'+run+'.'+ch+'.fragments.tsv'
        with open(infile) as f:
            for line in f:
                if line[0]=='#':
                    if i==0:
                        g.write(line)
                else:
                   g.write(line)


