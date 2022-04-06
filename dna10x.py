#! /usr/bin/python
import argparse
import gzip
import io
import os
from os.path import exists
from glob import glob
from pysam import AlignmentFile
from count_dna import addressct
import subprocess

def parse_user_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-bc','--bcl',required=True,help='Directory with base calling data.')
    parser.add_argument('-s','--samplesheet',required=True,help='Path to sample sheet.')
    parser.add_argument('-d','--directory',required=True,help='Output directory.')
    parser.add_argument('-b','--barcodes',required=True,help='Path to 10x barcode whitelist.')
    parser.add_argument('-t','--threads',required=True,type=int,help='Number of cpu threads.')
    parser.add_argument('-r','--reference',required=True,help='Path to reference genome.')
    parser.add_argument('-a','--alignment-score',type=int,required=True,help='Minimum alignment score.')
    parser.add_argument('-i','--insert-size',required=True,type=int,help='Maximum insert size.')
    parser.add_argument('-rc','--revcomp',help='Reverse complement cell barcodes.',action='store_true')
    parser.add_argument('-sf','--skip-fastq',action='store_true',help='Skip fastq generation if they already exist.')
    return parser

parser = parse_user_input()
ui = parser.parse_args()

samplesheet=ui.samplesheet
if exists(samplesheet):
    with open(samplesheet) as f:
        next(f)
        samples = {line.split(',')[1]:ui.directory for line in f}
else:
    print("Error: can't find sample sheet.")
    exit()

print('Launching Cell Ranger to make fastqs...')
threads=ui.threads
directory = ui.directory
bcl=ui.bcl
if not ui.skip_fastq:
    cmd = 'cellranger-atac mkfastq --run=%(bcl)s --id=%(directory)s --csv=%(samplesheet)s --project=%(directory)s' % vars()
    print(cmd)
    os.system(cmd)

reference = ui.reference
barcodes = ui.barcodes
for sample in samples:
    project = samples[sample]
    fastqpath = directory+'/outs/fastq_path/'+project+'/'+sample+'/*'
    fastqs = glob(fastqpath)
    R1list = [fq for fq in fastqs if fq.find('_R1_')>-1]
    R2list = [fq for fq in fastqs if fq.find('_R2_')>-1]
    R3list = [fq for fq in fastqs if fq.find('_R3_')>-1]
    R1list.sort()
    R2list.sort()
    R3list.sort()
    
    j=1
    procs=[]
    fastqouts=[]
    revcomp=0
    if ui.revcomp:
        revcomp=1
    for r1,r2,r3 in zip(R1list,R2list,R3list):
        fastqout = directory+'/outs/fastq_path/'+project+'/'+sample+'/'+sample+'_'+str(j)+'.fastq.gz'
        cmd='python call_demux.py -r1 %(r1)s -r2 %(r2)s -r3 %(r3)s -b %(barcodes)s -r %(revcomp)d| gzip > %(fastqout)s' % vars()
        print(cmd)
        p=subprocess.Popen(cmd,shell=True)
        procs.append(p)
        fastqouts.append(fastqout)
        j+=1
    p_exit = [p.wait() for p in procs]

    bamout=directory+'/outs/fastq_path/'+project+'/'+sample+'/'+sample+'.bam'
    fqs = ' '.join(fastqouts)
    cmd = "bwa mem -p -C -M -t %(threads)d %(reference)s '<zcat %(fqs)s' | samtools view -Sb - > %(bamout)s" % vars()
    print(cmd)
    os.system(cmd)
    address=directory+'/outs/fastq_path/'+project+'/'+sample+'/'+sample+'.address.txt.gz'    
    with gzip.open(address, 'wb') as g:
        with AlignmentFile(bamout,'rb') as f:
            pair = 0
            plist=[]
            for read in f:
                score = int(read.get_tag('AS'))
                isize = int(read.tlen)
                if score>ui.alignment_score and abs(isize)<ui.insert_size:   
                    readid = ':'.join(read.qname.split(':')[3:7])
                    cbc = read.get_tag('BC') 
                    ch = f.getrname(read.reference_id)
                    p1 = read.reference_start
                    if pair==0:
                        plist=[readid,cbc,ch,p1,score,isize]
                        pair=1
                    else:
                        if readid==plist[0]:
                            newline = readid+'\t'+cbc+'\t'+ch+'\t'+str(plist[3])+'\t'+str(p1)+'\t'+str(plist[4])+'\t'+str(score)+'\t'+str(plist[5])+'\n'
                            g.write(newline.encode())
                            pair=0
                        else:
                            plist=[readid,cbc,ch,p1,score,isize]
    counts = directory+'/outs/fastq_path/'+project+'/'+sample+'/'+sample+'.addressct.txt'
    stats = directory+'/outs/fastq_path/'+project+'/'+sample+'/'+sample+'.stats.txt'
    addressct(address,counts,stats)


