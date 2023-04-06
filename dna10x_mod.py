#! /usr/bin/python
import argparse
import gzip
import io
import os
from os.path import exists
from glob import glob
import sys
from pysam import AlignmentFile
from count_dna import chrfragments,chrfragments_output,fragments,fragments_and_corrected_bam
import subprocess
from address_mod import contig_address
from error_correct import cbccorrect,revcomp
from functools import partial
from multiprocessing import Pool
import re


def parse_user_input():
	parser = argparse.ArgumentParser()
	parser.add_argument('-bc','--bcl',required=False,help='Directory with base calling data.')
	parser.add_argument('-s','--samplesheet',required=True,help='Path to sample sheet.')
	parser.add_argument('-d','--directory',required=True,help='Output directory.')
	parser.add_argument('-b','--barcodes',required=True,help='Path to 10x barcode whitelist.')
	parser.add_argument('-t','--threads',required=True,type=int,help='Number of cpu threads.')
	parser.add_argument('-r','--reference',required=True,help='Path to reference genome.')
	parser.add_argument('-i','--insert-size',required=True,type=int,help='Maximum insert size.')
	parser.add_argument('-m','--minscore',required=True,type=float,help='Alignment score threshold (as a fraction of read length).')
	parser.add_argument('-p','--barcode-position',type=int,required=True,help='1-indexed cell barcode position in read 2.')
	parser.add_argument('-rc','--revcomp',help='Reverse complement cell barcodes.',action='store_true')
	parser.add_argument('-sf','--skip-fastq',action='store_true',help='Skip fastq generation if they already exist.')
	parser.add_argument('-sa','--skip-align',action='store_true',help='Skip alignment with bwa mem if bam already exists.')
	parser.add_argument('-c','--cutadapt',action='store_true',help='Run cutadapt on interleaved fastq.')
	parser.add_argument('-ad','--adapter',required=False,help='Adapter sequence to be trimmed by cutadapt.')
	parser.add_argument('-sd','--skip-demux',action='store_true',help='Skip production of barcoded fastqs.')
	parser.add_argument('-chr','--chromosomes', help='list of chromosomes to process')
	return parser

parser = parse_user_input()
ui = parser.parse_args()

samplesheet=ui.samplesheet
if exists(samplesheet):
	with open(samplesheet) as f:
		next(f)
		samples = dict(x.strip().split(',')[1:3] for x in f)
else:
	print("Error: can't find sample sheet.")
	exit()

threads=int(ui.threads)
directory = ui.directory
bcl=ui.bcl
if not ui.skip_fastq:
	print('Launching Cell Ranger to make fastqs...')
	cmd = 'cellranger-atac mkfastq --run=%(bcl)s --id=%(directory)s --csv=%(samplesheet)s --project=%(directory)s' % vars()
	print(cmd)
	# os.system(cmd)

reference = ui.reference
barcodes = ui.barcodes
pos = ui.barcode_position
for sample in samples:
	print("sample: ", sample)
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
	bamout=os.path.join(directory,'outs','fastq_path',project,sample,sample+'.bam')
	bamsorted = os.path.join(directory,'outs','fastq_path',project,sample,sample+'.sorted.bam')
	if not ui.skip_demux:
		for r1,r2,r3 in zip(R1list,R2list,R3list):
			fastqout = directory+'/outs/fastq_path/'+project+'/'+sample+'/'+sample+'_'+str(j)+'.fastq.gz'
			cmd='python call_demux.py -r1 %(r1)s -r2 %(r2)s -r3 %(r3)s -p %(pos)d| gzip > %(fastqout)s' % vars()
			print(cmd)
			p=subprocess.Popen(cmd,shell=True)
			procs.append(p)
			fastqouts.append(fastqout)
			j+=1
		p_exit = [p.wait() for p in procs]
	else:
		for r1,r2,r3 in zip(R1list,R2list,R3list):
			fastqout = directory+'/outs/fastq_path/'+project+'/'+sample+'/'+sample+'_'+str(j)+'.fastq.gz'
			fastqouts.append(fastqout)
			j+=1
	if not ui.skip_align:
		fqs = ' '.join(fastqouts)
		if not ui.cutadapt:
			cmd = "bwa mem -p -C -M -t %(threads)d %(reference)s '<zcat %(fqs)s' | samtools view -Sb - > %(bamout)s" % vars()
			print(cmd)
			# os.system(cmd)
		else:
			adapter = ui.adapter
			cmd = "zcat %(fqs)s | cutadapt -a %(adapter)s -A %(adapter)s --interleaved --cores=%(threads)s -o - - | bwa mem -p -C -M -t %(threads)s %(reference)s - | samtools view -Sb - > %(bamout)s" % vars()
			print(cmd)
			# os.system(cmd)

	# follow symlink if any
	if not os.path.exists(bamout):
		sys.exist('BAM file does not exist: {}'.format(bamout))
	
	addressfile=directory+'/outs/fastq_path/'+project+'/'+sample+'/'+sample+'.address.txt.gz'	
	fragfile=directory+'/outs/fastq_path/'+project+'/'+sample+'/'+sample+'.fragments.tsv'
	if ui.chromosomes is not None and os.path.exists(ui.chromosomes):
		with open(ui.chromosomes) as f:
			chrs=[line.split()[0] for line in f]
	else:
		if reference.endswith('.gz'):
			fai = re.sub('.gz$', '.fai', reference)
		else:
			fai = reference+'.fai'
		if os.path.exists(reference+'.fai'):
			chrs = [line.split()[0] for line in open(fai)]
		else:
			chrs=[line.split()[0][1::] for line in open(reference) if line[0]=='>']
	chrs=sorted(chrs)
	bclen=16
	if exists(barcodes):
		if ui.revcomp:
			bcset = set([revcomp(line.split()[0]) for line in open(barcodes)])
		else:
			bcset = set([line.split()[0] for line in open(barcodes)])
	else:
		print("Error: can't find barcode whitelist file.")
		exit()   
	addressfile = directory+'/outs/fastq_path/'+project+'/'+sample+'/'+sample
	ch_cbcdict,ch_qcbcdict,cbcfreq_dict=contig_address(addressfile,chrs,bamout,bcset,ui.insert_size,ui.minscore)
	for ch in chrs:
		newcbcs=cbccorrect(ch_cbcdict[ch],ch_qcbcdict[ch],cbcfreq_dict)
		fragfile=directory+'/outs/fastq_path/'+project+'/'+sample+'/'+sample+'.'+ch+'.fragments.tsv'
		addressfile = directory+'/outs/fastq_path/'+project+'/'+sample+'/'+sample+'.'+ch+'.address.txt.gz'
		if os.path.exists(bamsorted):
			fragments_and_corrected_bam(sample,reference,addressfile,fragfile,chrs,newcbcs,bamsorted,ch)
		else:
			fragments(sample,reference,addressfile,fragfile,chrs,newcbcs)
	# collect per chromosome fragments into a single file
	fragfile = directory+'/outs/fastq_path/'+project+'/'+sample+'/'+sample+'.fragments.tsv'
	with open(fragfile,'w') as g:
		for i,ch in enumerate(chrs):
			infile = directory+'/outs/fastq_path/'+project+'/'+sample+'/'+sample+'.'+ch+'.fragments.tsv'
			with open(infile) as f:
				for line in f:
					if line[0]=='#':
						if i==0:
							g.write(line)
					else:
						g.write(line)
			cmd='rm %(infile)s' % vars()
			print(cmd)
			os.system(cmd)
	cmd = 'samtools merge -n -o ' + \
			directory+'/outs/fastq_path/'+project+'/'+sample+'/'+sample+'.cbc_corrected.bam ' + \
			directory+'/outs/fastq_path/'+project+'/'+sample+'/'+sample+'.sorted.chr*corrected.bam && ' + \
			'rm '+directory+'/outs/fastq_path/'+project+'/'+sample+'/'+sample+'.sorted.chr*corrected.bam'
	print(cmd)
	os.system(cmd)

	
