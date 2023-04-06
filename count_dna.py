#! /usr/bin/python
import gzip
import io
import os
import numpy as np
from os.path import exists
import pysam

def chrfragments(addressfile,newcbcs):
    fdict={}
    with io.BufferedReader(gzip.open(addressfile,'rb')) as f:
        for i,line in enumerate(f):
            dlist = line.decode().split()
            address='_'.join([dlist[3],dlist[4],newcbcs[i]])
            if address not in fdict:
                fdict[address]=1
            else:
                fdict[address]+=1
    return fdict

def chrfragments_output(sample,fragfile,ref,frags,chrs):
    with open(fragfile,'w') as g:
        g.write('# id=%(sample)s\n' % vars())
        g.write('# description=\n')
        g.write('#\n')
        g.write('# pipeline_name=dna10x\n')
        g.write('# pipeline_version=dna10x-1.0.0\n')
        g.write('#\n')
        g.write('# reference_path=%(ref)s\n' % vars())
        g.write('# reference_fasta_hash=\n')
        g.write('# reference_gtf_hash=\n')
        g.write('# reference_version=\n')
        g.write('# mkref_version=\n')
        g.write('#\n')
        for ch in chrs:
            g.write('# primary_contig=%(ch)s\n' % vars())
        for ch,fdict in zip(chrs,frags):
            p1s=np.array([int(k.split('_')[0]) for k in fdict])
            ind=np.argsort(p1s)
            srtkey=np.array(list(fdict.keys()))[ind]
            for sk in srtkey:
                address=sk.split('_')
                st=ch+'\t'+address[0]+'\t'+address[1]+'\t'+address[2]+'\t'+str(fdict[sk])+'\n'
                g.write(st)



def fragments(sample,ref,addressfile,fragfile,chrs,newcbcs):
    frags={ch:{} for ch in chrs}
    with io.BufferedReader(gzip.open(addressfile,'rb')) as f:
        for i,line in enumerate(f):
            dlist = line.decode().split()
            ch=dlist[3]
            address='_'.join([dlist[4],dlist[5],newcbcs[i]])
            
            if address not in frags[ch]:
                frags[ch][address]=1
            else:
                frags[ch][address]+=1
    
    with open(fragfile,'w') as g:
        g.write('# id=%(sample)s\n' % vars())
        g.write('# description=\n')
        g.write('#\n')
        g.write('# pipeline_name=dna10x\n')
        g.write('# pipeline_version=dna10x-1.0.0\n')
        g.write('#\n')
        g.write('# reference_path=%(ref)s\n' % vars())
        g.write('# reference_fasta_hash=\n')
        g.write('# reference_gtf_hash=\n')
        g.write('# reference_version=\n')
        g.write('# mkref_version=\n')
        g.write('#\n')
        for ch in chrs:
            g.write('# primary_contig=%(ch)s\n' % vars())
        for ch in chrs:
            fdict=frags[ch]
            p1s=np.array([int(k.split('_')[0]) for k in fdict])
            ind=np.argsort(p1s)
            srtkey=np.array(list(fdict.keys()))[ind]
            for sk in srtkey:
                address=sk.split('_')
                st=ch+'\t'+address[0]+'\t'+address[1]+'\t'+address[2]+'\t'+str(fdict[sk])+'\n'
                g.write(st)

def fragments_and_corrected_bam(sample,ref,addressfile,fragfile,chrs,newcbcs,bam_sorted,chromosome):
    infile = pysam.AlignmentFile(bam_sorted, 'rb', require_index=True)
    outfile = pysam.AlignmentFile(bam_sorted.replace('.bam', '.'+chromosome+'.cbc_corrected.bam').replace('.sorted', ''), "wb", template=infile) 
    
    frags={ch:{} for ch in chrs}
    with io.BufferedReader(gzip.open(addressfile,'rb')) as f:
        for i,line in enumerate(f):
            dlist = line.decode().split()
            ch=dlist[3]
            address='_'.join([dlist[4],dlist[5],newcbcs[i]])
            for read in infile.fetch(contig=ch, start=int(dlist[4])-10, stop=int(dlist[5])+10):
                readid = ':'.join(read.qname.split(':')[3:7])
                if dlist[0] == readid:
                    if dlist[1] == read.get_tag('CR'):
                        if dlist[2] == read.get_tag('CY'):
                            read.set_tag('CB', newcbcs[i])
                            outfile.write(read)
            if address not in frags[ch]:
                frags[ch][address]=1
            else:
                frags[ch][address]+=1
    outfile.close()
    infile.close()

    with open(fragfile,'w') as g:
        g.write('# id=%(sample)s\n' % vars())
        g.write('# description=\n')
        g.write('#\n')
        g.write('# pipeline_name=dna10x\n')
        g.write('# pipeline_version=dna10x-1.0.0\n')
        g.write('#\n')
        g.write('# reference_path=%(ref)s\n' % vars())
        g.write('# reference_fasta_hash=\n')
        g.write('# reference_gtf_hash=\n')
        g.write('# reference_version=\n')
        g.write('# mkref_version=\n')
        g.write('#\n')
        for ch in chrs:
            g.write('# primary_contig=%(ch)s\n' % vars())
        for ch in chrs:
            fdict=frags[ch]
            p1s=np.array([int(k.split('_')[0]) for k in fdict])
            ind=np.argsort(p1s)
            srtkey=np.array(list(fdict.keys()))[ind]
            for sk in srtkey:
                address=sk.split('_')
                st=ch+'\t'+address[0]+'\t'+address[1]+'\t'+address[2]+'\t'+str(fdict[sk])+'\n'
                g.write(st)