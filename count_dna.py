#! /usr/bin/python
import gzip
import io
import numpy as np
import pysam
import os


def chrfragments(addressfile, newcbcs):
    fdict = {}
    with io.BufferedReader(gzip.open(addressfile, 'rb')) as f:
        for i, line in enumerate(f):
            dlist = line.decode().split()
            address = '_'.join([dlist[3], dlist[4], newcbcs[i]])
            if address not in fdict:
                fdict[address] = 1
            else:
                fdict[address] += 1
    return fdict


def chrfragments_output(sample, fragfile, ref, frags, chrs):
    with open(fragfile, 'w') as g:
        g.write(f"# id={sample}\n")
        g.write('# description=\n')
        g.write('#\n')
        g.write('# pipeline_name=dna10x\n')
        g.write('# pipeline_version=dna10x-1.0.0\n')
        g.write('#\n')
        g.write(f"# reference_path={ref}\n")
        g.write('# reference_fasta_hash=\n')
        g.write('# reference_gtf_hash=\n')
        g.write('# reference_version=\n')
        g.write('# mkref_version=\n')
        g.write('#\n')
        for ch in chrs:
            g.write(f"# primary_contig={ch}\n")
        for ch, fdict in zip(chrs, frags):
            p1s = np.array([int(k.split('_')[0]) for k in fdict])
            ind = np.argsort(p1s)
            srtkey = np.array(list(fdict.keys()))[ind]
            for sk in srtkey:
                address = sk.split('_')
                st = f"{ch}\t{address[0]}\t{address[1]}\t{address[2]}\t{str(fdict[sk])}\n"
                g.write(st)


def fragments(sample, ref, addressfile, fragfile, chrs, newcbcs):
    frags = {ch: {} for ch in chrs}
    with io.BufferedReader(gzip.open(addressfile, 'rb')) as f:
        for i, line in enumerate(f):
            dlist = line.decode().split()
            ch = dlist[3]
            address = '_'.join([dlist[4], dlist[5], newcbcs[i]])

            if address not in frags[ch]:
                frags[ch][address] = 1
            else:
                frags[ch][address] += 1

    with open(fragfile, 'w') as g:
        g.write(f"# id={sample}\n")
        g.write('# description=\n')
        g.write('#\n')
        g.write('# pipeline_name=dna10x\n')
        g.write('# pipeline_version=dna10x-1.0.0\n')
        g.write('#\n')
        g.write(f"# reference_path={ref}\n")
        g.write('# reference_fasta_hash=\n')
        g.write('# reference_gtf_hash=\n')
        g.write('# reference_version=\n')
        g.write('# mkref_version=\n')
        g.write('#\n')
        for ch in chrs:
            g.write(f"# primary_contig={ch}\n")
        for ch in chrs:
            fdict = frags[ch]
            p1s = np.array([int(k.split('_')[0]) for k in fdict])
            ind = np.argsort(p1s)
            srtkey = np.array(list(fdict.keys()))[ind]
            for sk in srtkey:
                address = sk.split('_')
                st = f"{ch}\t{address[0]}\t{address[1]}\t{address[2]}\t{str(fdict[sk])}\n"
                g.write(st)


def fragments_and_corrected_bam(sample, ref, addressfile, fragfile, chrs, newcbcs, bam_sorted, chromosome, bcset):
    bam_outdir = os.path.join(os.path.dirname(bam_sorted), 'by_chr_bam')
    if not os.path.exists(bam_outdir):
        os.mkdir(bam_outdir)
    out_bam_filename = os.path.basename(bam_sorted).replace('.sorted', '').replace('.bam',
                                                                                   f".{chromosome}.cbc_corrected.bam")
    infile = pysam.AlignmentFile(bam_sorted, 'rb', require_index=True)
    outfile = pysam.AlignmentFile(os.path.join(bam_outdir, out_bam_filename), "wb", template=infile)

    frags = {ch: {} for ch in chrs}
    with io.BufferedReader(gzip.open(addressfile, 'rb')) as f:
        for i, line in enumerate(f):
            dlist = line.decode().split()
            ch = dlist[3]
            address = '_'.join([dlist[4], dlist[5], newcbcs[i]])
            start = int(dlist[4]) - 10
            if start < 0:
                start = 0
            stop = int(dlist[5]) + 10
            if stop >= infile.get_reference_length(ch):
                stop = infile.get_reference_length(ch) - 1
            for read in infile.fetch(contig=ch, start=start, stop=stop):
                readid = ':'.join(read.qname.split(':')[3:7])
                if dlist[0] == readid:
                    if dlist[1] == read.get_tag('CR'):
                        if dlist[2] == read.get_tag('CY'):
                            if newcbcs[i] in bcset:
                                read.set_tag('CB', newcbcs[i])
                            outfile.write(read)
            if address not in frags[ch]:
                frags[ch][address] = 1
            else:
                frags[ch][address] += 1
    outfile.close()
    infile.close()

    with open(fragfile, 'w') as g:
        g.write(f"# id={sample}\n")
        g.write('# description=\n')
        g.write('#\n')
        g.write('# pipeline_name=dna10x\n')
        g.write('# pipeline_version=dna10x-1.0.0\n')
        g.write('#\n')
        g.write(f"# reference_path={ref}\n")
        g.write('# reference_fasta_hash=\n')
        g.write('# reference_gtf_hash=\n')
        g.write('# reference_version=\n')
        g.write('# mkref_version=\n')
        g.write('#\n')
        for ch in chrs:
            g.write(f"# primary_contig={ch}\n")
        for ch in chrs:
            fdict = frags[ch]
            p1s = np.array([int(k.split('_')[0]) for k in fdict])
            ind = np.argsort(p1s)
            srtkey = np.array(list(fdict.keys()))[ind]
            for sk in srtkey:
                address = sk.split('_')
                st = f"{ch}\t{address[0]}\t{address[1]}\t{address[2]}\t{str(fdict[sk])}\n"
                g.write(st)
