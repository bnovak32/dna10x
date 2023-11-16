#! /usr/bin/python
import argparse
import os
from glob import glob
from count_dna import fragments, fragments_and_corrected_bam
import subprocess
from address_mod import contig_address
from error_correct_mod import cbccorrect, revcomp
import re
from itertools import zip_longest


def grouper(iterable, n, *, incomplete='fill', fillvalue=None):
    """Collect data into non-overlapping fixed-length chunks or blocks"""
    # grouper('ABCDEFG', 3, fillvalue='x') --> ABC DEF Gxx
    # grouper('ABCDEFG', 3, incomplete='strict') --> ABC DEF ValueError
    # grouper('ABCDEFG', 3, incomplete='ignore') --> ABC DEF
    args = [iter(iterable)] * n
    if incomplete == 'fill':
        return zip_longest(*args, fillvalue=fillvalue)

    if incomplete == 'strict':
        return zip(*args, strict=True)
    if incomplete == 'ignore':
        return zip(*args)
    else:
        raise ValueError('Expected fill, strict, or ignore')


def create_frag_file_for_chr(ch, chrs, ch_cbcdict, ch_qcbcdict, cbcfreq_dict, sample,
                             reference, bamsorted, sample_dir, bcset):
    newcbcs = cbccorrect(ch_cbcdict[ch], ch_qcbcdict[ch], cbcfreq_dict)
    fragfile = os.path.join(sample_dir, 'fragments', f"{sample}.{ch}.fragments.tsv")
    addressfile = os.path.join(sample_dir, 'addresses', f"{sample}.{ch}.address.txt.gz")
    if os.path.exists(bamsorted):
        fragments_and_corrected_bam(sample, reference, addressfile, fragfile, chrs, newcbcs, bamsorted, ch, bcset)
    else:
        fragments(sample, reference, addressfile, fragfile, chrs, newcbcs)


def create_frag_files_for_all_chr(chrs, ch_cbcdict, ch_qcbcdict, cbcfreq_dict, sample, reference, bamsorted,
                                  sample_dir):
    for ch in chrs:
        create_frag_file_for_chr(ch, chrs, ch_cbcdict, ch_qcbcdict, cbcfreq_dict, sample, reference, bamsorted,
                                 sample_dir)


def parse_user_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-bc', '--bcl', required=False, help='Directory with base calling data.')
    parser.add_argument('-s', '--samplesheet', required=True, help='Path to sample sheet.')
    parser.add_argument('-d', '--directory', required=True, help='Output directory.')
    parser.add_argument('-b', '--barcodes', required=True, help='Path to 10x barcode whitelist.')
    parser.add_argument('-t', '--threads', required=True, type=int, help='Number of cpu threads.')
    parser.add_argument('-r', '--reference', required=True, help='Path to reference genome.')
    parser.add_argument('-i', '--insert-size', required=True, type=int, help='Maximum insert size.')
    parser.add_argument('-m', '--minscore', required=True, type=float,
                        help='Alignment score threshold (as a fraction of read length).')
    parser.add_argument('-p', '--barcode-position', type=int, required=True,
                        help='1-indexed cell barcode position in read 2.')
    parser.add_argument('-rc', '--revcomp', help='Reverse complement cell barcodes.', action='store_true')
    parser.add_argument('-sf', '--skip-fastq', action='store_true', help='Skip fastq generation if they already exist.')
    parser.add_argument('-sa', '--skip-align', action='store_true',
                        help='Skip alignment with bwa mem if bam already exists.')
    parser.add_argument('-c', '--cutadapt', action='store_true', help='Run cutadapt on interleaved fastq.')
    parser.add_argument('-ad', '--adapter', required=False, help='Adapter sequence to be trimmed by cutadapt.')
    parser.add_argument('-sd', '--skip-demux', action='store_true', help='Skip production of barcoded fastqs.')
    return parser


def main():

    my_parser = parse_user_input()
    ui = my_parser.parse_args()

    samplesheet = ui.samplesheet
    if os.path.exists(samplesheet):
        with open(samplesheet) as f:
            next(f)
            samples = {line.split(',')[1]: ui.directory for line in f}
    else:
        print("Error: can't find sample sheet.")
        exit()

    threads = int(ui.threads)
    directory = ui.directory
    bcl = ui.bcl
    if not ui.skip_fastq:
        print('Launching Cell Ranger to make fastqs...')
        cmd = f"cellranger-atac mkfastq --run={bcl} --id={directory} --csv={samplesheet} --project={directory}"
        print(cmd)
        os.system(cmd)

    reference = ui.reference
    barcodes = ui.barcodes
    pos = ui.barcode_position

    if os.path.exists(barcodes):
        if ui.revcomp:
            bcset = set([revcomp(line.split()[0]) for line in open(barcodes)])
        else:
            bcset = set([line.split()[0] for line in open(barcodes)])
    else:
        print("Error: can't find barcode whitelist file.")
        exit()

    for sample in samples:
        project = samples[sample]
        sample_dir = os.path.join(directory, "outs", "fastq_path", project, sample)
        os.mkdir(os.path.join(sample_dir, "addresses"))
        os.mkdir(os.path.join(sample_dir, "by_chr_bam"))
        os.mkdir(os.path.join(sample_dir, "fragments"))
        fastqpath = os.path.join(sample_dir, '*_R[123]_*.fastq.gz')
        fastqs = glob(fastqpath)
        R1list = [fq for fq in fastqs if fq.find('_R1_') > -1]
        R2list = [fq for fq in fastqs if fq.find('_R2_') > -1]
        R3list = [fq for fq in fastqs if fq.find('_R3_') > -1]
        R1list.sort()
        R2list.sort()
        R3list.sort()

        j = 1
        procs = []
        fastqouts = []
        bamout = os.path.join(sample_dir, sample + '.bam')
        bamsorted = os.path.join(sample_dir, sample + '.sorted.bam')
        if not ui.skip_demux:
            for r1, r2, r3 in zip(R1list, R2list, R3list):
                fastqout = directory + '/outs/fastq_path/' + project + '/' + sample + '/' + sample + '_' + str(
                    j) + '.fastq.gz'
                cmd = 'python call_demux.py -r1 %(r1)s -r2 %(r2)s -r3 %(r3)s -p %(pos)d| gzip > %(fastqout)s' % vars()
                print(cmd)
                p = subprocess.Popen(cmd, shell=True)
                procs.append(p)
                fastqouts.append(fastqout)
                j += 1
            p_exit = [p.wait() for p in procs]
        else:
            for r1, r2, r3 in zip(R1list, R2list, R3list):
                fastqout = sample_dir + sample + '_' + str(j) + '.fastq.gz'
                fastqouts.append(fastqout)
                j += 1
        if not ui.skip_align:
            fqs = ' '.join(fastqouts)
            print(fastqouts)
            if not ui.cutadapt:
                cmd = f"bwa mem -p -C -M -t {threads} {reference} '<zcat {fqs}' | samtools view -Sb - > {bamout}"
                print(cmd)
                os.system(cmd)
            else:
                adapter = ui.adapter
                cmd = "zcat {fqs} | cutadapt -a {adapter} -A {adapter} --interleaved --cores={threads} -o - - "
                cmd += "| bwa mem -p -C -M -t {threads} {reference} - | samtools view -Sb - > {bamout}"
                print(cmd)
                os.system(cmd)

        addressfile = os.path.join(sample_dir, 'addresses', f"{sample}.address.txt.gz")
        fragfile = os.path.join(sample_dir, 'fragments', f"{sample}.fragments.tsv")
        if reference.endswith('.gz'):
            fai = re.sub('.gz$', '.fai', reference)
        else:
            fai = reference + '.fai'
        if os.path.exists(reference + '.fai'):
            chrs = [line.split()[0] for line in open(fai)]
        else:
            chrs = [line.split()[0][1::] for line in open(reference) if line[0] == '>']
        chrs = sorted(chrs)
        bclen = 16

        addressfile = os.path.join(sample_dir, 'addresses', sample)
        ch_cbcdict, ch_qcbcdict, cbcfreq_dict = contig_address(addressfile, chrs, bamout, bclen, bcset, ui.insert_size,
                                                               ui.minscore)

        for ch in chrs:
            create_frag_file_for_chr(ch, chrs, ch_cbcdict, ch_qcbcdict, cbcfreq_dict, sample, reference, bamsorted,
                                     sample_dir)

        # collect per chromosome fragments into a single file
        fragfile = os.path.join(sample_dir, f"{sample}.fragments.tsv")
        with open(fragfile, 'w') as g:
            for i, ch in enumerate(chrs):
                infile = os.path.join(sample_dir, 'fragments', f"{sample}.{ch}.fragments.tsv")
                with open(infile) as f:
                    for line in f:
                        if line[0] == '#':
                            if i == 0:
                                g.write(line)
                        else:
                            g.write(line)
            # cmd='rm %(infile)s' % vars()
            # print(cmd)
            # os.system(cmd)
        output_filename = os.path.join(sample_dir, f"{sample}.cbc_corrected.bam")
        input_filepattern = os.path.join(sample_dir, 'by_chr_bam', f"{sample}.*.cbc_corrected.bam")
        cmd = f"samtools merge -o {output_filename} {input_filepattern}"

        #  && ' + \
        # 'rm '+directory+'/outs/fastq_path/'+project+'/'+sample+'/'+sample+'.sorted.chr*corrected.bam'
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    main()
