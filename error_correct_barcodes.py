import argparse
import csv
import os
import sys

from pysam import AlignmentFile

from address_mod import get_raw_cbcfreq_dict
from error_correct_mod import correct_cbc, revcomp


def get_cbc_stats(input_bam, output_dir):
    outstats_path = os.path.join(output_dir,
                                 os.path.basename(input_bam).replace('.bam', '.cbc_stats.csv'))

    cbc_dict = {}
    num_total_reads = 0
    num_reads_with_cbc = 0
    header = ['all_reads', 'mouse_reads', 'all_mapped_reads', 'mapq_30', 'mapq30_mm', 'mapq30_uq', 'mapq30_uq_mm',
              'map30_uq_nuc', 'mapq30_uq_mm_nuc']
    with AlignmentFile(input_bam, 'rb') as inbam:
        for read in inbam:
            num_total_reads += 1
            if read.has_tag('CB'):
                num_reads_with_cbc += 1
                cbc = read.get_tag('CB')
                if cbc in cbc_dict:
                    cbc_dict[cbc][0] += 1
                else:
                    # [all_reads, mouse_reads, all_mapped_reads, mapq_30, mapq30_mm, notdup, notdup_mm, chrM, chrM_mm]
                    cbc_dict[cbc] = [1, 0, 0, 0, 0, 0, 0, 0, 0]
                if read.is_mapped:
                    if read.reference_name.startswith('mm') or read.reference_name.startswith('mouse'):
                        cbc_dict[cbc][1] += 1
                    cbc_dict[cbc][2] += 1
                    if read.mapq >= 30 and read.is_proper_pair:
                        cbc_dict[cbc][3] += 1
                        if not read.is_duplicate:
                            cbc_dict[cbc][5] += 1
                            if 'chrM' not in read.reference_name:
                                cbc_dict[cbc][7] += 1
                        if read.reference_name.startswith('mm') or read.reference_name.startswith('mouse'):
                            cbc_dict[cbc][4] += 1
                            if not read.is_duplicate:
                                cbc_dict[cbc][6] += 1
                            if 'chrM' not in read.reference_name:
                                cbc_dict[cbc][8] += 1

    with open(outstats_path, 'w', newline='') as csvfile:
        statswriter = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_MINIMAL)
        statswriter.writerow(header)
        for cbc, stats in cbc_dict.items():
            stats.insert(0, cbc)
            statswriter.writerow(stats)


def correct_cbc_in_bam(input_bam, cbcfreq_dict, output_dir, barcode_whitelist, stats=False, min_mapq=30, min_tlen=100):
    inbam = AlignmentFile(input_bam, 'rb')
    outbam = AlignmentFile(os.path.join(output_dir, os.path.basename(input_bam).replace('.bam', '.cbc_corrected.bam')),
                           'wb', template=inbam)
    cbc_dict = {}
    num_corrected = 0
    num_total = 0

    for read in inbam:
        raw_cbc = read.get_tag('CR')
        qcbc = read.get_tag('CY')
        is_corrected, cbc = correct_cbc(raw_cbc, qcbc, cbcfreq_dict)
        if cbc in barcode_whitelist:
            read.set_tag('CB', cbc)
            # count read pairs (only read 1, no supplementary, no secondary alignments
            # all_reads, all_mapped_reads, mm_mapped, mapped_filtered, mm_mapped_filtered, 
            if read.is_read1 and not read.is_supplementary and not read.is_secondary:
                num_total += 1
                is_mouse = read.reference_name.startswith('mm')
                if is_corrected:
                    num_corrected += 1
                if cbc in cbc_dict:
                    cbc_dict[cbc][0] += 1
                else:
                    # [all_read_pairs,mapped,mm_mapped,filtered,mm_filtered,filtered_noChrM,mm_filtered_noChrM,filtered_short,filtered_long,mm_filtered_short,mm_filtered_long
                    cbc_dict[cbc] = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                if read.is_mapped:
                    cbc_dict[cbc][1] += 1
                    if is_mouse:
                        cbc_dict[cbc][2] += 1
                    if read.mapq >= min_mapq and read.is_proper_pair and not read.is_duplicate:
                        cbc_dict[cbc][3] += 1
                        if 'chrM' in read.reference_name:
                            cbc_dict[cbc][5] += 1
                        if is_mouse:
                            cbc_dict[cbc][4] += 1
                            if 'chrM' in read.reference_name:
                                cbc_dict[cbc][6] += 1
                        if read.template_length >= min_tlen:
                            cbc_dict[cbc][7] += 1
                            if is_mouse:
                                cbc_dict[cbc][9] += 1
                        else:
                            cbc_dict[cbc][8] += 1
                            if is_mouse:
                                cbc_dict[cbc][10] += 1

        outbam.write(read)

    if stats:
        header = f"#min_mapq={min_mapq}, tlen_cutoff={min_tlen}\n"
        header += "barcode,all_read_pairs,mapped,mm_mapped,filtered,mm_filtered,filtered_noChrM," \
                  "mm_filtered_noChrM,filtered_short,filtered_long,mm_filtered_short,mm_filtered_long\n"
        stats_filename = os.path.basename(input_bam).replace('.bam', '.cbc_stats.txt')
        with open(os.path.join(output_dir, stats_filename), 'w') as stats:
            stats.write(header)
            for cbc, counts in cbc_dict.items():
                stats.write(f"{cbc},{','.join(counts)}\n")

    inbam.close()
    outbam.close()
    return num_corrected, num_total


def parse_user_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--directory', required=True, help='Input directory.')
    parser.add_argument('-b', '--barcodes', required=True, help='Path to 10x barcode whitelist.')
    parser.add_argument('-s', '--samplesheet', required=True, help='Path to sample sheet.')
    parser.add_argument('-o', '--output-directory', required=False, help='Output directory')
    parser.add_argument('-rc', '--revcomp', action='store_true', help='Reverse complement cell barcodes.')
    parser.add_argument('-os', '--output-bc-stats', action='store_true',
                        help='Output barcode frequency file and stats.')
    parser.add_argument('-cs', '--compute-stats-only', action='store_true',
                        help='Compute stats from input BAM files (no error correction, '
                             'assumes CB tags already added to BAM)')
    return parser


def main():

    my_parser = parse_user_input()
    ui = my_parser.parse_args()

    samplesheet = ui.samplesheet
    if os.path.exists(samplesheet):
        with open(samplesheet) as f:
            next(f)
            samples = dict(x.strip().split(',')[1:3] for x in f)
    else:
        sys.exit("Error: can't find sample sheet.")

    barcodes = ui.barcodes

    if os.path.exists(barcodes):
        if ui.revcomp:
            bcset = set([revcomp(line.split()[0]) for line in open(barcodes)])
        else:
            bcset = set([line.split()[0] for line in open(barcodes)])
    else:
        sys.exit("Error: can't find barcode whitelist file.")

    indir = ui.directory
    if not os.path.exists(indir) and not os.path.isdir(indir):
        sys.exit(f"input directory does not exist: {indir}")

    if ui.output_directory:
        outdir = ui.output_directory
    else:
        outdir = indir

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if not os.path.isdir(outdir):
        sys.exit(f"output directory is not a folder: {outdir}")

    for sample in samples:
        print(f"Sample: {sample}")
        bam_file = os.path.join(indir, sample, f"{sample}.bam")

        if ui.compute_stats_only:
            bam_file = os.path.join(indir, f"{sample}.cbc_corrected.dedup.bam")
            if not os.path.exists(bam_file):
                bam_file = os.path.join(indir, f"{sample}.cbc_corrected.sorted.bam")
                if not os.path.exists(bam_file):
                    sys.exit(f"Can't find a valid BAM file for sample: {sample}. BAM file should be named "
                             f"<sample>.cbc_corrected.dedup.bam or <sample>.cbc_corrected.sorted.bam")
            print(f"BAM file: {bam_file}")
            get_cbc_stats(bam_file, outdir)
        else:

            raw_cbcfreq_dict = get_raw_cbcfreq_dict(bam_file)

            outbam_path = correct_cbc_in_bam(bam_file, raw_cbcfreq_dict, outdir, bcset, ui.output_bc_stats)

            if ui.output_bc_stats:
                with open(os.path.join(outdir, f"{sample}.raw_cbc_freqdict.csv"), 'w') as f:
                    for bc, freq in raw_cbcfreq_dict.items():
                        f.write(f"{bc},{freq:.6f}\n")


if __name__ == '__main__':
    main()
