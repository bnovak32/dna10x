#! /usr/bin/python
import argparse
from demux_atac import demux


def parse_user_input():
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-r1', '--read1-fastq', required=True,
                           help='Path to read 1 fastq file.')
    my_parser.add_argument('-r2', '--read2-fastq', required=True,
                           help='Path to read 2 fastq file.')
    my_parser.add_argument('-r3', '--read3-fastq', required=True,
                           help='Path to read 3 fastq file.')
    my_parser.add_argument('-p', '--barcode-position', type=int, required=True,
                           help='1-indexed position of barcode in read 2.')
    my_parser.add_argument('-sra', action=argparse.BooleanOptionalAction, type=bool,
                           required=False, default=False, help='Data came from SRA')
    return my_parser


def main():
    parser = parse_user_input()
    ui = parser.parse_args()

    pos = ui.barcode_position - 1
    demux(ui.read1_fastq, ui.read2_fastq, ui.read3_fastq, pos, sra=ui.sra)


if __name__ == '__main__':
    main()
