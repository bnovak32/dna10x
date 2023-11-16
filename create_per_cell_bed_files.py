import os
import argparse
from pathlib import Path
from pysam import AlignmentFile


def parse_user_input():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-if', '--input-fragments-file', help='Input fragents.tsv (sorted by Barcode).')
    group.add_argument('-ib', '--input-bam-file', help='Input BAM file, sorted by barcode tag.')
    parser.add_argument('-o', '--output-directory', required=True, help="Output directory")
    parser.add_argument('-c', '--min-count', default=1000, help="Minimum number of reads required to output a barcode.")

    return parser


def main():
    my_parser = parse_user_input()
    ui = my_parser.parse_args()

    cbc = ""
    header = ""
    out = None
    count = 0

    is_bam = False

    if ui.input_bam_file is not None:
        filepath = Path(ui.input_bam_file)
        is_bam = True

    else:
        filepath = Path(ui.input_fragments_file)
    file_base = filepath.stem
    file_ext = filepath.suffix
    outdir = Path(ui.output_directory)
    if not filepath.exists() and not filepath.is_file():
        print("Error: can't find input file: ", filepath)
        exit(1)

    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=False)

    if is_bam:
        with AlignmentFile(ui.input_bam_file, 'rb') as inbam:
            print(inbam.filename)

            for read in inbam:
                if not read.has_tag('CB'):
                    print('Error: BAM file has no CB barcode tags.')
                    exit(1)
                if read.get_tag('CB') != cbc:
                    if len(cbc) != 0 and out is not None:
                        out.close()
                        print(f"Barcode: {cbc}\tCount: {count}")
                        if count < ui.min_count:
                            os.remove(outdir.joinpath(f"{file_base}.{cbc}.bam"))
                        count = 0
                    cbc = read.get_tag('CB')
                    out = AlignmentFile(outdir.joinpath(f"{file_base}.{cbc}.bam"), 'wb', template=inbam)
                    out.write(read)
                    count += 1
                else:
                    out.write(read)
                    count += 1
    else:
        with open(filepath) as inbam:
            for line in inbam:
                if not line.startswith('#'):
                    words = line.strip().split('\t')
                    if cbc != words[3]:
                        if len(cbc) != 0 and out is not None:
                            out.close()
                            print(f"Barcode: {cbc}\tCount: {count}")
                            if count < ui.min_count:
                                os.remove(outdir.joinpath(f"{file_base}.{cbc}{file_ext}"))
                            count = 0
                        cbc = words[3]
                        out = open(outdir.joinpath(f"{file_base}.{cbc}{file_ext}"), 'w')
                        out.write(header)
                        out.write(line)
                        count += 1
                    else:
                        out.write(line)
                        count += 1
                else:
                    header += line


if __name__ == '__main__':
    main()
