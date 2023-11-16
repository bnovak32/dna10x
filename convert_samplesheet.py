import json
import argparse


def convert_10x_samplesheet(sheet_10x, indexlist_json, max_lanes=1, sample_project=''):
    index_dict = {}
    with open(indexlist_json) as ilj:
        for index in json.load(ilj):
            index_dict[index[0]] = index[1]
    with open(sheet_10x, 'r') as f:
        with open(sheet_10x.replace('.csv', '.bcl2fastq.csv'), 'w') as fo:
            for line in f:
                if line.startswith('Lane'):
                    fo.write('[Data]\n')
                    fo.write("Lane,Sample_ID,Sample_Name,Sample_Project,index,index2\n")
                else:
                    lanestr, sample, indexcode = line.strip().split(',')
                    for barcode in index_dict[indexcode]:
                        if lanestr == '*':
                            lanes = range(1, max_lanes+1)
                        elif '-' in lanestr:
                            lane_range = lanestr.split('-')
                            lanes = range(lane_range[0], lane_range[1]+1)
                        else:
                            lanes = [lanestr]
                        for lane in lanes:
                            outl = ','.join([str(lane), sample, sample, sample_project, barcode, ''])
                            fo.write(f"{outl}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sample-sheet', type=str, required=True,
                        help='Input samplesheet in 10x format for mkfastq')
    parser.add_argument('-b', '--barcode-mapping-file', type=str, required=True,
                        help='csv file mapping 10x index codes to index sequences')
    parser.add_argument('-l', '--max-lanes', type=int, choices=[1, 2, 4, 8], default=1,
                        help='maximum number of lanes on instrument')
    parser.add_argument('-p', '--sample-project', type=str, default='',
                        help='project name')

    args = parser.parse_args()
    convert_10x_samplesheet(args.sample_sheet, args.barcode_mapping_file, args.max_lanes, args.sample_project)
    print(args)


if __name__ == '__main__':
    main()


