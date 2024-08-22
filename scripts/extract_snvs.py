#######################
#### import module ####
#######################
import os
import sys
import getopt

import re
import pysam
from itertools import chain
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

script_name = os.path.basename(__file__)

##########
# ReadMe #
##########
readME = f'''
USAGE:\t{script_name} -i <inputfile> -o <output directory>

Arguments:
\t-i, --input           input bam file
\t-o, --outputdir       output directory

Optional:
\t-b, --basequal        the minimum quality score of a base, eg: 30, default: 30
\t-r, --region          the expected variant region, eg: 100,200, default: 0,0
\t-p, --prefix          output prefix

Flags:
\t-h, --help            help information
\t-v, --version         version information
''' 

###############
# subfunction #
###############
def versions():
    verStr = f"Program:\t{script_name}\nVersion:\t1.0"
    print(verStr)

def usage_info():
    versions()
    print(readME)

def get_snvs(cigar_str: str, md_tag: str, mapping_start: int, read_seq: str, read_qual: str):
    parts_cigar = re.findall(r'(\d+|[A-Z])', cigar_str)
    parts_md = re.findall(r'(\d+|[A-Z]|\^[A-Z]+)', md_tag)

    read_pos = 0
    first_softclip_end = 0
    insertion_ranges = []
    for i in range(len(parts_cigar)):
        match parts_cigar[i]:
            case 'I':
                if parts_cigar[i-1].isdigit():
                    insertion_ranges.append((read_pos + 1, read_pos + int(parts_cigar[i-1])))
                    read_pos = read_pos + int(parts_cigar[i-1])
            case 'M':
                if parts_cigar[i-1].isdigit():
                    read_pos = read_pos + int(parts_cigar[i-1])
            case 'S':
                if parts_cigar[i-1].isdigit() and i == 1:
                    read_pos = read_pos + int(parts_cigar[i-1])
                    if first_softclip_end == 0:
                        first_softclip_end = read_pos

    read_pos = first_softclip_end
    map_pos = mapping_start
    snvs = []
    for md in parts_md:
        if md.isdigit():
            read_pos = read_pos + int(md)
            read_pos = adjust_read_pos(read_pos, insertion_ranges)
            map_pos = map_pos + int(md)
        elif md.isalpha():
            read_pos = read_pos + 1
            map_pos = map_pos + 1
            snvs.append((map_pos, md, read_seq[read_pos-1], read_qual[read_pos-1]))
        elif md.startswith('^'):
            map_pos = map_pos + len(md) - 1

    return snvs

def adjust_read_pos(read_pos: int, insertion_ranges: list):
    if insertion_ranges:
        for range in insertion_ranges:
            if read_pos >= range[0] and read_pos <= range[1]:
                return range[1] - range[0] + 1 + read_pos
                break
        return read_pos
    else:
        return read_pos

def plot_frequency(cov_file: str, png_file: str, var_region: str):
    df = pd.read_csv(cov_file, sep = '\t', header = 0)

    positions = df['pos'].tolist()

    # Calculate frequencies and percentages excluding 'ref'
    A_freq = df.apply(lambda row: row['A'] if row['ref'] != 'A' else 0, axis=1)
    C_freq = df.apply(lambda row: row['C'] if row['ref'] != 'C' else 0, axis=1)
    G_freq = df.apply(lambda row: row['G'] if row['ref'] != 'G' else 0, axis=1)
    T_freq = df.apply(lambda row: row['T'] if row['ref'] != 'T' else 0, axis=1)

    total_counts = df[['A', 'C', 'G', 'T']].sum(axis=1)
    A_pct = A_freq / total_counts * 100
    C_pct = C_freq / total_counts * 100
    G_pct = G_freq / total_counts * 100
    T_pct = T_freq / total_counts * 100

    total_pct = [a + t + c + g for a, t, c, g in zip(A_pct, T_pct, C_pct, G_pct)]

    heatmap_data = np.array([total_pct, A_pct, C_pct, G_pct, T_pct])

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(50, 10), sharex=True)

    ax1.plot(np.arange(len(positions)), A_freq, marker='o', markersize=3, linestyle='-', label='A')
    ax1.plot(np.arange(len(positions)), C_freq, marker='o', markersize=3, linestyle='-', label='C')
    ax1.plot(np.arange(len(positions)), G_freq, marker='o', markersize=3, linestyle='-', label='G')
    ax1.plot(np.arange(len(positions)), T_freq, marker='o', markersize=3, linestyle='-', label='T')
    ax1.set_xticks(np.arange(len(positions)))
    ax1.set_xticklabels(positions, rotation=90, fontsize=8)
    ax1.set_ylabel('Frequency', fontsize=16)
    ax1.set_title('Frequency of Nucleotides (A, C, G, T) excluding Ref by Position', fontsize=16)
    ax1.legend()

    #ax2.plot(positions, A_pct, marker='o', markersize=3, linestyle='-', label='A')
    #ax2.plot(positions, C_pct, marker='o', markersize=3, linestyle='-', label='C')
    #ax2.plot(positions, G_pct, marker='o', markersize=3, linestyle='-', label='G')
    #ax2.plot(positions, T_pct, marker='o', markersize=3, linestyle='-', label='T')
    #ax2.set_ylabel('Percentage')
    #ax2.set_xlabel('Position (pos)')
    #ax2.set_title('Percentage of Nucleotides (A, C, G, T) excluding Ref by Position')
    #ax2.legend()

    im = ax2.imshow(heatmap_data, cmap='Reds', aspect='auto', interpolation='none', vmin=0, vmax=10)

    # Customize the plot
    ax2.set_xticks(np.arange(len(positions)))
    ax2.set_xticklabels(positions, rotation=90, fontsize=8)
    ax2.set_yticks(np.arange(5))
    ax2.set_yticklabels(['Total', 'A', 'C', 'G', 'T'], fontsize=16)
    ax2.set_title('Percentage of Nucleotides (A, C, G, T) excluding Ref by Position', fontsize=16)

    [region_start, region_end] = var_region.split(',')
    region_start = int(region_start)
    region_end = int(region_end)
    if (region_start != 0 | region_end != 0) & (region_start < region_end):
        region_start_idx = 0
        region_end_idx = 0
        for index, pos in enumerate(positions):
            if region_start >= pos:
                region_start_idx = index
            if region_end >= pos:
                region_end_idx = index
        ax2.axvspan(region_start_idx + 1, region_end_idx + 1, color = 'grey', alpha = 0.2)

    # Add colorbar
    cbar = fig.colorbar(im, ax=ax2, orientation='horizontal', pad=0.2)
    cbar.set_label('Percentage (%)')
    cbar.set_ticks([0, 5, 10])

    plt.savefig(png_file)


#################
# main function #
#################

def main(argvs):
    inputFile = ''
    outputDir = ''
    baseQual = 30
    varRegion = "0,0"
    outputPrefix = ''

    try:
        opts, args = getopt.getopt(argvs,
                                   "vhi:o:b:r:p:",
                                   ["version", "help", "input=", "outputdir=", "basequal=", "region=", "prefix="])
        if len(opts) == 0:
            usage_info()
            sys.exit(2)
    except getopt.GetoptError as err:
        print(f"Error: {str(err)}")
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-v", "--version"):
            verbose = True
            versions()
            sys.exit()
        elif opt in ("-h", "--help"):
            usage_info()
            sys.exit()
        elif opt in ("-i", "--input"):
            inputFile = arg
        elif opt in ("-o", "--outputdir"):
            outputDir = arg
        elif opt in ("-b", "--basequal"):
            baseQual = arg
        elif opt in ("-r", "--region"):
            varRegion = arg            
        elif opt in ("-p", "--prefix"):
            outputPrefix = arg
        else:
            assert False, "unhandled option"

    if inputFile == "": 
        sys.exit("Please use the correct arguments, option -i missing")
    if outputDir == "": 
        sys.exit("Please use the correct arguments, option -o missing")

    if outputPrefix == '':
        fileName = os.path.basename(inputFile)
        outputPrefix = fileName.split('.')[0]

    # open bam
    bam = pysam.AlignmentFile(inputFile, "rb")

    # collect all the snvs in all the reads
    all_snvs = []
    chr_id = bam.references[0]
    for record in bam:
        if record.has_tag("MD"):
            all_snvs.append(get_snvs(record.cigarstring, 
                                     record.get_tag("MD"), 
                                     record.reference_start, 
                                     record.query_sequence,
                                     record.query_qualities))

    # collect all the genomic positions of all the snvs
    all_snvs_uniq = list(set(list(chain.from_iterable(all_snvs))))
    all_refs_uniq = [(item[0], item[1]) for item in all_snvs_uniq]
    all_refs_uniq = list(set(all_refs_uniq))
    all_refs_uniq.sort(key=lambda x: x[0])

    # collect all the coverages of all the snvs
    snvs_cov = []
    for pos, ref in all_refs_uniq:
        # quality_threshold is the minimum quality score (in phred) a base has to reach to be counted.
        snv_cov = bam.count_coverage(chr_id, pos-1, pos, quality_threshold=int(baseQual))
        snvs_cov.append([cov[0] for cov in snv_cov])
    all_snvs_uniq_cov = dict(zip(all_refs_uniq, snvs_cov))

    # close bam
    bam.close()

    outputFile = outputDir + "/" + outputPrefix + ".snv_cov.txt"
    with open(outputFile, 'w') as file:
        file.write("pos\tref\tA\tC\tG\tT\n")
        for key, value in all_snvs_uniq_cov.items():
            line = f"{key[0]}\t{key[1]}\t" + "\t".join(map(str, value))
            file.write(line + "\n")

    pngFile = outputDir + "/" + outputPrefix + ".snv_cov.png"
    plot_frequency(outputFile, pngFile, varRegion)

###############
# program run #
###############

if __name__ == "__main__":
	main(sys.argv[1:])