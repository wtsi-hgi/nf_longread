#######################
#### import module ####
#######################
import os
import sys
import getopt
import pysam
from collections import Counter

script_name = os.path.basename(__file__)

##########
# ReadMe #
##########
readME = f'''
USAGE:\t{script_name} -i <input file> -s <start position> -e <end position> -o <output directory>

Arguments:
\t-i, --input           input bam file
\t-s, --startpos        barcode start pos in the reference
\t-e, --endpos          barcode end pos in the reference
\t-o, --outputdir       output directory

Optional:
\t-b, --barcode         the barcode sequence, eg: NNNNATNNNNATNNNN
\t-q, --qualcut         the cutoff of base quality, default: 10
\t-n, --numcut          the cutoff of the number of the low quality bases allowed in the barcode, default: 3
\t-c, --countcut        the cutoff of the number of reads supporting the barcode, default: 5

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

def fix_mismatch(seq: str, indexes_A: list, indexes_T: list, indexes_C: list, indexes_G: list):
    fixed_seq = ''
    for index in indexes_A:
        fixed_seq = seq[:index] + 'A' + seq[index + 1:]
        seq = fixed_seq
    for index in indexes_T:
        fixed_seq = seq[:index] + 'T' + seq[index + 1:]
        seq = fixed_seq
    for index in indexes_C:
        fixed_seq = seq[:index] + 'C' + seq[index + 1:]
        seq = fixed_seq
    for index in indexes_G:
        fixed_seq = seq[:index] + 'G' + seq[index + 1:]  
        seq = fixed_seq         
    return seq

#################
# main function #
#################
cigar_chars = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B']

def main(argvs):
    inputFile = ''
    startPos = ''
    endPos = ''
    outputDir = ''
    barcodeTemplate = ''
    qualCutoff = 10
    numCutoff = 3
    countCutoff = 5

    try:
        opts, args = getopt.getopt(argvs,
                                   "vhi:s:e:o:b:q:n:c:",
                                   ["version", "help", "input=", "startpos=", "endpos=", "outputdir=", "barcode=", "qualcut=", "numcut=", "countcut="])
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
        elif opt in ("-s", "--startpos"):
            startPos = int(arg)
        elif opt in ("-e", "--endpos"):
            endPos = int(arg)
        elif opt in ("-o", "--outputdir"):
            outputDir = arg
        elif opt in ("-b", "--barcode"):
            barcodeTemplate = arg
        elif opt in ("-q", "--qualcut"):
            qualCutoff = int(arg)
        elif opt in ("-n", "--numcut"):
            numCutoff = int(arg)
        elif opt in ("-c", "--countcut"):
            countCutoff = int(arg)
        else:
            assert False, "unhandled option"

    if inputFile == "": 
        sys.exit("Please use the correct arguments, option -i missing")
    if startPos == "": 
        sys.exit("Please use the correct arguments, option -s missing")
    if endPos == "": 
        sys.exit("Please use the correct arguments, option -e missing")
    if outputDir == "": 
        sys.exit("Please use the correct arguments, option -o missing")

    barcodeLen = endPos - startPos + 1

    if barcodeTemplate != '':
        indexes_A = []
        indexes_T = []
        indexes_C = []
        indexes_G = []
        for i, char in enumerate(barcodeTemplate):
            if char == 'A':
                indexes_A.append(i)
            elif char == 'T':
                indexes_T.append(i)
            elif char == 'C':
                indexes_C.append(i)
            elif char == 'G':
                indexes_G.append(i)

    fileName = os.path.basename(inputFile)
    if "__" not in fileName:
        print("Not mutant file, skip")
        sys.exit() 

    target = fileName.split('.')[0]
    var_start = int(target.split('__')[1])
    var_positions = [var_start, var_start + 1, var_start + 2]
    var_positions_str = ','.join(str(pos) for pos in var_positions)

    samobj = pysam.AlignmentFile(inputFile, "rb")
    var_barcode_list = []
    for read in samobj.fetch(target, startPos, endPos):
        if read.is_unmapped:
            continue

        barcode_bases = []
        barcode_qualities = []
        barcode_seq = ''
        
        variant_bases = []
        variant_qualities = []
        variant_seq = ''

        aligned_pairs = read.get_aligned_pairs(matches_only=True)
        for read_pos, ref_pos in aligned_pairs:
            if ref_pos is not None and startPos <= (ref_pos + 1) <= endPos:
                barcode_bases.append(read.query_sequence[read_pos])
                barcode_qualities.append(read.query_qualities[read_pos])
            if (ref_pos + 1) in var_positions:
                variant_bases.append(read.query_sequence[read_pos])
                variant_qualities.append(read.query_qualities[read_pos])
        barcode_seq = ''.join(barcode_bases)
        variant_seq = ''.join(variant_bases)

        # check barcode length and variant length
        if len(barcode_seq) != barcodeLen or len(variant_seq) != len(var_positions):
            continue

        # check barcode quality and variant quality
        barcode_bad_bases = sum(q <= qualCutoff for q in barcode_qualities)
        var_bad_bases = sum(q <= qualCutoff for q in variant_qualities)
        if barcode_bad_bases > numCutoff or var_bad_bases > 0:
            continue

        if barcodeTemplate != '':
            barcode_seq = fix_mismatch(barcode_seq, indexes_A, indexes_T, indexes_C, indexes_G)

        var_barcode_list.append((variant_seq, barcode_seq))

    var_barcode_counts = Counter(var_barcode_list)
    var_barcode_counts_filter = {item: count for item, count in var_barcode_counts.items() if count > countCutoff}

    with open(outputDir + "/" + target + ".barcodes.txt", 'w') as file:
        for item, count in var_barcode_counts_filter.items():
            file.write(f"{var_positions_str}\t{item[0]}\t{item[1]}\t{count}\n")

###############
# program run #
###############

if __name__ == "__main__":
	main(sys.argv[1:])