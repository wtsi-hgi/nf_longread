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
USAGE:\t{script_name} -i <inputfile> -s <start position> -e <end position> -o <output directory>

Arguments:
\t-i, --input           input bam file
\t-s, --startpos        target start pos based on alignment
\t-e, --endpos          target end pos based on alignment
\t-o, --outputdir       output directory

Optional:
\t-t, --target          target name
\t-b, --barcode         the barcode sequence, eg: NNNNATNNNNATNNNN

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

def correct_read(read: str, cigars: tuple):
    newread = ''
    startpos = 0
    endpos = 0
    for cigar in cigars:
        match cigar_chars[cigar[0]]:
            case 'M':
                endpos = startpos + cigar[1]
                newread = newread + read[startpos:endpos]
                startpos = endpos
            case 'I':
                endpos = startpos + cigar[1]
                startpos = endpos
            case 'D':
                newread = newread + '-' * cigar[1]
    return newread

def fix_mismatch(reads: str, indexes_A: list, indexes_T: list, indexes_C: list, indexes_G: list):
    newreads = []
    for read in reads:
        newread = ''
        for index in indexes_A:
            newread = read[:index] + 'A' + read[index + 1:]
            read = newread
        for index in indexes_T:
            newread = read[:index] + 'T' + read[index + 1:]
            read = newread
        for index in indexes_C:
            newread = read[:index] + 'C' + read[index + 1:]
            read = newread
        for index in indexes_G:
            newread = read[:index] + 'G' + read[index + 1:]   
            read = newread         
        newreads.append(newread)
    return newreads

#################
# main function #
#################
cigar_chars = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B']

def main(argvs):
    inputFile = ''
    startPos = ''
    endPos = ''
    outputDir = ''
    target = ''
    barcodeSeq = ''

    try:
        opts, args = getopt.getopt(argvs,
                                   "vhi:s:e:o:t:b:",
                                   ["version", "help", "input=", "startpos=", "endpos=", "outputdir=", "target=", "barcode="])
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
        elif opt in ("-t", "--target"):
            target = arg
        elif opt in ("-b", "--barcode"):
            barcodeSeq = arg
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

    if target == '':
        fileName = os.path.basename(inputFile)
        target = fileName.split('.')[0]

    if barcodeSeq != '':
        indexes_A = []
        indexes_T = []
        indexes_C = []
        indexes_G = []
        for i, char in enumerate(barcodeSeq):
            if char == 'A':
                indexes_A.append(i)
            elif char == 'T':
                indexes_T.append(i)
            elif char == 'C':
                indexes_C.append(i)
            elif char == 'G':
                indexes_G.append(i)

    match_reads = []
    mismatch_reads = []
    samobj = pysam.AlignmentFile(inputFile, "rb")
    for read in samobj.fetch(target, startPos, endPos):
        if read.infer_query_length() == barcodeLen:
            match_reads.append(read.query_sequence)
        else:
            mismatch_reads.append(correct_read(read.query_sequence, read.cigartuples))
    samobj.close()

    all_reads = match_reads + mismatch_reads

    if barcodeSeq != '':
        all_reads_fixed = fix_mismatch(all_reads, indexes_A, indexes_T, indexes_C, indexes_G)
    else:
        all_reads_fixed = all_reads

    all_reads_clean = list(filter(lambda x: '-' not in x, all_reads_fixed))
    all_reads_counts = Counter(all_reads_clean)

    with open(outputDir + "/" + target + ".barcodes.txt", 'w') as file:
        for read, count in all_reads_counts.items():
            file.write(target + "\t" + read + "\t" + str(count) + "\n")

###############
# program run #
###############

if __name__ == "__main__":
	main(sys.argv[1:])