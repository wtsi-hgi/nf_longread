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
\t-q, --qualcut         the cutoff of base quality, default: 20
\t-n, --numcut          the cutoff of the number of the low quality bases allowed in the barcode, default: 3
\t-c, --countcut        the cutoff of the number of reads supporting the barcode, default: 10

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

def correct_read(read: str, cigars: tuple, quality: str):
    newread = ''
    newquality=''
    startpos = 0
    endpos = 0
    for cigar in cigars:
        match cigar_chars[cigar[0]]:
            case 'M':
                endpos = startpos + cigar[1]
                newread = newread + read[startpos:endpos]
                newquality = newquality + ''.join(chr(q + 33) for q in quality[startpos:endpos])
                startpos = endpos
            case 'I':
                endpos = startpos + cigar[1]
                startpos = endpos
            case 'D':
                newread = newread + '-' * cigar[1]
                newquality = newquality + '-' * cigar[1]
    return newread, newquality

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

def base_quality(reads: str, qualities: str, qual_cutoff: int, num_cutoff: int):
    filtered_reads=[]
    for read, qual in zip(reads, qualities):
        # converting quality score
        quality_score = [(i if isinstance(i, int) else ord(i) - 33) for i in qual]
        # count bases above quality cutoff
        count_base= sum(q <= qual_cutoff for q in quality_score)
        if count_base <= num_cutoff:
            filtered_reads.append(read)
    return filtered_reads

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
    qualCutoff = 20
    numCutoff = 3
    countCutoff =10

    try:
        opts, args = getopt.getopt(argvs,
                                   "vhi:s:e:o:t:b:q:n:c:",
                                   ["version", "help", "input=", "startpos=", "endpos=", "outputdir=", "target=", "barcode=", "qualcut=", "numcut=", "countcut="])
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
    match_qualities = []
    mismatch_qualities = []
    samobj = pysam.AlignmentFile(inputFile, "rb")
    for read in samobj.fetch(target, startPos, endPos):
        if read.infer_query_length() == barcodeLen:
            match_reads.append(read.query_sequence)
            match_qualities.append(read.query_qualities)
        else:
            corrected_read, corrected_quality = correct_read(read.query_sequence, read.cigartuples, read.query_qualities)
            mismatch_reads.append(corrected_read)
            mismatch_qualities.append(corrected_quality)
    samobj.close()

    all_reads = match_reads + mismatch_reads
    all_qualities = match_qualities + mismatch_qualities

    if barcodeSeq != '':
        all_reads_fixed = fix_mismatch(all_reads, indexes_A, indexes_T, indexes_C, indexes_G)
    else:
        all_reads_fixed = all_reads

    # filter reads by barcode quality
    all_reads_filter = base_quality(all_reads_fixed, all_qualities, qualCutoff, numCutoff)
    # clean reads by deletions in barcodes
    all_reads_clean = list(filter(lambda x: '-' not in x, all_reads_filter ))
    # count the unique barcodes
    all_reads_counts = Counter(all_reads_clean)
    all_reads_counts_filter = {read: count for read, count in all_reads_counts.items() if count > countCutoff}

    with open(outputDir + "/" + target + ".barcodes.txt", 'w') as file:
        for read, count in all_reads_counts_filter.items():
            file.write(target + "\t" + str(read) + "\t" + str(count) + "\n")

###############
# program run #
###############

if __name__ == "__main__":
	main(sys.argv[1:])