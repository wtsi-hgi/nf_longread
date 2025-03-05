#######################
#### import module ####
#######################
import os
import sys
import getopt

import re
import pysam
from Bio import SeqIO
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed

script_name = os.path.basename(__file__)

##########
# ReadMe #
##########
readME = f'''
USAGE:\t{script_name} -i <input file> -s <barcode start> -e <barcode end> -o <output dir> -g <gene start> -l <gene length> -f <fasta file>

Arguments:
\t-i, --input           input bam file
\t-s, --startpos        barcode start position in the reference
\t-e, --endpos          barcode end position in the reference
\t-o, --outputdir       output directory
\t-g, --genestart       gene start position in the reference
\t-l, --genelen         gene length
\t-f, --fasta           reference fasta file

Optional:
\t-b, --barcode         the barcode sequence, eg: NNNNATNNNNATNNNN
\t-q, --qualcut         the cutoff of base quality, default: 10
\t-n, --numcut          the cutoff of the number of the low quality bases allowed in the barcode, default: 3
\t-c, --countcut        the cutoff of the number of reads supporting the barcode, default: 5
\t-t, --thread          the number of threads, default: 4

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

# Fix the mismatch in the barcode sequence by template sequence
# Parameters:
# -- seq (str): barcode sequence
# -- indexes_A (list): the list of all the indexes of As in the template sequence
# -- indexes_T (list): the list of all the indexes of Ts in the template sequence
# -- indexes_C (list): the list of all the indexes of Cs in the template sequence
# -- indexes_G (list): the list of all the indexes of Gs in the template sequence
# Returns:
# -- the fixed barcode sequence
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

# Parse the CS tag for mutations occurring after the given start_pos.
# Parameters:
# -- read_obj (pysam.AlignedSegment): The read object in pysam
# -- gene_start (int): Only process positions after this value.
# -- gene_end (int): Only process positions before this value.
# -- ref_seq (str): The reference sequence
# -- qual_cut (int): The quality cutoff for base quality
# Returns:
# -- tuple: Three lists for substitutions, insertions, and deletions.
def parse_cs(read_obj: pysam.AlignedSegment, gene_start: int, gene_end: int, ref_seq: str, qual_cut: int):
    current_position = 0
    query_index = 0
    substitutions = []
    insertions = []
    deletions = []

    # Patterns to identify matches, insertions, deletions, and substitutions in the CS tag
    # Match blocks  --> :52
    # Substitutions --> *ag, *tc
    # Insertions    --> +ATGC
    # Deletions     --> -GC
    cs_tag = read_obj.get_tag("cs")
    pattern = re.compile(r"(:\d+)"
                         r"|(\*\w{1,2})"
                         r"|(\+[ATCGNatcgn]+)"
                         r"|(-[ATCGNatcgn]+)")
    matches = pattern.findall(cs_tag)

    base_qual = read_obj.query_qualities
    for segment in matches:
        segment = next(x for x in segment if x)
        
        # matches
        if segment.startswith(":"): 
            length = int(segment[1:])
            current_position += length
            query_index += length
        # substitutions
        elif segment.startswith("*"):
            if gene_start <= current_position <= gene_end and base_qual[query_index] >= qual_cut:
                ref, alt = segment[1:][0], segment[1:][1]
                substitutions.append(f"{current_position}:{ref}:{alt}:{base_qual[query_index]}")
            current_position += 1
            query_index += 1
        # insertions
        elif segment.startswith("+"):  
            if gene_start <= current_position <= gene_end:
                alt = segment[1:]
                avg_qual = sum(base_qual[query_index:query_index + len(alt)]) // len(alt)
                if avg_qual >= qual_cut:
                    insertions.append(f"{current_position}:INS:{alt}:{avg_qual}")
            query_index += len(segment[1:])
        # deletions
        elif segment.startswith("-"):  
            if gene_start <= current_position <= gene_end:
                ref = segment[1:]
                deletions.append(f"{current_position}:{ref}:DEL:99")
            current_position += len(segment[1:])

    # Create the mutant sequence, be careful about the relative position
    # 1. gene start offset
    # 2. first --> substitution
    # 3. second --> deletion, replace base with "-"
    # 4. third --> insertion
    # 5. remove "-" in the mutant seq
    ref_gene_seq = ref_seq[(gene_start-1):gene_end] 
    mutant_gene_list = list(ref_gene_seq)

    # ref_gene_seq string becomes mutant_gene_list list
    # the first base starts at index 0
    # apply substitutions
    for sub_entry in substitutions:
        pos, ref, alt, _ = sub_entry.split(":")
        gene_pos_index = int(pos) - gene_start
        mutant_gene_list[gene_pos_index] = alt

    # apply deletions
    for del_entry in deletions:
        pos, ref, _, _ = del_entry.split(":")
        del_len = len(ref)
        gene_pos_index = int(pos) - gene_start
        mutant_gene_list[gene_pos_index:(gene_pos_index + del_len)] = "-" * del_len

    # apply insertions
    # in the reverse order, start from the end to avoid indexing issues
    for ins_entry in sorted(insertions, key=lambda x: int(x.split(":")[0]), reverse=True):
        pos, _, alt, _ = ins_entry.split(":")
        pos = int(pos)
        gene_pos_index = int(pos) - gene_start
        mutant_gene_list.insert((gene_pos_index + 1), alt)

    mutant_gene_seq = "".join(mutant_gene_list)

    return substitutions, insertions, deletions, mutant_gene_seq

# multi-threading processing reads
def process_read(read, targetID, startPos, endPos, geneStart, geneEnd, refSeq, qualCutoff, barcodeTemplate, indexes_A, indexes_T, indexes_C, indexes_G, barcodeLen, numCutoff):
    #---------------------------#
    # read alignment validation #
    #---------------------------#
    if read.is_unmapped:
        return(read.query_name, "no barcode", "read unmapped", "F")

    if read.reference_start > geneStart and read.reference_end < geneEnd:
        return(read.query_name, "no barcode", "read clipped", "F")

    # Note: not perfect
    # Reason: reads can be still aligned with mismatches to the target genes
    # Solution: need the full references with all the variants
    # But: too many references, not sure if it is a problem to process
    # Consider: need to examine reads with low mismatches in the gene regions

    #-------------------#
    # barcode detection #
    #-------------------#
    barcode_bases = []
    barcode_qualities = []
    barcode_seq = ''
    
    # get barcode bases and qualities
    aligned_pairs = read.get_aligned_pairs(matches_only = False, with_seq = False)
    for index, (read_pos, ref_pos) in enumerate(aligned_pairs):
            # fix insertions
            if ref_pos is None:
                subindex = index - 1
                while aligned_pairs[subindex][1] is None:
                    subindex = subindex - 1
                if startPos <= aligned_pairs[subindex][1] <= endPos:
                    barcode_bases.append(read.query_sequence[read_pos])
                    barcode_qualities.append(read.query_qualities[read_pos])
            else:
                # skip deletions
                if read_pos is not None and startPos <= (ref_pos + 1) <= endPos:
                    barcode_bases.append(read.query_sequence[read_pos])
                    barcode_qualities.append(read.query_qualities[read_pos])
    barcode_seq = ''.join(barcode_bases)
    
    # check barcode length
    if len(barcode_seq) > (barcodeLen + 1) or len(barcode_seq) < (barcodeLen - 1):
        return(read.query_name, barcode_seq, "barcode length", "F")

    # check barcode quality and variant quality
    barcode_bad_bases = sum(q <= qualCutoff for q in barcode_qualities)
    if barcode_bad_bases >= numCutoff:
        return(read.query_name, barcode_seq, "barcode low quality", "F")

    # fix the mismatches in the barcode if possible
    if barcodeTemplate != '' and len(barcode_seq) == barcodeLen:
        barcode_seq = fix_mismatch(barcode_seq, indexes_A, indexes_T, indexes_C, indexes_G)

    return(read.query_name, barcode_seq, targetID, "P")

# Fix the mismatch in the barcode sequence by template sequence
# Parameters:
# -- codon_str (str): codon strting, 3 bases
# -- pos (int): which one needs to be replaced, 1 or 2 or 3
# -- alt (str): the base to be replaced with
def replace_base(codon_str: str, pos: int, alt: str):
    codon_list = list(codon_str)
    codon_list[pos - 1] = alt
    return(''.join(codon_list))

#################
# main function #
#################
def main(argvs):
    inputFile = ''
    startPos = ''
    endPos = ''
    outputDir = ''
    geneStart = ''
    geneLength = ''
    fastaFile = ''
    barcodeTemplate = ''
    qualCutoff = 10
    numCutoff = 3
    countCutoff = 5
    numThreads = 4

    try:
        opts, args = getopt.getopt(argvs,
                                   "vhi:s:e:o:g:l:f:b:q:n:c:t:",
                                   ["version", "help", 
                                    "input=", "startpos=", "endpos=", "outputdir=", "genestart=", "genelen=", "fasta=", 
                                    "barcode=", "qualcut=", "numcut=", "countcut=", "thread="])
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
        elif opt in ("-g", "--genestart"):
            geneStart = int(arg)
        elif opt in ("-l", "--genelen"):
            geneLength = int(arg)
        elif opt in ("-f", "--fasta"):
            fastaFile = arg
        elif opt in ("-b", "--barcode"):
            barcodeTemplate = arg
        elif opt in ("-q", "--qualcut"):
            qualCutoff = int(arg)
        elif opt in ("-n", "--numcut"):
            numCutoff = int(arg)
        elif opt in ("-c", "--countcut"):
            countCutoff = int(arg)
        elif opt in ("-t", "--thread"):
            numThreads = int(arg)
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
    if geneStart == "": 
        sys.exit("Please use the correct arguments, option -g missing")
    if geneLength == "": 
        sys.exit("Please use the correct arguments, option -l missing")
    if fastaFile == "": 
        sys.exit("Please use the correct arguments, option -f missing")

    #----------------#
    # initialization #
    #----------------#
    fileName = os.path.basename(inputFile)
    targetID = fileName.split('.')[0]

    barcodeLen = endPos - startPos + 1
    geneEnd = geneStart + geneLength - 1

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

    refSeq = ""
    for record in SeqIO.parse(fastaFile, "fasta"):
        if record.id == targetID:
            refSeq = record.seq
    if refSeq == "":
        sys.exit("Cannot find the reference in the fasta file, please check the file name.")

    #-----------------------------------#
    # process each read in the bam file #
    #-----------------------------------#
    samobj = pysam.AlignmentFile(inputFile, "rb")

    output_pass = outputDir + "/" + targetID + ".barcodes_pass.txt"
    output_fail = outputDir + "/" + targetID + ".barcodes_fail.txt"
 
    var_barcode_list = []
    with open(output_pass, 'w') as pass_file, open(output_fail, 'w') as fail_file:
        with ThreadPoolExecutor(max_workers = numThreads) as executor:
            for read in samobj.fetch(targetID, startPos, endPos):
                future = executor.submit(process_read, 
                                         read, targetID, startPos, endPos, geneStart, geneEnd, refSeq, qualCutoff, 
                                         barcodeTemplate, indexes_A, indexes_T, indexes_C, indexes_G, barcodeLen, numCutoff)
                result = future.result()
                if result[-1] == "P":
                    pass_file.write("\t".join(map(str, result[:-1])) + "\n")
                    if result[2] != "NA":
                        var_barcode_list.append((result[1], result[2]))
                else:
                    fail_file.write("\t".join(map(str, result[:-1])) + "\n")

###############
# program run #
###############

if __name__ == "__main__":
	main(sys.argv[1:])