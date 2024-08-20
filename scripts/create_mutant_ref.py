#######################
#### import module ####
#######################
import os
import sys
import getopt

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

script_name = os.path.basename(__file__)

##########
# ReadMe #
##########
readME = f'''
USAGE:\t{script_name} -i <input file> -s <start position> -e <end position> -o <output file>

Arguments:
\t-i, --input           input fasta file
\t-s, --startpos        target start pos based on alignment
\t-e, --endpos          target end pos based on alignment
\t-o, --output          output fasta file

Optional:
\t-t, --target          target name

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

def get_nnk_codons():
    bases = ['A', 'T', 'G', 'C']
    nnk_codons = []
    for n1 in bases:
        for n2 in bases:
            for k in ['G', 'T']:
                nnk_codons.append(n1 + n2 + k)
    return nnk_codons

def generate_mutations(sequence, start, end, mutated_codons):
    mutated_sequences = []

    start = max(0, start - 1)
    
    if end > len(sequence):
        print(f"Warning: the end position ({end}) exceeds the length of seqeunce ({len(sequence)})!")
        end = len(sequence)

    for i in range(start, end, 3):
        if i + 3 <= len(sequence):
            original_codon = sequence[i:i+3]
            for mutated_codon in mutated_codons:
                if mutated_codon != original_codon:   
                    mutated_sequence = sequence[:i] + mutated_codon + sequence[i+3:]
                    mutated_sequences.append({
                        "mutant": mutated_codon,
                        "sequence": mutated_sequence,
                        "position": i + 1
                    })

    return mutated_sequences

#################
# main function #
#################

def main(argvs):
    target = ''

    try:
        opts, args = getopt.getopt(argvs,
                                   "vhi:s:e:o:t:",
                                   ["version", "help", "input=", "startpos=", "endpos=", "output=", "target="])
        if len(opts) == 0:
            usage_info()
            sys.exit(2)
    except(getopt.GetoptError):
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
        elif opt in ("-o", "--output"):
            outputFile = arg
        elif opt in ("-t", "--target"):
            target = arg
        else:
            assert False, "unhandled option"

    if inputFile == "": sys.exit("Please use the correct arguments, option -i missing")
    if startPos == "": sys.exit("Please use the correct arguments, option -i missing")
    if endPos == "": sys.exit("Please use the correct arguments, option -i missing")
    if outputFile == "": sys.exit("Please use the correct arguments, option -o missing")

    if target == '':
        fileName = os.path.basename(inputFile)
        target = fileName.split('.')[0]

    sequences = []
    for record in SeqIO.parse(inputFile, "fasta"):
        sequences.append({
            "id": record.id, 
            "sequence": str(record.seq)
        }) 

    if len(sequences) != 1:
        sys.exit("Multiple records are found in the input fasta file! Only one is needed!")

    mutated_codons = get_nnk_codons()
    mutated_sequences = generate_mutations(sequences[0]["sequence"], startPos, endPos, mutated_codons)

    records = []
    for idx, seq in enumerate(mutated_sequences):
        header = f"{target}__{seq['position']}__{seq['mutant']}"
        record = SeqRecord(Seq(seq['sequence']), id=header, description="")
        records.append(record)
    
    with open(outputFile, "w") as file:
        SeqIO.write(records, file, "fasta")

###############
# program run #
###############

if __name__ == "__main__":
	main(sys.argv[1:])
