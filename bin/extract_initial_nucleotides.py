#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from Bio.Seq import Seq


parser = argparse.ArgumentParser(description='''Extract initial n bases of each sequence in a (multi)fasta''', add_help=True)
parser.add_argument('-f', '--fast', help="Path to fasta OR fastq file", required=True)
parser.add_argument('-n', '--num', help="Number of initial bases to extract", required=True, type=int)
parser.add_argument('-o', '--out', help="Output file name", required=True)
parser.add_argument('-s', '--side', help="from which end shoud the sequence be extracted 5' or 3', B for both,C to convert into FASTA", required=False, choices=['5', '3', 'C', 'B'], default='B')
parser.add_argument('-v', '--verbose', dest='verbose', action='store_false', help="Verbose output", default=False)
args = parser.parse_args()

fasta = args.fast
num = args.num
out_file = args.out

lost_sequences = 0
with open(fasta, "r") as f, open(out_file, "w") as o:
    if fasta.endswith(".fastq") or fasta.endswith(".fq"):
        file_type = "fastq"
    else:
        file_type = "fasta"
    
    for record in SeqIO.parse(f, file_type):
        seq_id = record.id
        if args.verbose:
            print(f"Processing sequence {seq_id}...")
        sequence = str(record.seq)
        sequence_len = len(sequence)
        current_num = num
        if sequence_len < current_num:
            print(f"Warning: Sequence {seq_id} is shorter than {num} bases. Extracting full length of {sequence_len} bases.")
            current_num = sequence_len
        if args.side == '5':    
            initial_n = sequence[:num]                
        elif args.side == '3':
            initial_n = sequence[-num:]                
        elif args.side == 'C':
            initial_n = sequence
        elif args.side == 'B':
            seq_id = record.id + "_5prime"
            initial_n = sequence[:num] 
            o.write(f">{seq_id}\n{initial_n}\n")
            seq_id = record.id + "_3prime"
            initial_n = sequence[-num:]
            o.write(f">{seq_id}\n{initial_n}\n")
        else:
            print("Error: Invalid side option. Use '5', '3', 'B', or 'C'.")
            continue
            

        if args.side != 'B' and args.side != 'O':
            o.write(f">{seq_id}\n{initial_n}\n")


if args.verbose:
    print("DONE!")


