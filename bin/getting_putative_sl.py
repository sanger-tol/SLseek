#!/usr/bin/env python3

import argparse
from Bio import SeqIO 
import ahocorasick 
import math
from Bio.Seq import Seq

# get files from terminal using argparse
parser = argparse.ArgumentParser(description="Retrieve putative SLs - highly expression kmers.")
parser.add_argument("-t", "--transcripts", required=True, help="transcripts extracted ends.")
parser.add_argument("-k", "--kmers_table", required=True, help="tsv table with 1st column as kmer sequence and 2nd as its coverage.")
parser.add_argument("-a", "--abund_thrd", required=False, type =int, help="Abundance thershold for kmers")
parser.add_argument("-m", "--max_distance", required=False, type =int, help="maximum distance between highly abundant kmers", default = 3)
parser.add_argument("-b", "--extra_border", required=False, type = int, help="Extra upstream and downstream nucleotides to main kmer chain", default=0)
parser.add_argument("-y", "--size_limit_max", required=False, type = int, help="If longer than size limit, putative sl will not be recorded", default=40)
parser.add_argument("-x", "--size_limit_min", required=False, type = int, help="If smaller than size limit, putative sl will not be recorded", default=15)
parser.add_argument("-e", "--entropy_lim", required=False, type = float, help="If entropy lower, remove sequence", default=1.40)
parser.add_argument("-n", "--nucleotide_limit", required=False, type = float, help="Percentage of same nucleotide permited in kmers. If above it, kmer not informative, than remove it", default=0.70)
parser.add_argument("-o", "--output", required=True, help="Basename to output.")
args = parser.parse_args()

# Assigning arguments to variables
transcripts = args.transcripts
kmers_table = args.kmers_table
abund_thrd = args.abund_thrd
max_distance = args.max_distance # lim = args.limit
extra_border = args.extra_border
yslim = args.size_limit_max
xslim = args.size_limit_min # slimm
entropy_lim = args.entropy_lim
nucleotide_limit = args.nucleotide_limit
output = args.output

# Removing low abundant kmers and removing uninformative kmers
def is_biased(seq, threshold):
    seq = seq.upper()
    length = len(seq)
    for base in "ATGC":
        if seq.count(base)/length >= threshold:
            return True 
    return False

# Building kmer ahocorasick automaton
kmers_track = set()
kmer_automaton = ahocorasick.Automaton()
for line in open(kmers_table):
    seq, cov = line.strip().split("\t")
    if int(cov) >= abund_thrd:
        if not is_biased(seq, nucleotide_limit):
            seq_revcomp = str(Seq(seq).reverse_complement())
            if seq not in kmers_track:
                kmers_track.add(Seq)
                kmer_automaton.add_word(seq, (cov, seq))
            if seq_revcomp not in kmers_track:
                kmers_track.add(seq_revcomp)
                kmer_automaton.add_word(seq_revcomp, (cov, seq_revcomp))
kmer_automaton.make_automaton()
kmers_track.clear()

# Calculation entropy
def entropy(seq):
    freq = {base: seq.count(base)/len(seq) for base in "ACGT"}
    return -sum(p * math.log2(p) for p in freq.values() if p > 0)

idx_name = 0
output_fasta = f'{output}.fasta'
kmer_len = 0
with open(output_fasta, "w") as out_f:
    # iterate through each FASTA record (transcript end)
    for record in SeqIO.parse(transcripts, "fasta"):
        transcript_seq = str(record.seq)
        # 'pos' will hold lists of end positions for kmers grouped by proximity
        # Example after grouping: pos = [[12,15,17], [55,57], ...]
        pos = []
        cov_sum = []
        group_id = 0

        # iterate over matches from the automaton; automaton.iter yields (end_index, value)
        for end_pos, (cov, kmer) in kmer_automaton.iter(str(record.seq)):
            if kmer_len == 0:
                kmer_len = len(kmer)

            # If no groups yet, start first group with this end_pos
            if pos == []:
                pos.append([end_pos])
                cov_sum.append(int(cov))
            # If current end_pos is within max_distance of the last recorded end in current group,
            # append it to the same group. 
            elif abs(pos[group_id][-1] - end_pos) <= max_distance:
                pos[group_id].append(end_pos)
                cov_sum[group_id] += int(cov)
            else:
                # start a new group for this distant hit
                group_id += 1
                pos.append([end_pos])
                cov_sum.append(int(cov))

        subindex = 0
        for group, coverage in zip(pos, cov_sum):
            start = group[0] - kmer_len + 1 - extra_border
            end = group[-1] + extra_border
            putative_sl = transcript_seq[start:end]

            n_supporting_kmers = len(group)
            avg_cov = coverage/n_supporting_kmers
            if len(putative_sl) <= yslim and len(putative_sl) >= xslim:
                ent = entropy(putative_sl)
                if ent >= entropy_lim:
                    out_f.write(f">putativeSL.{idx_name}_{record.id}.{subindex}_avgCov;{avg_cov}_entropy;{ent}_supportingKmers;{n_supporting_kmers}\n{putative_sl}\n")
                    idx_name += 1
                    subindex += 1
                    


    


